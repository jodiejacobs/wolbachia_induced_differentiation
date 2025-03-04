# Scanpy Clustering and DGET Functional Enrichment Analysis
# Written 02/04/2025
# Run with:
#    python scanpy_plot_projected_umap.py -d bulk_adata.h5ad -o output_folder -a "subclustering" --mem 1024

import scanpy as sc
import pandas as pd
import resource
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import bbknn
import anndata
import requests
import time  # For API rate limiting

# DGET API URL
DGET_API_URL = "https://dget.hms.harvard.edu/api/v1/query"

# Argument parsing
parser = argparse.ArgumentParser(description="Scanpy Bulk to Single-cell Integration with Functional Enrichment")

parser.add_argument("--bulk_adata_path", "-d", type=str, required=True, help="Path to the bulk AnnData file")
parser.add_argument("--output_folder", "-o", type=str, required=True, help="Folder to save the output results")
parser.add_argument("--mem", "-m", type=int, required=False, default=1024, help="Memory limit in GB")
parser.add_argument("--annotation", "-a", type=str, required=False, default="subclustering", help="Color based on annotation type")

args = parser.parse_args()

bulk_adata_path = args.bulk_adata_path
output_folder = args.output_folder
mem = args.mem
annotation = args.annotation

# Create output directory
os.makedirs(output_folder, exist_ok=True)
os.makedirs(f"{output_folder}/plots", exist_ok=True)

# Set working directory and figure directory
os.chdir(output_folder)
sc.settings.figdir = output_folder  

# Set memory limit
mem_limit = mem * 1024 * 1024 * 1024  # Convert GB to bytes
try:
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))
except ValueError as e:
    print(f"Error setting memory limit: {e}")

print(f"Memory limit set to: {mem_limit / (1024 ** 3)} GB")

sc.settings.verbosity = 0  
sc.settings.set_figure_params(dpi=600, frameon=False, facecolor='white', format='pdf', vector_friendly=True)

def preprocess(adata):
    """Preprocess the AnnData object."""
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=6)

def neighbors_rank(adata):
    """Compute neighbors and clustering."""
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.paga(adata)
    sc.tl.umap(adata, init_pos='paga')
    sc.pl.umap(adata, color=["leiden"], size=25, alpha=0.5)

def identify_marker_genes(adata, annotation):
    """Identify marker genes for cell types or clusters."""
    sc.tl.rank_genes_groups(adata, annotation, method='wilcoxon')
    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).iloc[:5]

    # Save marker genes
    marker_genes_file = os.path.join(output_folder, "marker_genes.csv")
    marker_genes.to_csv(marker_genes_file, index=False)

    # Generate marker gene plot
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=annotation,
        n_genes=4,
        values_to_plot="logfoldchanges", 
        cmap='bwr',  
        vmin=-4, vmax=4, min_logfoldchange=3
    )
    plt.savefig('marker_genes_by_tissue.pdf')
    plt.close()

    return marker_genes

### **🔹 DGET Query Function**
def query_dget(genes, species="dmel"):
    """
    Query DGET API to get predicted cell types for a list of genes.
    :param genes: List of gene symbols.
    :param species: Species ('dmel' for Drosophila melanogaster).
    :return: List of predicted cell types.
    """
    payload = {
        "species": species,
        "genes": genes,
        "top_n": 5  # Get the top 5 matching cell types
    }

    try:
        response = requests.post(DGET_API_URL, json=payload, timeout=30)
        response.raise_for_status()  # Raise error for bad status codes (4xx, 5xx)
        data = response.json()
        
        if "predictions" in data:
            cell_types = [entry["cell_type"] for entry in data["predictions"]]
            return cell_types[:1]  # Return the top predicted cell type
        else:
            return ["Unknown"]

    except requests.exceptions.RequestException as e:
        print(f"❌ DGET query failed: {e}")
        return ["Unknown"]

### **🔹 Run DGET for Clusters**
def run_dget(cell_type_gene_lists, species="dmel"):
    """
    Query DGET API for each cluster's marker genes and retrieve predicted cell types.
    :param cell_type_gene_lists: Dictionary where keys are clusters, values are marker genes.
    :param species: Organism ('dmel' for Drosophila).
    :return: Dictionary of assigned cell types.
    """
    assigned_cell_types = {}

    for cluster, gene_list in cell_type_gene_lists.items():
        if len(gene_list) == 0:
            continue  # Skip empty clusters

        print(f"🔍 Querying DGET for cluster: {cluster}")
        
        predicted_types = query_dget(gene_list, species)

        # Assign top predicted cell type
        assigned_cell_types[cluster] = predicted_types[0]

        # Respect API rate limits (avoid 429 errors)
        time.sleep(1)

    return assigned_cell_types

### **🔹 Assign Predicted Cell Types**
def assign_cell_types_to_clusters(adata, assigned_types):
    """Assign cell type labels to clusters based on DGET predictions."""
    adata.obs['assigned_cell_type'] = adata.obs[annotation].map(assigned_types).fillna("Unknown")

    # Save assigned cell types
    cell_type_file = os.path.join(output_folder, "assigned_cell_types_dget.csv")
    adata.obs[['assigned_cell_type']].to_csv(cell_type_file)

    print(f"✅ Assigned cell types saved to: {cell_type_file}")
    return adata

# **🚀 Load and process data**
bulk_adata = sc.read_h5ad(bulk_adata_path)

# Assign cell type annotations
cellid = bulk_adata.obs['Sample'].index
tissue = [s[:-8] for s in cellid]
bulk_adata.obs[annotation] = tissue

# Perform PCA and UMAP projection
neighbors_rank(bulk_adata)

# Identify marker genes
marker_genes = identify_marker_genes(bulk_adata, annotation)

# Convert marker genes into a dictionary
cell_type_gene_lists = {
    cluster: marker_genes[cluster].dropna().tolist() for cluster in marker_genes.columns
}

# 🚀 **Run DGET API for Cell Type Assignment**
assigned_cell_types = run_dget(cell_type_gene_lists)

# 🔹 **Assign cell types to clusters**
combined_adata = assign_cell_types_to_clusters(bulk_adata, assigned_cell_types)

print("🎉 Analysis complete! DGET cell type assignments are saved.")