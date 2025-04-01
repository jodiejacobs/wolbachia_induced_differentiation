# Scanpy Clustering and GProfiler Functional Enrichment Analysis
# Written 02/04/2025
# Run with:
#    python scanpy_plot_projected_umap.py -d bulk_adata.h5ad -r ref_adata.h5ad -o output_folder -a "subclustering" --mem 1024

import scanpy as sc
import pandas as pd
import resource
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import bbknn
import anndata
from gprofiler import GProfiler  # Import GProfiler

# Initialize g:Profiler API
gp = GProfiler(return_dataframe=True)

# Argument parsing
parser = argparse.ArgumentParser(description="Scanpy Bulk to Single-cell Integration with Functional Enrichment")

parser.add_argument("--bulk_adata_path", "-d", type=str, required=True,
                    help="Path to the bulk AnnData file")

parser.add_argument("--ref_adata_path", "-r", type=str, required=True, 
                    help="Path to the reference AnnData file")

parser.add_argument("--output_folder", "-o", type=str, required=True, 
                    help="Folder to save the output results")

parser.add_argument("--mem", "-m", type=int, required=False, 
                    default=1024, 
                    help="Memory limit in GB")

parser.add_argument("--annotation", "-a", type=str, required=False, 
                    default="subclustering",  
                    help="Color based on annotation type")

args = parser.parse_args()

bulk_adata_path = args.bulk_adata_path
ref_adata_path = args.ref_adata_path
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

## Color map 
color_dict = {
    'JW18DOX': '#C7E576',
    'JW18wMel': '#241E4E',
    'S2DOX': '#F4E04D',
    'S2wMel': '#587792'
}

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

def project_umap(adata, ref):
    """Integrate bulk and reference datasets, perform batch correction and UMAP projection."""
    adata.obs['dataset'] = 'Bulk'
    ref.obs['dataset'] = 'Single-cell'

    # Ensure no zero values
    adata.X = np.clip(adata.X, 1e-10, None)
    ref.X = np.clip(ref.X, 1e-10, None)

    # Concatenate datasets
    combined_adata = anndata.concat(
        [adata, ref], 
        label='dataset',  
        keys=['Bulk', 'Single-cell'], 
        join="outer"
    )

    # Merge metadata
    combined_adata.uns.update(adata.uns)
    combined_adata.uns.update(ref.uns)

    # Perform batch correction with BBKNN
    bbknn.bbknn(combined_adata, batch_key='dataset')

    # Perform clustering and visualization
    neighbors_rank(combined_adata)

    # Save UMAP plots
    sc.pl.umap(combined_adata, color=['dataset'], save='combine-dataset.pdf')
    sc.pl.umap(combined_adata, color=[annotation], save='combined-tissue.pdf')

    return combined_adata

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

def run_gprofiler(cell_type_gene_lists, species="dmelanogaster"):
    """
    Run functional enrichment analysis on gene lists using g:Profiler.
    :param cell_type_gene_lists: Dictionary where keys are cell types and values are lists of genes.
    :param species: Organism name for g:Profiler (default: *Drosophila melanogaster*).
    :return: DataFrame with enriched pathways.
    """
    enriched_results = {}

    for cell_type, gene_list in cell_type_gene_lists.items():
        if len(gene_list) == 0:
            continue  # Skip empty gene lists
        
        print(f"Running g:Profiler for cell type: {cell_type}")
        gprofiler_results = gp.profile(organism=species, query=gene_list)

        if not gprofiler_results.empty:
            enriched_results[cell_type] = gprofiler_results

    return enriched_results

# Load and preprocess data
bulk_adata = sc.read_h5ad(bulk_adata_path)
ref_adata = sc.read_h5ad(ref_adata_path)
preprocess(ref_adata)

# Add cell type annotations for bulk data
cellid = bulk_adata.obs['Sample'].index
tissue = [s[:-8] for s in cellid]
bulk_adata.obs[annotation] = tissue

# Filter to shared genes
shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
bulk_adata = bulk_adata[:, shared_genes].copy()
ref_adata = ref_adata[:, shared_genes].copy()

# Perform PCA and UMAP projection
neighbors_rank(bulk_adata)
combined_adata = project_umap(bulk_adata, ref_adata)

# Identify marker genes
marker_genes = identify_marker_genes(combined_adata, annotation)

# Convert marker genes into cell type gene lists
cell_type_gene_lists = {}
for cell_type in marker_genes.columns:
    cell_type_gene_lists[cell_type] = marker_genes[cell_type].dropna().tolist()

# Run GProfiler enrichment analysis
gprofiler_results = run_gprofiler(cell_type_gene_lists)

# Save g:Profiler results
for cell_type, df in gprofiler_results.items():
    enrichment_file = os.path.join(output_folder, f"gprofiler_{cell_type}.csv")
    df.to_csv(enrichment_file, index=False)
    print(f"Saved g:Profiler enrichment for {cell_type} to {enrichment_file}")

print("GProfiler functional enrichment analysis complete!")