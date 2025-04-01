# Scanpy Clustering Plotting

# Written 02/04/2025
# Run with:
    # python /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/04022025_scanpy_plot_projected_umap.py \
    #    -d /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad \
    #    -r /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/myeloid_cho_et_al_2020/allnew20210215_circulation.combined.indep_harmony.h5ad \
    #    -o /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/blood_atlas/20250317 \
    #    -a "subclustering" \
    #    --mem 1024


# Run with:
    # python /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/04022025_scanpy_plot_projected_umap.py \
    #    -d /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad \
    #    -r /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/fca_adata.h5ad \
    #    -o /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/fca_atlas/20250317 \
    #    -a "tissue" \
    #    --mem 2024


# Run with:
    # python /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/04022025_scanpy_plot_projected_umap.py \
    #    -d /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad \
    #    -r /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/embryo_adata_dense.h5ad \
    #    -o /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/embryo_atlas/20250317 \
    #    -a "cell_type" \
    #    --mem 2024

import scanpy as sc
import pandas as pd
import resource
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import bbknn
import anndata

# Argument parsing
parser = argparse.ArgumentParser(description="Scanpy Bulk to Single-cell Integration")

parser.add_argument("--bulk_adata_path", "-d", type=str, required=True,
                    default="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad", 
                    help="Path to the bulk AnnData file")

parser.add_argument("--ref_adata_path", "-r", type=str, required=True, 
                    default="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/myeloid_cho_et_al_2020/allnew20210215_circulation.combined.indep_harmony.h5ad",
                    help="Path to the reference AnnData file")

parser.add_argument("--output_folder", "-o", type=str, required=True, 
                    default="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/blood_atlas",
                    help="Folder to save the output results")

parser.add_argument("--mem", "-m", type=int, required=False, 
                    default=1024, 
                    help="Memory limit in GB")

parser.add_argument("--annotation", "-a", type=str, required=False, 
                    default="subclustering", # 'subclustering' for myloid blood atlas, 'tissue' for fly cell atlas, 'annotation' for embryo cell atlas
                    help="Color based on annotation type")

args = parser.parse_args()

bulk_adata_path = args.bulk_adata_path
ref_adata_path = args.ref_adata_path
output_folder = args.output_folder
mem = args.mem
annotation = args.annotation

# Create output directory if it does not exist
os.makedirs(output_folder, exist_ok=True)

# Set working directory to output folder
os.chdir(output_folder)
sc.settings.figdir = output_folder  # Set Scanpy's figure directory

# Create /plots folder if it does not exist
os.makedirs(f'{output_folder}/plots', exist_ok=True)


# Set memory limit
mem_limit = mem * 1024 * 1024 * 1024  # Convert GB to bytes
try:
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))
except ValueError as e:
    print(f"Error setting memory limit: {e}")

print(f"Memory limit set to: {mem_limit / (1024 ** 3)} GB")

sc.settings.verbosity = 0  
sc.settings.set_figure_params(dpi=600, frameon=False, facecolor='white', format='pdf', vector_friendly=True)

# Color map to match final figures
color_dict={
    'JW18DOX':'#87de87', # green
    'JW18wMel':'#00aa44',  # dark green
    'S2DOX':'#ffb380', # orange
    'S2wMel':'#d45500' # dark orange

}

def preprocess(adata):
    """Preprocess the AnnData object."""
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=6)

def neighbors_rank(adata):
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
    sc.tl.umap(adata, init_pos='paga')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=["leiden"], size = 25, alpha = 0.5)

    #find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

def plot_bulk_pca(adata):
    # Generate PCA of the data
    sc.tl.pca(adata, svd_solver='arpack')

    # # Assuming 'bulk_adata' is already loaded, PCA has been performed, and 'in_vitro' cell types are in bulk_adata.obs
    sc.tl.pca(bulk_adata, svd_solver='arpack')

    # Extract PCA results
    pca_coordinates = bulk_adata.obsm["X_pca"][:, :2]  # Get the first two PCs
    variance_ratio = bulk_adata.uns['pca']['variance_ratio']
    cell_types = bulk_adata.obs['in_vitro']

    # Create a figure and an axes object
    fig, ax = plt.subplots()

    # Plot each cell type in its specified color from the color_dict
    for cell_type, color in color_dict.items():
        indices = cell_types == cell_type
        ax.scatter(pca_coordinates[indices, 0], pca_coordinates[indices, 1], label=cell_type, alpha=0.5, color=color)

    # Adding variance explained to the axes labels
    ax.set_xlabel(f'PC1 ({variance_ratio[0]*100:.2f}% variance)')
    ax.set_ylabel(f'PC2 ({variance_ratio[1]*100:.2f}% variance)')

    # Remove background grid and ticks for a cleaner look
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # Add a legend to help identify the cell types
    ax.legend(title='Cell Type')

    # Save the figure
    plt.savefig('bulk_pca.pdf', dpi=300)    

def project_umap(adata, ref):
    # Add a column to identify the dataset type (bulk or single-cell)
    adata.obs['dataset'] = 'Bulk'
    ref.obs['dataset'] = 'Single-cell'

    # Ensure .raw attributes don't interfere with concatenation
    if adata.raw is not None:
        adata.raw = None
    if ref.raw is not None:
        ref.raw = None

    # Ensure no zero values before log transformation
    adata.X = np.clip(adata.X, 1e-10, None)
    ref.X = np.clip(ref.X, 1e-10, None)

    # # Use anndata.concat() instead of the deprecated concatenate()
    # combined_adata = anndata.concat(
    #     [adata, ref], 
    #     label='dataset',  # Adds a 'dataset' column to obs
    #     keys=['Bulk', 'Single-cell'], 
    #     join="outer"  # Keeps all genes, filling missing values with NaN
    # )

    # Combine datasets 
    # # combined_adata = adata.concatenate(ref, batch_key='dataset')
    combined_adata = adata.concatenate(ref, batch_key='dataset', join='outer')
    
    # Ensure variable metadata is retained
    combined_adata.var = adata.var.copy()

    # Merge .uns attributes if present
    combined_adata.uns.update(adata.uns)
    combined_adata.uns.update(ref.uns)

    # # Perform batch correction with BBKNN
    bbknn.bbknn(combined_adata, batch_key='dataset')

    # Perform UMAP
    neighbors_rank(combined_adata)

    # Plot UMAP, save to file
    sc.pl.umap(combined_adata, color=['dataset'], save='combine-dataset.pdf')
    # Visualization
    sc.pl.umap(combined_adata, color=[annotation], save='combined-tissue.pdf')

    #Annotate UMAP with larger markers
    # Plot the UMAP
    sc.pl.umap(combined_adata, color=[annotation], show=False)

    # Get the UMAP coordinates
    umap_coords = combined_adata.obsm['X_umap']

    # Get your sample labels from the data
    labels = combined_adata.obs['Sample']

    # Iterate over each point and add a label if it's not NA
    for idx, label in enumerate(labels):
        if pd.notna(label):  # Check if the label is not NA

            plt.plot(umap_coords[idx, 0], umap_coords[idx, 1], color=color_dict[label[:-8]], marker='o', markersize=3, alpha=0.5)
        # Remove background grid and ticks for a cleaner look
    plt.grid(False)

    # Adjust layout to fit the legend outside
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leaves space for legend on the right
    plt.savefig('combined_dataset_samples_and_tissue.pdf', dpi=600, bbox_inches="tight")
    plt.close() 

    # Save the plot
    plt.savefig('marker_gene_clustering.pdf', dpi=600)
    plt.close()

    # Plot Leiden clustering
    sc.pl.umap(combined_adata, color=['leiden'], save='leiden_clustering.pdf')

    # Identify marker genes
    identify_marker_genes(combined_adata, annotation)

def identify_marker_genes(adata, annotation):
    sc.tl.rank_genes_groups(adata, annotation, method='wilcoxon') #Find marker genes by tissue instead of by leiden clustering (done earlier)
    marker_genes = adata.uns['rank_genes_groups']['names']

    # Get the top N genes for each cluster
    n_top_genes = 5
    top_genes = pd.DataFrame(marker_genes).iloc[:n_top_genes]

    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=annotation,  # Use tissue labels for grouping instead of 'leiden'
        n_genes=4,
        values_to_plot="logfoldchanges", cmap='bwr', #changed from 'viridis'  
        vmin=-4,
        vmax=4,
        min_logfoldchange=3,
        colorbar_title='log fold change'
    )
    plt.savefig('marker_genes_by_tissue.pdf')
    plt.close()


# Load data
bulk_adata = sc.read_h5ad(bulk_adata_path)
ref_adata = sc.read_h5ad(ref_adata_path)
preprocess(ref_adata)

# Add a column to identify the "celltype" in the bulk data
cellid=bulk_adata.obs['Sample'].index
tissue =  [s[:-8] for s in cellid]
bulk_adata.obs[annotation]=tissue

# Filter to shared genes and explicitly make copies, this modification saves memory with large datasets
shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
bulk_adata = bulk_adata[:, shared_genes].copy()  # <-- Fix: Explicitly copy
ref_adata = ref_adata[:, shared_genes].copy()  # <-- Fix: Explicitly copy

# Plot Bulk Data PCA    
neighbors_rank(bulk_adata)
plot_bulk_pca(bulk_adata)

# Project UMAP
project_umap(bulk_adata, ref_adata)

