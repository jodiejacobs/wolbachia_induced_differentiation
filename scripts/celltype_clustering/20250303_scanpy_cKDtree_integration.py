# Scanpy Clustering Plotting

# Written 02/04/2025
# Run with:
    # python /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/04022025_scanpy_plot_projected_umap.py \
    #    -d /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad \
    #    -r /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/myeloid_cho_et_al_2020/allnew20210215_circulation.combined.indep_harmony.h5ad \
    #    -o /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/blood_atlas \
    #    -a "subclustering" \
    #    --mem 1024

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import resource
import argparse
import os
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
import anndata as ad
import scanpy.external as sce
from scipy.stats import percentileofscore
import bbknn

bulk_adata_path ="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad"
ref_adata_path="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/myeloid_cho_et_al_2020/allnew20210215_circulation.combined.indep_harmony.h5ad"
output_folder="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/blood_atlas"
mem=1024
annotation="subclustering"


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

## Color map 
color_dict={
    'JW18DOX':'#87de87',
    'JW18wMel':'#00aa44',
    'S2DOX':'#ffb380',
    'S2wMel':'#d45500'

}

def subsample_celltypes(adata, annotation):
    """Subsample cell types to ensure even representation."""
    # Get the cell type counts
    celltype_counts = adata.obs[annotation].value_counts()

    # Get the minimum cell type count
    min_count = celltype_counts.min()

    # Create a list to store the subsampled AnnData objects
    subsampled_adatas = []

    # Iterate over each cell type
    for celltype in celltype_counts.index:
        # Get the indices for the current cell type
        celltype_indices = adata.obs[annotation] == celltype

        # Subsample the current cell type
        subsampled_adata = adata[celltype_indices].copy()
        subsampled_adata = subsampled_adata[np.random.choice(subsampled_adata.shape[0], min_count, replace=False)]

        # Append the subsampled AnnData object to the list
        subsampled_adatas.append(subsampled_adata)

    # Concatenate the subsampled AnnData objects
    subsampled_adata = anndata.AnnData.concatenate(*subsampled_adatas, batch_key='subsampled', index_unique='-')

    return subsampled_adata


# Load datasets
print("Loading data...")
bulk_adata = sc.read_h5ad(bulk_adata_path)
ref_adata = sc.read_h5ad(ref_adata_path)

 #Before merging datasets, explicitly remove .raw attributes:
if ref_adata.raw is not None:
    ref_adata.raw = None
if bulk_adata.raw is not None:
    bulk_adata.raw = None

# Ensure gene intersection between bulk and reference data
shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
bulk_adata = bulk_adata[:, shared_genes].copy()
ref_adata = ref_adata[:, shared_genes].copy()

# Subsample reference data
ref_adata = subsample_celltypes(ref_adata, annotation)  

# Add a column to identify the "celltype" in the bulk data
cellid=bulk_adata.obs['Sample'].index
tissue =  [s[:-8] for s in cellid]
bulk_adata.obs[annotation]=tissue

# Preprocess reference dataset
print("Preprocessing reference dataset...")
ref_adata.var_names_make_unique()
ref_adata.obs_names_make_unique()
sc.pp.filter_cells(ref_adata, min_genes=200)
sc.pp.filter_genes(ref_adata, min_cells=6)

# Merge reference and bulk datasets
print("Merging datasets...")
bulk_adata.obs["dataset"] = "Unknown"
ref_adata.obs["dataset"] = "Reference"
combined_adata = ad.concat([ref_adata, bulk_adata], join="outer", merge="first")

combined_adata

# Check again after log transformation
combined_adata.X = np.nan_to_num(combined_adata.X, nan=0, posinf=0, neginf=0)

# # # Apply batch correction
sce.pp.mnn_correct(combined_adata, batch_key="dataset")

# Perform batch correction with BBKNN
# bbknn.bbknn(combined_adata, batch_key='dataset')



# kNN-based classification
def kNN_classifier_with_stats(combined_adata, ref_label_key, k):
    """
    Classifies unknown cells based on their k nearest neighbors in the reference dataset.
    Reports neighbor indices, distances, cell type assignments, and statistical confidence.

    Returns:
    - DataFrame with predicted labels, neighbor information, and z-scores.
    """

    print("Performing kNN-based classification...")

    ref_indices = combined_adata.obs["dataset"] == "Reference"
    unknown_indices = combined_adata.obs["dataset"] == "Unknown"

    # Extract gene expression space
    X_ref = combined_adata[ref_indices].X.toarray() if hasattr(combined_adata.X, "toarray") else combined_adata[ref_indices].X
    X_unknown = combined_adata[unknown_indices].X.toarray() if hasattr(combined_adata.X, "toarray") else combined_adata[unknown_indices].X

    # Construct k-d tree for nearest neighbor searching
    tree = cKDTree(X_ref)
    distances, indices = tree.query(X_unknown, k=k)

    # Retrieve reference cell labels
    ref_labels = combined_adata.obs.loc[ref_indices, ref_label_key].values

    results = []
    for i, (dists, neighbor_indices) in enumerate(zip(distances, indices)):
        neighbor_labels = ref_labels[neighbor_indices]
        assigned_label = pd.Series(neighbor_labels).mode()[0]

        # Compute probability
        cell_type_counts = pd.Series(neighbor_labels).value_counts()
        expected_prob = cell_type_counts / k
        observed_prob = expected_prob.loc[assigned_label]
        z_score = (observed_prob - expected_prob.mean()) / expected_prob.std()

        # Compute p-value with permutation test
        p_value = compute_p_value(neighbor_labels, assigned_label, k)

        results.append({
            "Unknown_Cell": combined_adata.obs.index[unknown_indices][i],
            "Predicted_Label": assigned_label,
            "Neighbor_Cells": list(combined_adata.obs.index[ref_indices][neighbor_indices]),
            "Neighbor_Types": list(neighbor_labels),
            "Distances": list(dists),
            "Z-score": z_score,
            "P-value": p_value
        })

    return pd.DataFrame(results)

# p-value computation via permutation test
def compute_p_value(neighbor_labels, assigned_label, k, num_permutations=1000):
    """
    Compute p-value by shuffling labels and checking how often
    the assigned label appears by chance.
    """
    simulated_counts = []
    
    for _ in range(num_permutations):
        shuffled_labels = np.random.permutation(neighbor_labels)
        simulated_counts.append((shuffled_labels == assigned_label).sum() / k)

    observed_prob = (np.array(neighbor_labels) == assigned_label).sum() / k
    p_value = (100 - percentileofscore(simulated_counts, observed_prob)) / 100

    return p_value


def neighbors_rank(adata,k):
    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=20)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=10, directed=False)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
    sc.tl.umap(adata, init_pos='paga')
    sc.tl.umap(adata)
    # sc.pl.umap(adata, color=["leiden"], size = 25, alpha = 0.5,show=False)

    #find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

def plot_UMAP(combined_adata, annotation,k):
    # Perform UMAP
    neighbors_rank(combined_adata,k)

    # Plot UMAP, save to file
    # sc.pl.umap(combined_adata, color=['dataset'], save='combine-dataset.pdf')
    # Visualization
    # sc.pl.umap(combined_adata, color=[annotation], save='combined-tissue.pdf')

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

            plt.plot(umap_coords[idx, 0], umap_coords[idx, 1], color=color_dict[label[:-8]], marker='*', markersize=8, alpha=0.5)
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

