# %%
# Completed 2/01/2025, run with a memory limit of 1024 GB using sbatch scanpy_bootstrap.sh
# Scanpy Analysis with Single Bootstrapping Iteration Using concat and bbknn functions 
import scanpy as sc
import pandas as pd
import numpy as np
import resource
import argparse
import os
import anndata
import bbknn
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# %%


bulk_adata_path="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad"
ref_adata_path="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/s_fca_biohub_all_wo_blood_10x.h5ad"
output_folder="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/figures/fca_atlas/"
mem=1024
annotation="annotation"
bs="20250304"

# Create output directory if it does not exist
os.makedirs(output_folder, exist_ok=True)

# Set working directory to output folder
os.chdir(output_folder)
sc.settings.figdir = output_folder  # Set Scanpy's figure directory

# Set memory limit
mem_limit = mem * 1024 * 1024 * 1024  # Convert GB to bytes
try:
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))
except ValueError as e:
    print(f"Error setting memory limit: {e}")

print(f"Memory limit set to: {mem_limit / (1024 ** 3)} GB")

sc.settings.verbosity = 0  
sc.settings.set_figure_params(dpi=600, frameon=False, facecolor='white', format='pdf', vector_friendly=True)



# %%
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

# %%
# Load data
bulk_adata = sc.read_h5ad(bulk_adata_path)
ref_adata = sc.read_h5ad(ref_adata_path)

preprocess(ref_adata)

# Filter to shared genes and explicitly make copies
shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
bulk_adata = bulk_adata[:, shared_genes].copy()
ref_adata = ref_adata[:, shared_genes].copy()

# %%
bulk_adata

# %%
ref_adata

# %%
## Color map 
color_dict={
    'JW18DOX':'#87de87',
    'JW18wMel':'#00aa44',
    'S2DOX':'#ffb380',
    'S2wMel':'#d45500'
}


def kNN_classifier(combined_adata, k=20):
    # kNN-based Label Transfer for Bulk Cells
    print("Performing kNN-based label transfer...")

    # Extract reference and bulk indices
    ref_indices = combined_adata.obs['dataset'] == 'Single-cell'
    bulk_indices = combined_adata.obs['dataset'] == 'Bulk'

    # Create k-d tree for reference cells in UMAP space
    tree = cKDTree(combined_adata.obsm['X_umap'][ref_indices])

    # Query the k nearest neighbors for bulk cells
    distances, indices = tree.query(combined_adata.obsm['X_umap'][bulk_indices], k=20)

    # Assign the majority vote of the nearest neighbors as labels for bulk data
    bulk_labels = pd.DataFrame(combined_adata.obs.loc[ref_indices, annotation].values[indices])

    # Fix: Use `mode()` from pandas to correctly assign categorical labels
    assigned_labels = bulk_labels.mode(axis=1)[0]  # Get the most frequent value

    # Assign labels to bulk cells
    combined_adata.obs.loc[bulk_indices, annotation] = assigned_labels.values

    # Extract just the unlabeled dataset (Bulk) and return results
    output_df = pd.DataFrame({
        'cell_type': combined_adata.obs.loc[bulk_indices, annotation]
    }, index=combined_adata.obs.loc[bulk_indices].index)

    return output_df


# %%
adata = bulk_adata.copy()
ref=ref_adata.copy()

# %%
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

# %%
ref = subsample_celltypes(ref, 'tissue')  # Subsample reference data
ref

# %%
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

# %%
adata

# %%

# Concatenate datasets while preserving metadata
combined_adata = anndata.concat(
    [bulk_adata, ref_adata], 
    label='dataset',  
    keys=['Bulk', 'Single-cell'], 
    join="outer",
    merge="same"
)


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
sc.pl.umap(combined_adata, color=['dataset'], save=f'_combine-dataset_{bs}.pdf')
# Visualization, annotation
sc.pl.umap(combined_adata, color=[annotation], save=f'_combined-{annotation}_{bs}.pdf')

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
plt.savefig(f'combined_dataset_samples_and_tissue_{bs}.pdf', dpi=600, bbox_inches="tight")
plt.close() 

# Save the plot
plt.savefig(f'marker_gene_clustering_{bs}.pdf', dpi=600)
plt.close()

# Plot Leiden clustering
sc.pl.umap(combined_adata, color=['leiden'], save=f'leiden_clustering_{bs}.pdf')

bootstrap_result = kNN_classifier(combined_adata)
    


# %%

# Run a single bootstrap iteration
# bootstrap_result = single_bootstrap_iteration(bulk_adata, ref_adata, annotation)

# Save the results
bootstrap_results_file = os.path.join(output_folder, f"bootstrap_{bs}.csv")
bootstrap_result.to_csv(bootstrap_results_file, index=True)

print(f"Bootstrap values saved to: {bootstrap_results_file}")

print("Bootstrap analysis complete!")


