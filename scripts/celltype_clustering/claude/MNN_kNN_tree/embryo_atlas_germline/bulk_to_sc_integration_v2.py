# Final working version of the bulk to single-cell integration pipeline
import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import percentileofscore
import warnings
import logging
import sys
import bbknn
import resource
import gc

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Suppress specific warnings
warnings.filterwarnings("ignore", category=FutureWarning)
sc.settings.verbosity = 1


bulk_path='/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad'
ref_path='/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/combined_germline_sg_trachea.h5ad'
output_dir='/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/MNN_kNN_tree/embryo_atlas_germline'
annotation_key='subtypes'
k_neighbors=None
num_permutations=1000
seed=42
mem = 1024

# Color map to match final figures
color_dict={
    'JW18DOX':'#87de87', # green
    'JW18wMel':'#00aa44',  # dark green
    'S2DOX':'#ffb380', # orange
    'S2wMel':'#d45500' # dark orange

}

"""Set up output directory and plotting parameters."""
np.random.seed(seed)

# Create output directories
os.makedirs(output_dir, exist_ok=True)
plots_dir = os.path.join(output_dir, 'plots')
os.makedirs(plots_dir, exist_ok=True)

# Set up log file
log_file = os.path.join(output_dir, 'integration_log.txt')
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

# Set scanpy settings
sc.settings.figdir = plots_dir
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(10, 8), facecolor='white')

# Set memory limit
mem_limit = mem * 1024 * 1024 * 1024  # Convert GB to bytes
try:
    resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))
except ValueError as e:
    print(f"Error setting memory limit: {e}")

print(f"Memory limit set to: {mem_limit / (1024 ** 3)} GB")

sc.settings.verbosity = 0  
sc.settings.set_figure_params(dpi=600, frameon=False, facecolor='white', format='pdf', vector_friendly=True)


def subsample_celltypes(adata, annotation_key):
    """Subsample cell types to ensure even representation."""
    # Get the cell type counts
    celltype_counts = adata.obs[annotation_key].value_counts()

    # Get the minimum cell type count
    min_count = celltype_counts.min()

    # Create a list to store the subsampled AnnData objects
    subsampled_adatas = []

    # Iterate over each cell type
    for celltype in celltype_counts.index:
        # Get the indices for the current cell type
        celltype_indices = adata.obs[annotation_key] == celltype

        # Subsample the current cell type
        subsampled_adata = adata[celltype_indices].copy()
        subsampled_adata = subsampled_adata[np.random.choice(subsampled_adata.shape[0], min_count, replace=False)]

        # Append the subsampled AnnData object to the list
        subsampled_adatas.append(subsampled_adata)

    # Concatenate the subsampled AnnData objects
    subsampled_adata = ad.AnnData.concatenate(*subsampled_adatas, batch_key='subsampled', index_unique='-')

    return subsampled_adata
    
def load_and_validate_data(bulk_path, ref_path):
    """Load and validate input AnnData objects."""
    logger.info("Loading data files...")
    
    # Load bulk data
    try:
        bulk_adata = sc.read_h5ad(bulk_path)
        logger.info(f"Bulk dataset loaded: {bulk_adata.shape} (samples × genes)")
    except Exception as e:
        logger.error(f"Error loading bulk data: {e}")
        raise
    
    # Load reference data
    try:
        ref_adata = sc.read_h5ad(ref_path)
        logger.info(f"Reference dataset loaded: {ref_adata.shape} (cells × genes)")
    except Exception as e:
        logger.error(f"Error loading reference data: {e}")
        raise
    
    # Ensure unique gene names
    bulk_adata.var_names_make_unique()
    ref_adata.var_names_make_unique()
    
    # Check for shared genes
    shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
    if len(shared_genes) == 0:
        logger.error("No shared genes between bulk and reference datasets!")
        raise ValueError("No shared genes between datasets")
    else:
        logger.info(f"Number of shared genes: {len(shared_genes)}")
    
    return bulk_adata, ref_adata


def preprocess_data(bulk_adata, ref_adata, annotation_key):
    """Preprocess AnnData objects and prepare for integration."""
    logger.info("Preprocessing datasets...")
    
    # Make copies to avoid modifying the originals
    bulk = bulk_adata.copy()
    ref = ref_adata.copy()

    sc.pp.log1p(ref)
    sc.pp.normalize_total(ref, target_sum=1e4)
    
    sc.pp.log1p(bulk)
    sc.pp.normalize_total(bulk, target_sum=1e4)
    
    # Add dataset labels for batch correction
    bulk.obs["dataset"] = "bulk"
    ref.obs["dataset"] = "reference"
    
    # Find shared genes
    shared_genes = bulk.var_names.intersection(ref.var_names)
    logger.info(f"Using {len(shared_genes)} shared genes")
    
    # Subset to shared genes
    bulk_subset = bulk[:, shared_genes].copy()
    ref_subset = ref[:, shared_genes].copy()
    
    # Check if annotation key exists in reference data
    if annotation_key not in ref_subset.obs.columns:
        logger.error(f"Annotation key '{annotation_key}' not found in reference data.")
        available_keys = list(ref_subset.obs.columns)
        logger.error(f"Available keys: {available_keys}")
        raise KeyError(f"Annotation key '{annotation_key}' not found in reference data")
    
    return bulk_subset, ref_subset


def integrate_datasets(bulk_adata, ref_adata):
    """Integrate bulk and single-cell datasets."""
    logger.info("Integrating datasets...")
    
    # Concatenate datasets
    combined = ad.concat([ref_adata, bulk_adata], join="outer", merge="first")
    logger.info(f"Combined dataset shape: {combined.shape}")
    
    # Since data is already normalized and log-transformed, we skip those steps
    logger.info("Data is already normalized and log-transformed")
    
    # Ensure no NaN values that could cause issues
    if sparse.issparse(combined.X):
        combined.X = sparse.csr_matrix(np.nan_to_num(combined.X.toarray(), nan=0, posinf=0, neginf=0))
    else:
        combined.X = np.nan_to_num(combined.X, nan=0, posinf=0, neginf=0)
    

    # Apply BBKNN batch correction
    logger.info("Applying BBKNN batch correction...")

    # Perform batch correction with BBKNN
    bbknn.bbknn(combined, batch_key='dataset')

    return combined



def compute_p_value(neighbor_labels, assigned_label, k, num_permutations=1000):
    """
    Compute p-value by shuffling labels and checking how often
    the assigned label appears by chance.
    IMPORTANT: This function is kept for backward compatibility but is no longer used.
    The p-value calculation is now done directly in the kNN_classifier function with
    distance weighting.
    """
    # Convert to numpy array for efficient operations
    neighbor_labels = np.array(neighbor_labels)
    simulated_counts = []
    
    for _ in range(num_permutations):
        shuffled_labels = np.random.permutation(neighbor_labels)
        simulated_counts.append((shuffled_labels == assigned_label).sum() / k)
    
    observed_prob = (neighbor_labels == assigned_label).sum() / k
    p_value = (100 - percentileofscore(simulated_counts, observed_prob)) / 100
    
    return p_value

def determine_optimal_k(ref_adata, annotation_key):
    """Determine the optimal k value for kNN classification."""
    num_ref_cells = ref_adata.shape[0]
    
    # Get the size of the smallest class
    try:
        min_class_size = ref_adata.obs[annotation_key].value_counts().min()
    except:
        logger.warning("Could not compute min class size, using default")
        min_class_size = 100
    
    # Calculate potential k values:
    # 1. Square root of number of reference cells
    # 2. 10% of smallest class size
    k_sqrt = int(np.sqrt(num_ref_cells))
    k_10pct = int(min_class_size * 0.1)
    
    # Use the smaller of the two values, but ensure k is at least 5
    k = max(5, min(k_sqrt, k_10pct))
    
    logger.info(f"Automatically determined k = {k} (sqrt(n) = {k_sqrt}, 10% of smallest class = {k_10pct})")
    
    return k

# def kNN_classifier(combined_adata, ref_label_key, k, num_permutations=1000):
#     """
#     Classify bulk cells based on their k nearest neighbors in the reference dataset.
#     Each bulk sample is processed independently.
#     """
#     logger.info(f"Performing kNN classification with k={k}...")
    
#     # Identify reference and bulk cells
#     ref_indices = combined_adata.obs["dataset"] == "reference"
#     bulk_indices = combined_adata.obs["dataset"] == "bulk"
    
#     # Get indices as arrays
#     ref_idx = np.where(ref_indices)[0]
#     bulk_idx = np.where(bulk_indices)[0]
    
#     logger.info(f"Reference cells: {len(ref_idx)}, Bulk samples: {len(bulk_idx)}")
    
#     # Extract data for classification
#     try:
#         # Handle sparse matrices if needed
#         if sparse.issparse(combined_adata.X):
#             X_ref = combined_adata.X[ref_idx].toarray()
#             X_bulk = combined_adata.X[bulk_idx].toarray()
#         else:
#             X_ref = combined_adata.X[ref_idx]
#             X_bulk = combined_adata.X[bulk_idx]
        
#         # Handle any NaNs or infs
#         X_ref = np.nan_to_num(X_ref, nan=0, posinf=0, neginf=0)
#         X_bulk = np.nan_to_num(X_bulk, nan=0, posinf=0, neginf=0)
        
#         # Get reference cell labels
#         ref_labels = combined_adata.obs.loc[ref_indices, ref_label_key].values
        
#         # Build kd-tree for efficient nearest neighbor search
#         tree = cKDTree(X_ref)
        
#         results = []
        
#         # Process each bulk sample INDEPENDENTLY
#         for i, bulk_idx_i in enumerate(bulk_idx):
#             # Get current bulk sample
#             current_sample = X_bulk[i].reshape(1, -1)
            
#             # Query KD-tree for this specific sample
#             distances, neighbor_idx = tree.query(current_sample, k=k)
            
#             # Convert from 2D to 1D arrays
#             distances = distances[0]
#             neighbor_idx = neighbor_idx[0]
            
#             # Get labels of k nearest neighbors for this bulk sample
#             neighbor_labels = ref_labels[neighbor_idx]
            
#             # Determine the most frequent label (majority vote)
#             unique_labels, counts = np.unique(neighbor_labels, return_counts=True)
#             assigned_label = unique_labels[np.argmax(counts)]
#             max_count = counts[np.argmax(counts)]
            
#             # Calculate confidence score (percentage of neighbors with the assigned label)
#             confidence = max_count / k
            
#             # Compute p-value with permutation test
#             p_value = compute_p_value(neighbor_labels, assigned_label, k, num_permutations)
            
#             # Store results
#             results.append({
#                 "Bulk_Sample": combined_adata.obs.index[bulk_idx_i],
#                 "Predicted_Label": assigned_label,
#                 "Confidence": confidence,
#                 "P_value": p_value,
#                 "Nearest_Neighbors": list(combined_adata.obs.index[ref_idx[neighbor_idx]]),
#                 "Neighbor_Labels": list(neighbor_labels),
#                 "Neighbor_Distances": list(distances)
#             })
        
#         # Create results DataFrame
#         results_df = pd.DataFrame(results)
        
#         # Add predicted labels to the combined dataset
#         for i, idx in enumerate(bulk_idx):
#             combined_adata.obs.loc[combined_adata.obs.index[idx], ref_label_key] = results_df.iloc[i]["Predicted_Label"]
        
#         logger.info("kNN classification completed")
#         return results_df, combined_adata
    
#     except Exception as e:
#         logger.error(f"kNN classification failed: {e}")
#         import traceback
#         logger.error(traceback.format_exc())
#         raise


def kNN_classifier(combined_adata, ref_label_key, k, num_permutations=1000):
    """
    Classify bulk cells based on their k nearest neighbors in the reference dataset.
    Each bulk sample is processed independently and assigned a cell type based on its own nearest neighbors.
    This implementation takes into account distances to neighbors when calculating confidence scores.
    """
    logger.info(f"Performing kNN classification with k={k}...")
    
    # Verify the annotation key exists and print the first few values
    if ref_label_key not in combined_adata.obs.columns:
        logger.error(f"Annotation key '{ref_label_key}' not found in combined data")
        logger.error(f"Available columns: {list(combined_adata.obs.columns)}")
        
        # Try to find a suitable column as fallback
        categorical_cols = [col for col in combined_adata.obs.columns 
                          if combined_adata.obs[col].dtype.name == 'category' 
                          or combined_adata.obs[col].dtype == 'object']
        
        cell_type_cols = [col for col in categorical_cols 
                         if 'cell' in col.lower() or 'type' in col.lower() 
                         or 'cluster' in col.lower()]
        
        if cell_type_cols:
            fallback_key = cell_type_cols[0]
            logger.warning(f"Using '{fallback_key}' as fallback annotation key")
            ref_label_key = fallback_key
        elif categorical_cols:
            fallback_key = categorical_cols[0]
            logger.warning(f"Using '{fallback_key}' as fallback annotation key")
            ref_label_key = fallback_key
        else:
            raise KeyError(f"Annotation key '{ref_label_key}' not found and no alternative available")
    
    # Print some values from the annotation column to verify it's accessible
    ref_mask = combined_adata.obs['dataset'] == 'reference'
    sample_values = combined_adata.obs.loc[ref_mask, ref_label_key].iloc[:5].tolist()
    logger.info(f"Sample values from '{ref_label_key}': {sample_values}")
    
    # Identify reference and bulk cells
    ref_indices = combined_adata.obs["dataset"] == "reference"
    bulk_indices = combined_adata.obs["dataset"] == "bulk"
    
    # Get cell type labels for the reference dataset
    try:
        ref_adata = combined_adata[ref_indices]
        ref_labels = ref_adata.obs[ref_label_key].values
        logger.info(f"Successfully extracted {len(ref_labels)} reference labels")
    except Exception as e:
        logger.error(f"Error extracting reference labels: {e}")
        raise
    
    # Get indices as arrays
    ref_idx = np.where(ref_indices)[0]
    bulk_idx = np.where(bulk_indices)[0]
    
    logger.info(f"Reference cells: {len(ref_idx)}, Bulk samples: {len(bulk_idx)}")
    
    # Results list to store classification results
    results = []
    
    # Extract embeddings for both reference and bulk data
    if 'X_umap' in combined_adata.obsm:
        logger.info("Using UMAP embeddings for kNN search")
        embedding_key = 'X_umap'
    elif 'X_pca' in combined_adata.obsm:
        logger.info("Using PCA embeddings for kNN search")
        embedding_key = 'X_pca'
    else:
        # Fall back to gene expression space
        logger.info("No dimensionality reduction found, will use gene expression space")
        embedding_key = None
    
    # Process each bulk sample independently
    for i, bulk_idx_i in enumerate(bulk_idx):
        logger.info(f"Processing bulk sample {i+1}/{len(bulk_idx)}: {combined_adata.obs.index[bulk_idx_i]}")
        
        # Get embeddings for this bulk sample and all reference cells
        if embedding_key is not None:
            # Using dimensionality reduction embeddings
            bulk_sample_embedding = combined_adata.obsm[embedding_key][bulk_idx_i].reshape(1, -1)
            ref_embeddings = combined_adata.obsm[embedding_key][ref_idx]
        else:
            # Using gene expression space
            if sparse.issparse(combined_adata.X):
                bulk_sample_embedding = combined_adata.X[bulk_idx_i].toarray().reshape(1, -1)
                ref_embeddings = combined_adata.X[ref_idx].toarray()
            else:
                bulk_sample_embedding = combined_adata.X[bulk_idx_i].reshape(1, -1)
                ref_embeddings = combined_adata.X[ref_idx]
        
        # Handle any NaNs or infs
        bulk_sample_embedding = np.nan_to_num(bulk_sample_embedding, nan=0, posinf=0, neginf=0)
        ref_embeddings = np.nan_to_num(ref_embeddings, nan=0, posinf=0, neginf=0)
        
        # Build kd-tree for efficient nearest neighbor search
        tree = cKDTree(ref_embeddings)
        
        # Query KD-tree to find k nearest neighbors for this bulk sample
        distances, neighbor_indices = tree.query(bulk_sample_embedding, k=k)
        
        # Reshape to 1D arrays
        distances = distances.flatten()
        neighbor_indices = neighbor_indices.flatten()
        
        # Compute distance weights (closer neighbors get higher weights)
        # Add small epsilon to avoid division by zero
        epsilon = 1e-10
        distance_weights = 1.0 / (distances + epsilon)
        
        # Normalize weights to sum to 1
        distance_weights = distance_weights / np.sum(distance_weights)
        
        # Get labels of the k nearest neighbors for this bulk sample
        neighbor_labels = ref_labels[neighbor_indices]
        
        # Calculate weighted votes for each unique label
        unique_labels = np.unique(neighbor_labels)
        weighted_votes = {}
        
        for label in unique_labels:
            # Find indices where this label appears
            label_indices = np.where(neighbor_labels == label)[0]
            # Sum the weights for this label
            weighted_votes[label] = np.sum(distance_weights[label_indices])
        
        # Determine the label with the highest weighted vote
        assigned_label = max(weighted_votes, key=weighted_votes.get)
        max_weighted_vote = weighted_votes[assigned_label]
        
        # Calculate confidence score (percentage of weighted votes for the assigned label)
        confidence = max_weighted_vote
        
        # Compute p-value with permutation test (using distance-weighted voting)
        simulated_votes = []
        for _ in range(num_permutations):
            # Shuffle the labels
            shuffled_labels = np.random.permutation(neighbor_labels)
            # Compute weighted votes for shuffled labels
            shuffled_votes = {}
            for label in np.unique(shuffled_labels):
                label_indices = np.where(shuffled_labels == label)[0]
                shuffled_votes[label] = np.sum(distance_weights[label_indices])
            # Get maximum vote
            max_vote = max(shuffled_votes.values()) if shuffled_votes else 0
            simulated_votes.append(max_vote)
        
        # Calculate p-value as the fraction of simulated maximum votes greater than or equal to our observed vote
        p_value = np.mean(np.array(simulated_votes) >= max_weighted_vote)
        
        # Get indices of reference cells corresponding to the nearest neighbors
        neighbor_ref_indices = ref_idx[neighbor_indices]
        
        # Get the actual cell IDs for the nearest neighbors
        neighbor_cell_ids = [combined_adata.obs.index[idx] for idx in neighbor_ref_indices]
        
        # Store results for this bulk sample
        results.append({
            "Bulk_Sample": combined_adata.obs.index[bulk_idx_i],
            "Predicted_Label": assigned_label,
            "Confidence": confidence,
            "P_value": p_value,
            "Nearest_Neighbors": neighbor_cell_ids,
            "Neighbor_Labels": list(neighbor_labels),
            "Neighbor_Distances": list(distances),
            "Distance_Weights": list(distance_weights)
        })
        
        # Assign the predicted label directly to the original dataset
        combined_adata.obs.loc[combined_adata.obs.index[bulk_idx_i], ref_label_key] = assigned_label
        
        # Force garbage collection after processing each sample
        gc.collect()
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    logger.info(f"kNN classification completed for all {len(bulk_idx)} bulk samples")
    
    # Print a summary of the results
    logger.info("Classification summary:")
    for label, count in results_df['Predicted_Label'].value_counts().items():
        logger.info(f"  {label}: {count} bulk samples (mean confidence: {results_df[results_df['Predicted_Label'] == label]['Confidence'].mean():.4f})")
    
    # Print p-value summary
    logger.info("P-value summary:")
    logger.info(f"  Mean p-value: {results_df['P_value'].mean():.4f}")
    logger.info(f"  Median p-value: {results_df['P_value'].median():.4f}")
    logger.info(f"  Min p-value: {results_df['P_value'].min():.4f}")
    logger.info(f"  Max p-value: {results_df['P_value'].max():.4f}")
    
    return results_df, combined_adata



def visualize_integration(combined_adata, annotation_key, plots_dir):
    """Generate UMAP visualizations of the integrated data."""
    logger.info("Generating UMAP visualizations...")
    
    # Compute PCA
    sc.pp.pca(combined_adata, svd_solver='arpack')
    
    # Compute neighborhood graph
    sc.pp.neighbors(combined_adata, n_neighbors=15, n_pcs=30)
    
    # Compute UMAP embedding
    sc.tl.umap(combined_adata)
    
    # Save plots
    sc.pl.umap(combined_adata, color='dataset', title='Dataset (bulk vs reference)',
               save='_dataset.pdf')
    
    sc.pl.umap(combined_adata, color=annotation_key, title=f'Cell Types ({annotation_key})',
               save=f'_{annotation_key}.pdf')
    
    # Run leiden clustering
    sc.tl.leiden(combined_adata, resolution=0.8)
    sc.pl.umap(combined_adata, color='leiden', title='Leiden Clusters',
               save='_leiden.pdf')
    
    # Create a custom plot highlighting bulk samples
    # fig, ax = plt.subplots(figsize=(10 , 10))
    fig, ax = plt.subplots(figsize=(14, 10))
    plt.tight_layout(rect=[0, 0, 0.8, 1])  # This makes room for the legend
    # Plot reference cells in gray
    ref_mask = combined_adata.obs['dataset'] == 'reference'
    ax.scatter(
        combined_adata.obsm['X_umap'][ref_mask, 0],
        combined_adata.obsm['X_umap'][ref_mask, 1],
        c='lightgray', s=5, alpha=0.5, label='Reference'
    )
    
    # Plot bulk samples with distinct colors based on their predicted cell type
    bulk_mask = combined_adata.obs['dataset'] == 'bulk'
    bulk_cell_types = combined_adata.obs.loc[bulk_mask, annotation_key].astype('category')
    
    for ct in bulk_cell_types.cat.categories:
        ct_mask = (combined_adata.obs['dataset'] == 'bulk') & (combined_adata.obs[annotation_key] == ct)
        ax.scatter(
            combined_adata.obsm['X_umap'][ct_mask, 0],
            combined_adata.obsm['X_umap'][ct_mask, 1],
            s=100, alpha=0.9, label=f'Bulk - {ct}'
        )
    
    ax.set_title('UMAP - Bulk Samples Highlighted')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'bulk_samples_highlighted.pdf'), bbox_inches='tight')
    plt.close()
    
    logger.info("UMAP visualizations completed")
    
    return combined_adata

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

# Log start of processing
logger.info("Starting Scanpy bulk to single-cell integration pipeline")
logger.info(f"Bulk data: {bulk_path}")
logger.info(f"Reference data: {ref_path}")
logger.info(f"Output directory: {output_dir}")

# try:
# Load and validate data
bulk_adata, ref_adata = load_and_validate_data(bulk_path, ref_path)

# Preprocess data
bulk_processed, ref_processed = preprocess_data(bulk_adata, ref_adata, annotation_key)

# Sample even cell types across the reference
ref_processed=subsample_celltypes(ref_processed, annotation_key)

# Integrate datasets
combined_adata = integrate_datasets(bulk_processed, ref_processed)

# Determine optimal k if not provided
k = k_neighbors
if k is None:
    k = determine_optimal_k(ref_processed, annotation_key)

# Run kNN classification
results_df, annotated_adata = kNN_classifier(
    combined_adata,
    ref_label_key=annotation_key,
    k=k,
    num_permutations=num_permutations
)

# Generate visualizations
annotated_adata = visualize_integration(annotated_adata, annotation_key, plots_dir)

# Extract bulk annotations and save to files
bulk_samples = annotated_adata[annotated_adata.obs['dataset'] == 'bulk']
bulk_annotations = bulk_samples.obs[[annotation_key]]

# Save results
results_df.to_csv(os.path.join(output_dir, 'bulk_classification_results.csv'), index=False)
bulk_annotations.to_csv(os.path.join(output_dir, 'bulk_annotations.csv'))
annotated_adata.write_h5ad(os.path.join(output_dir, 'annotated_data.h5ad'))

logger.info("Integration pipeline completed successfully")
logger.info(f"Results saved to {output_dir}")
# Plot UMAP
logger.info("Plotting UMAP...")

# Plot Annotated UMAP:
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(combined_adata, color=[annotation_key], show=False, ax=ax)

# Get the UMAP coordinates
umap_coords = combined_adata.obsm['X_umap']

# Get your sample labels from the data
labels = combined_adata.obs['Sample']

# Collect handles and labels for the legend
legend_handles = []
legend_labels = []

# Iterate over each point and add a label if it's not NA
for idx, label in enumerate(labels):
    if pd.notna(label):  # Check if the label is not NA
        # Plot the point
        point = ax.plot(umap_coords[idx, 0], umap_coords[idx, 1], 
                color=color_dict[label[:-8]], 
                marker='o',            # '*' for star marker
                markersize=18,          # Increase size for better visibility
                markeredgecolor='white', # White outline
                markeredgewidth=0.5,   # Width of the outline
                alpha=1)[0]  # [0] to get the Line2D object
        
        # Only add to legend if this label hasn't been added yet
        if label[:-8] not in legend_labels:
            legend_handles.append(point)
            legend_labels.append(label[:-8])

# Remove grid
ax.grid(False)

# Create a separate figure for the legend
fig_legend = plt.figure(figsize=(6, 2))
fig_legend.legend(legend_handles, legend_labels, loc='center', ncol=3)
fig_legend.savefig(f'{output_dir}/legend.pdf', bbox_inches='tight', dpi=600)
plt.close(fig_legend)

# Save the main plot without legend
plt.savefig(f'{output_dir}/combined_dataset_samples_and_tissue.pdf', dpi=600, bbox_inches="tight")
plt.show()
plt.close()

def plot_enhanced_umap(combined_adata, annotation_key, plots_dir):
    """Generate enhanced UMAP visualizations that better show differences between groups."""
    # Create a figure with two subplots side by side
    fig, axes = plt.subplots(1, 2, figsize=(20, 20))
    
    # Plot 1: Reference cells with bulk cells highlighted by sample group
    ax = axes[0]
    
    # Plot reference cells in gray first
    ref_mask = combined_adata.obs['dataset'] == 'reference'
    ax.scatter(
        combined_adata.obsm['X_umap'][ref_mask, 0],
        combined_adata.obsm['X_umap'][ref_mask, 1],
        c='lightgray', s=5, alpha=0.3, label='Reference Cells'
    )
    
    # Get sample groups (assuming Sample column format like "JW18DOX221117-1")
    bulk_mask = combined_adata.obs['dataset'] == 'bulk'
    sample_groups = []
    for sample_name in combined_adata.obs.loc[bulk_mask, 'Sample']:
        # Extract the prefix (e.g., "JW18DOX", "JW18wMel", "S2DOX", "S2wMel")
        if pd.notna(sample_name):
            prefix = sample_name.split('221117')[0]  # Remove date and number suffix
            sample_groups.append(prefix)
        else:
            sample_groups.append('Unknown')
    
    # Add Sample group as a new column
    combined_adata.obs.loc[bulk_mask, 'SampleGroup'] = sample_groups
    
    # Plot bulk samples with colors based on their sample group
    for group, color in color_dict.items():
        group_mask = (combined_adata.obs['dataset'] == 'bulk') & (combined_adata.obs['SampleGroup'] == group)
        ax.scatter(
            combined_adata.obsm['X_umap'][group_mask, 0],
            combined_adata.obsm['X_umap'][group_mask, 1],
            c=color, s=100, alpha=0.9, label=f'{group}',
            edgecolors='white', linewidths=0.5
        )
    
    ax.set_title('UMAP - Samples Colored by Experimental Condition')
    ax.legend(loc='upper right')
    ax.grid(False)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    
    # Plot 2: Reference cells colored by cell type, bulk samples as larger points
    ax = axes[1]
    
    # Create a colormap for cell types
    cell_types = combined_adata.obs[annotation_key].cat.categories
    n_cell_types = len(cell_types)
    colors = plt.cm.tab20(np.linspace(0, 1, n_cell_types))
    
    # Plot reference cells colored by cell type
    for i, cell_type in enumerate(cell_types):
        ct_mask = (combined_adata.obs['dataset'] == 'reference') & (combined_adata.obs[annotation_key] == cell_type)
        ax.scatter(
            combined_adata.obsm['X_umap'][ct_mask, 0],
            combined_adata.obsm['X_umap'][ct_mask, 1],
            c=[colors[i]], s=10, alpha=0.6, label=f'{cell_type}'
        )
    
    # Plot bulk samples as larger points
    for i, cell_type in enumerate(cell_types):
        ct_mask = (combined_adata.obs['dataset'] == 'bulk') & (combined_adata.obs[annotation_key] == cell_type)
        ax.scatter(
            combined_adata.obsm['X_umap'][ct_mask, 0],
            combined_adata.obsm['X_umap'][ct_mask, 1],
            c=[colors[i]], s=150, alpha=1.0, edgecolors='black', linewidths=0.5
        )
    
    ax.set_title(f'UMAP - Cell Types ({annotation_key})')
    ax.legend(loc='upper right')
    ax.grid(False)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'enhanced_umap_visualization.pdf'), bbox_inches='tight', dpi=300)
    plt.close()
    
    # Create a visualization to verify KNN classification results
    verify_knn_visualization(combined_adata, annotation_key, plots_dir)
    
    return combined_adata

def verify_knn_visualization(combined_adata, annotation_key, plots_dir):
    """Create a visualization to verify KNN classification results for a few samples."""
    bulk_mask = combined_adata.obs['dataset'] == 'bulk'
    ref_mask = combined_adata.obs['dataset'] == 'reference'
    
    # Get sample groups
    sample_groups = combined_adata.obs.loc[bulk_mask, 'SampleGroup'].unique()
    
    # Select one sample from each group to visualize
    selected_samples = []
    for group in sample_groups:
        group_samples = combined_adata.obs.loc[bulk_mask & (combined_adata.obs['SampleGroup'] == group)].index
        if len(group_samples) > 0:
            selected_samples.append(group_samples[0])
    
    # Create a figure with subplots for each selected sample
    fig, axes = plt.subplots(len(selected_samples), 1, figsize=(12, 6*len(selected_samples)))
    if len(selected_samples) == 1:
        axes = [axes]  # Make axes iterable if only one subplot
    
    # For each selected sample
    for i, sample_id in enumerate(selected_samples):
        ax = axes[i]
        
        # Get sample info
        sample_idx = combined_adata.obs.index.get_loc(sample_id)
        sample_group = combined_adata.obs.loc[sample_id, 'SampleGroup']
        predicted_label = combined_adata.obs.loc[sample_id, annotation_key]
        
        # Plot all reference cells as background
        ax.scatter(
            combined_adata.obsm['X_umap'][ref_mask, 0],
            combined_adata.obsm['X_umap'][ref_mask, 1],
            c='lightgray', s=5, alpha=0.2
        )
        
        # Plot cells of the predicted cell type
        pred_mask = ref_mask & (combined_adata.obs[annotation_key] == predicted_label)
        ax.scatter(
            combined_adata.obsm['X_umap'][pred_mask, 0],
            combined_adata.obsm['X_umap'][pred_mask, 1],
            c='blue', s=20, alpha=0.5, label=f'Reference {predicted_label} cells'
        )
        
        # Plot the sample itself
        ax.scatter(
            combined_adata.obsm['X_umap'][sample_idx, 0],
            combined_adata.obsm['X_umap'][sample_idx, 1],
            c=color_dict[sample_group], s=100, alpha=1.0, 
            edgecolors='black', linewidths=1.0,
            label=f'{sample_id} ({sample_group})'
        )
        
        # Get the K nearest neighbors for this sample from the results file
        # Note: This would need the results DataFrame to be passed as an argument
        # For visualization purposes, we'll just use the 10 nearest reference cells based on UMAP distance
        
        # Calculate distances in UMAP space
        sample_umap = combined_adata.obsm['X_umap'][sample_idx]
        ref_umaps = combined_adata.obsm['X_umap'][ref_mask]
        
        # Calculate Euclidean distances
        distances = np.sqrt(np.sum((ref_umaps - sample_umap)**2, axis=1))
        
        # Get indices of 10 nearest neighbors
        nearest_indices = np.argsort(distances)[:10]
        nearest_indices = np.where(ref_mask)[0][nearest_indices]
        
        # Plot nearest neighbors
        ax.scatter(
            combined_adata.obsm['X_umap'][nearest_indices, 0],
            combined_adata.obsm['X_umap'][nearest_indices, 1],
            c='red', s=80, alpha=0.7, 
            edgecolors='black', linewidths=0.5,
            marker='*', label='Nearest neighbors (UMAP space)'
        )
        
        ax.set_title(f'Sample: {sample_id} - Group: {sample_group} - Predicted: {predicted_label}')
        ax.legend(loc='upper right')
        ax.grid(False)
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'knn_verification.pdf'), bbox_inches='tight', dpi=300)
    plt.close()

# Plot enhanced UMAP visualization
plot_enhanced_umap(combined_adata, annotation_key, plots_dir)

# Verify KNN classification results
verify_knn_visualization(combined_adata, annotation_key, plots_dir)
