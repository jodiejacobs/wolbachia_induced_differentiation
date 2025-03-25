# %%
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

# %%
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


# %%
bulk_path='/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad'
ref_path='/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/embryo_adata_dense.h5ad'
output_dir='/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/MNN_kNN_tree/embryo_atlas'
annotation_key='cell_type'
k_neighbors=20
num_permutations=1000
seed=42

# %%
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

# Define custom color palette for cell types
# custom_palette = sns.color_palette("husl", 100)  # Generate a large color palette

# %%

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
    
    # Apply batch correction using MNN
    logger.info("Applying MNN batch correction...")
    try:
        corrected = sce.pp.mnn_correct(combined, batch_key="dataset", return_only_var_genes=False)
        # MNN returns a tuple with the corrected AnnData as the first element
        corrected_adata = corrected[0]
        logger.info("MNN batch correction completed")
        print(type(corrected_adata))
        return corrected_adata
    except Exception as e:
        logger.error(f"MNN batch correction failed: {e}")
        logger.warning("Continuing without batch correction")
        return combined
    

def compute_p_value(neighbor_labels, assigned_label, k, num_permutations=1000):
    """
    Compute p-value by shuffling labels and checking how often
    the assigned label appears by chance.
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


def kNN_classifier(combined_adata, ref_label_key, k, num_permutations=1000):
    """
    Classify bulk cells based on their k nearest neighbors in the reference dataset.
    """
    logger.info(f"Performing kNN classification with k={k}...")
    
    # Identify reference and bulk cells
    ref_indices = combined_adata.obs["dataset"] == "reference"
    bulk_indices = combined_adata.obs["dataset"] == "bulk"
    
    # Get indices as arrays
    ref_idx = np.where(ref_indices)[0]
    bulk_idx = np.where(bulk_indices)[0]
    
    logger.info(f"Reference cells: {len(ref_idx)}, Bulk samples: {len(bulk_idx)}")
    
    # Extract data for classification
    try:
        # Handle sparse matrices if needed
        if sparse.issparse(combined_adata.X):
            X_ref = combined_adata.X[ref_idx].toarray()
            X_bulk = combined_adata.X[bulk_idx].toarray()
        else:
            X_ref = combined_adata.X[ref_idx]
            X_bulk = combined_adata.X[bulk_idx]
        
        # Handle any NaNs or infs
        X_ref = np.nan_to_num(X_ref, nan=0, posinf=0, neginf=0)
        X_bulk = np.nan_to_num(X_bulk, nan=0, posinf=0, neginf=0)
        
        # Build kd-tree for efficient nearest neighbor search
        tree = cKDTree(X_ref)
        distances, neighbor_idx = tree.query(X_bulk, k=k)
        
        # Get reference cell labels
        ref_labels = combined_adata.obs.loc[ref_indices, ref_label_key].values
        
        results = []
        for i, (dists, neighbors) in enumerate(zip(distances, neighbor_idx)):
            # Get labels of k nearest neighbors
            neighbor_labels = ref_labels[neighbors]
            
            # Determine the most frequent label (majority vote)
            unique_labels, counts = np.unique(neighbor_labels, return_counts=True)
            assigned_label = unique_labels[np.argmax(counts)]
            max_count = counts[np.argmax(counts)]
            
            # Calculate confidence score (percentage of neighbors with the assigned label)
            confidence = max_count / k
            
            # Compute p-value with permutation test
            p_value = compute_p_value(neighbor_labels, assigned_label, k, num_permutations)
            
            # Store results
            results.append({
                "Bulk_Sample": combined_adata.obs.index[bulk_idx[i]],
                "Predicted_Label": assigned_label,
                "Confidence": confidence,
                "P_value": p_value,
                "Nearest_Neighbors": list(combined_adata.obs.index[ref_idx[neighbors]]),
                "Neighbor_Labels": list(neighbor_labels),
                "Neighbor_Distances": list(dists)
            })
        
        # Create results DataFrame
        results_df = pd.DataFrame(results)
        
        # Add predicted labels to the combined dataset
        for i, idx in enumerate(bulk_idx):
            combined_adata.obs.loc[combined_adata.obs.index[idx], ref_label_key] = results_df.iloc[i]["Predicted_Label"]
        
        logger.info("kNN classification completed")
        return results_df, combined_adata
    
    except Exception as e:
        logger.error(f"kNN classification failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        raise


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
    fig, ax = plt.subplots(figsize=(10, 8))
    
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



# %%

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





# %%
# Integrate datasets
combined_adata = integrate_datasets(bulk_processed, ref_processed)

combined_adata=combined_adata[0]

type(combined_adata)

# %%

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

# except Exception as e:
# logger.error(f"Integration pipeline failed: {e}")
# import traceback
# logger.error(traceback.format_exc())
# sys.exit(1)



