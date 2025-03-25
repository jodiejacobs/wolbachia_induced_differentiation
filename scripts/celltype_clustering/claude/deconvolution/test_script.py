# %%

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
from scipy import sparse
from scipy.optimize import nnls
from scipy.stats import zscore, percentileofscore
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.preprocessing import normalize
import random
import warnings
from sklearn.utils import resample
from matplotlib.colors import LinearSegmentedColormap

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
ref_path='/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/fca_subset.h5ad'
output_dir='/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/output'
annotation_key='annotation'
n_markers=100
n_bootstrap=1
seed=42

# %%
"""Set up output directory and plotting parameters."""
np.random.seed(seed)
random.seed(seed)

# Create output directories
os.makedirs(output_dir, exist_ok=True)
plots_dir = os.path.join(output_dir, 'plots')
os.makedirs(plots_dir, exist_ok=True)

# Set up log file
log_file = os.path.join(output_dir, 'deconvolution_log.txt')
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

# Set scanpy settings
sc.settings.figdir = plots_dir
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(10, 8), facecolor='white')

# Define custom color palette for cell types
custom_palette = sns.color_palette("husl", 100)  # Generate a large color palette

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
    
def identify_cell_type_markers(adata, groupby, n_markers=100):
    """
    Identify marker genes for each cell type using a modified TF-IDF approach.
    Handles cell types with only one sample.
    
    Args:
        adata: AnnData object
        groupby: Column name in adata.obs for cell type annotations
        n_markers: Number of marker genes to select per cell type
    
    Returns:
        dict: Dictionary mapping cell types to marker genes with weights
    """
    logger.info(f"Identifying marker genes for each cell type using {groupby}...")
    
    # Get the list of all cell types
    cell_types = adata.obs[groupby].cat.categories.tolist()
    logger.info(f"Found {len(cell_types)} cell types")
    
    # Filter out cell types with only one sample
    valid_cell_types = []
    for cell_type in cell_types:
        mask = adata.obs[groupby] == cell_type
        if np.sum(mask) > 1:
            valid_cell_types.append(cell_type)
        else:
            logger.warning(f"Cell type '{cell_type}' has only one sample and will be excluded from differential testing")
    
    logger.info(f"Using {len(valid_cell_types)} cell types with more than one sample for statistical testing")
    
    # Create a subset with only valid cell types
    if len(valid_cell_types) < len(cell_types):
        valid_mask = adata.obs[groupby].isin(valid_cell_types)
        adata_subset = adata[valid_mask].copy()
    else:
        adata_subset = adata.copy()
    
    # Run standard Scanpy differential expression only if we have valid cell types
    if len(valid_cell_types) > 0:
        sc.tl.rank_genes_groups(adata_subset, groupby, method='wilcoxon')
    
    # Create a "document" for each cell type consisting of gene expression
    # Convert adata to dense format for cell type aggregation if needed
    if sparse.issparse(adata.X):
        adata_dense = adata.X.toarray()
    else:
        adata_dense = adata.X
    
    # Create cell type expression profiles (mean expression per cell type)
    cell_type_profiles = {}
    gene_names = adata.var_names.tolist()
    
    for cell_type in cell_types:  # Include ALL cell types here, not just valid ones
        mask = adata.obs[groupby] == cell_type
        # Skip if no cells for this type
        if not np.any(mask):
            logger.warning(f"No cells found for {cell_type}, skipping")
            continue
            
        # Calculate mean expression for this cell type
        cell_type_profiles[cell_type] = np.mean(adata_dense[mask, :], axis=0)
    
    # Create term-frequency matrix (cell types × genes)
    tf_matrix = np.zeros((len(cell_type_profiles), len(gene_names)))
    for i, cell_type in enumerate(cell_type_profiles):
        tf_matrix[i, :] = cell_type_profiles[cell_type]
    
    # Apply TF-IDF transformation
    tfidf = TfidfTransformer()
    tfidf_matrix = tfidf.fit_transform(tf_matrix)
    
    # Convert to dense if sparse
    if sparse.issparse(tfidf_matrix):
        tfidf_matrix = tfidf_matrix.toarray()
    
    # Extract top markers for each cell type
    markers_dict = {}
    cell_type_list = list(cell_type_profiles.keys())
    for i, cell_type in enumerate(cell_type_list):
        # Get gene scores for this cell type
        gene_scores = tfidf_matrix[i, :]
        
        # Sort genes by TF-IDF score
        sorted_indices = np.argsort(-gene_scores)  # Descending order
        
        # Take top n_markers genes
        top_indices = sorted_indices[:n_markers]
        
        # Store gene names and scores
        markers_dict[cell_type] = {
            'genes': [gene_names[idx] for idx in top_indices],
            'scores': [gene_scores[idx] for idx in top_indices]
        }
    
    logger.info(f"Identified {n_markers} marker genes for each of {len(markers_dict)} cell types")
    
    # Plot heatmap of top 10 marker genes per cell type
    try:
        plot_marker_heatmap(adata, markers_dict, groupby)
    except Exception as e:
        logger.warning(f"Could not create marker heatmap: {e}")
    
    return markers_dict


def plot_marker_heatmap(adata, markers_dict, groupby, n_top=10):
    """
    Create a heatmap of top marker genes per cell type.
    
    Args:
        adata: AnnData object
        markers_dict: Dictionary of marker genes per cell type
        groupby: Column name for cell type annotations
        n_top: Number of top genes to include per cell type
    """
    # Collect top n_top genes per cell type
    all_top_genes = []
    for cell_type in markers_dict:
        top_genes = markers_dict[cell_type]['genes'][:n_top]
        all_top_genes.extend(top_genes)
    
    # Remove duplicates while preserving order
    unique_top_genes = []
    for gene in all_top_genes:
        if gene not in unique_top_genes:
            unique_top_genes.append(gene)
    
    # Create AnnData object with just these genes
    if len(unique_top_genes) > 0:
        adata_markers = adata[:, unique_top_genes].copy()
        
        # Plot heatmap
        sc.pl.heatmap(
            adata_markers, var_names=unique_top_genes, 
            groupby=groupby, 
            standard_scale='var',  # Scale by gene
            cmap='viridis',
            swap_axes=True,
            show_gene_labels=True,
            figsize=(14, 10),
            dendrogram=True,
            save="_top_markers.pdf"
        )
        logger.info("Created marker gene heatmap")
    else:
        logger.warning("No marker genes identified for heatmap")


def create_signature_matrix(ref_adata, markers_dict, annotation_key, shared_genes):
    """
    Create a signature matrix from reference data.
    
    Args:
        ref_adata: Reference AnnData object
        markers_dict: Dictionary of marker genes per cell type
        annotation_key: Column name for cell type annotations
        shared_genes: List of genes shared between bulk and reference
    
    Returns:
        DataFrame: Signature matrix with genes as rows and cell types as columns
    """
    logger.info("Creating cell type signature matrix...")
    
    # Filter reference data to include only shared genes
    ref_subset = ref_adata[:, shared_genes].copy()
    
    # Get expression for each cell type
    cell_types = list(markers_dict.keys())
    
    # Create empty signature matrix
    signature_matrix = pd.DataFrame(index=shared_genes, columns=cell_types)
    
    # Fill signature matrix with average expression values
    for cell_type in cell_types:
        # Get cells of this type
        cells = ref_subset[ref_subset.obs[annotation_key] == cell_type]
        
        if cells.shape[0] == 0:
            logger.warning(f"No cells found for {cell_type}, using zeros")
            signature_matrix[cell_type] = 0
            continue
        
        # Calculate mean expression
        if sparse.issparse(cells.X):
            mean_expr = cells.X.mean(axis=0).A1
        else:
            mean_expr = cells.X.mean(axis=0)
        
        # Add to signature matrix
        signature_matrix[cell_type] = mean_expr
    
    # Apply marker gene weighting
    for cell_type in cell_types:
        if cell_type not in markers_dict:
            continue
            
        # Get marker genes for this cell type
        marker_genes = markers_dict[cell_type]['genes']
        marker_scores = markers_dict[cell_type]['scores']
        
        # Only use marker genes that are in shared genes
        valid_markers = []
        valid_scores = []
        for gene, score in zip(marker_genes, marker_scores):
            if gene in shared_genes:
                valid_markers.append(gene)
                valid_scores.append(score)
        
        # Apply weight to marker genes
        for gene, score in zip(valid_markers, valid_scores):
            # Emphasize this gene for this cell type by multiplying by score
            signature_matrix.at[gene, cell_type] *= (1 + score)
    
    # Normalize signature matrix (each cell type column sums to 1)
    signature_matrix = signature_matrix.apply(lambda x: x / x.sum() if x.sum() > 0 else x, axis=0)
    
    logger.info(f"Created signature matrix with {signature_matrix.shape[0]} genes and {signature_matrix.shape[1]} cell types")
    
    return signature_matrix


def deconvolve_samples(bulk_adata, signature_matrix, shared_genes):
    """
    Deconvolve bulk samples into cell type proportions using signature matrix.
    
    Args:
        bulk_adata: Bulk RNA-seq AnnData
        signature_matrix: Signature matrix DataFrame (genes × cell types)
        shared_genes: List of genes shared between bulk and reference
    
    Returns:
        DataFrame: DataFrame with deconvolution results (samples × cell types)
    """
    logger.info("Deconvolving bulk samples into cell type proportions...")
    
    # Filter bulk data to include only shared genes
    bulk_subset = bulk_adata[:, shared_genes].copy()
    
    # Prepare output DataFrame
    results = pd.DataFrame(index=bulk_subset.obs_names, columns=signature_matrix.columns)
    
    # Get signature matrix as numpy array (genes × cell types)
    S = signature_matrix.values
    
    # For each bulk sample
    for i, sample_id in enumerate(bulk_subset.obs_names):
        # Get expression vector for this sample
        if sparse.issparse(bulk_subset.X):
            b = bulk_subset.X[i].toarray().flatten()
        else:
            b = bulk_subset.X[i]
        
        # Solve non-negative least squares problem: min ||Sx - b||^2, s.t. x >= 0
        try:
            proportions, residual = nnls(S, b)
            
            # Normalize proportions to sum to 1
            if np.sum(proportions) > 0:
                proportions = proportions / np.sum(proportions)
            
            # Store results
            results.loc[sample_id] = proportions
            
            # Log progress for every 10th sample
            if (i + 1) % 10 == 0 or i == 0 or i == len(bulk_subset.obs_names) - 1:
                logger.info(f"Deconvolved {i+1}/{len(bulk_subset.obs_names)} bulk samples")
                
        except Exception as e:
            logger.error(f"Error deconvolving sample {sample_id}: {e}")
            results.loc[sample_id] = np.nan
    
    logger.info("Deconvolution completed")
    
    return results


def bootstrap_confidence_intervals(bulk_adata, signature_matrix, shared_genes, n_bootstrap=500):
    """
    Calculate confidence intervals for deconvolution results using bootstrapping.
    
    Args:
        bulk_adata: Bulk RNA-seq AnnData
        signature_matrix: Signature matrix DataFrame
        shared_genes: List of shared genes
        n_bootstrap: Number of bootstrap iterations
    
    Returns:
        tuple: (Deconvolution results, lower bound, upper bound) DataFrames
    """
    logger.info(f"Calculating confidence intervals using {n_bootstrap} bootstrap iterations...")
    
    # Filter bulk data to include only shared genes
    bulk_subset = bulk_adata[:, shared_genes].copy()
    
    # Initialize results storage
    all_results = []
    
    # Original deconvolution results
    original_results = deconvolve_samples(bulk_subset, signature_matrix, shared_genes)
    all_results.append(original_results)
    
    # Bootstrap iterations
    for i in range(n_bootstrap):
        try:
            # Resample genes with replacement
            bootstrap_genes = resample(shared_genes, replace=True, n_samples=len(shared_genes))
            
            # Remove duplicate genes (ensure unique indices for reindexing)
            bootstrap_genes = list(dict.fromkeys(bootstrap_genes))
            
            # Create bootstrapped signature matrix (only including resampled genes)
            bootstrap_sig = signature_matrix.loc[bootstrap_genes].copy()
            
            # Ensure index is unique
            if not bootstrap_sig.index.is_unique:
                logger.warning(f"Duplicate indices found in bootstrap {i}, using unique genes only")
                bootstrap_sig = bootstrap_sig.loc[~bootstrap_sig.index.duplicated(keep='first')]
            
            # Run deconvolution
            bootstrap_results = deconvolve_samples(bulk_subset[:, bootstrap_genes], bootstrap_sig, bootstrap_genes)
            all_results.append(bootstrap_results)
            
            # Log progress
            if (i + 1) % 50 == 0:
                logger.info(f"Completed {i+1}/{n_bootstrap} bootstrap iterations")
                
        except Exception as e:
            logger.warning(f"Error in bootstrap iteration {i}: {e}")
            continue
    
    # Calculate confidence intervals (2.5th and 97.5th percentiles)
    # Convert all DataFrames to numpy arrays with the same shape
    sample_names = original_results.index
    cell_types = original_results.columns
    
    # Initialize arrays for storing results
    result_arrays = []
    
    for result_df in all_results:
        # Reindex to ensure consistent shape
        try:
            aligned_df = result_df.reindex(index=sample_names, columns=cell_types, fill_value=0)
            result_arrays.append(aligned_df.values)
        except Exception as e:
            logger.warning(f"Error aligning bootstrap result: {e}")
            continue
    
    if len(result_arrays) == 0:
        logger.error("No valid bootstrap results!")
        # Return original results with same bounds
        return original_results, original_results.copy(), original_results.copy()
    
    # Stack arrays
    stacked_results = np.stack(result_arrays, axis=0)
    
    # Calculate percentiles
    lower_bound = np.percentile(stacked_results, 2.5, axis=0)
    upper_bound = np.percentile(stacked_results, 97.5, axis=0)
    
    # Convert to DataFrames
    lower_df = pd.DataFrame(
        lower_bound, 
        index=sample_names, 
        columns=cell_types
    )
    
    upper_df = pd.DataFrame(
        upper_bound, 
        index=sample_names, 
        columns=cell_types
    )
    
    logger.info("Confidence interval calculation completed")
    
    return original_results, lower_df, upper_df


def calculate_significance(deconv_results, n_permutations=1000):
    """
    Calculate statistical significance of cell type proportions.
    
    Args:
        deconv_results: Deconvolution results DataFrame
    
    Returns:
        DataFrame: P-values for cell type proportions
    """
    logger.info(f"Calculating significance using {n_permutations} permutations...")
    
    # Initialize p-value DataFrame
    pvalues = pd.DataFrame(index=deconv_results.index, columns=deconv_results.columns)
    
    # For each sample-cell type combination
    for sample in deconv_results.index:
        # Get observed proportions
        obs_proportions = deconv_results.loc[sample].values
        
        # Perform permutation test
        for _ in range(n_permutations):
            # Shuffle proportions
            shuffled = np.random.permutation(obs_proportions)
            
            # For each cell type, count how often shuffled value >= observed
            for i, cell_type in enumerate(deconv_results.columns):
                if 'permutation_counts' not in pvalues.loc[sample, cell_type]:
                    pvalues.at[sample, cell_type] = {'permutation_counts': 0}
                
                if shuffled[i] >= obs_proportions[i]:
                    pvalues.at[sample, cell_type]['permutation_counts'] += 1
    
    # Calculate final p-values
    for sample in pvalues.index:
        for cell_type in pvalues.columns:
            count = pvalues.at[sample, cell_type]['permutation_counts']
            pvalues.at[sample, cell_type] = count / n_permutations
    
    logger.info("Significance calculation completed")
    
    return pvalues


def create_pseudobulk_validation(ref_adata, annotation_key, n_samples=20, min_cell_types=3, max_cell_types=10):
    """
    Create synthetic "pseudobulk" samples from single-cell data for validation.
    
    Args:
        ref_adata: Reference AnnData object
        annotation_key: Column name for cell type annotations
        n_samples: Number of pseudobulk samples to create
        min_cell_types: Minimum number of cell types per pseudobulk
        max_cell_types: Maximum number of cell types per pseudobulk
    
    Returns:
        tuple: (AnnData with pseudobulk samples, DataFrame with true proportions)
    """
    logger.info(f"Creating {n_samples} pseudobulk samples for validation...")
    
    # Get unique cell types
    cell_types = ref_adata.obs[annotation_key].cat.categories.tolist()
    
    # Initialize storage for pseudobulk samples and true proportions
    pseudobulk_X = []
    true_props = []
    
    # Create pseudobulk samples
    for i in range(n_samples):
        # Randomly select number of cell types to include
        n_types = random.randint(min_cell_types, min(max_cell_types, len(cell_types)))
        
        # Randomly select cell types
        selected_types = random.sample(cell_types, n_types)
        
        # Generate random proportions
        props = np.random.dirichlet(np.ones(n_types))
        
        # Initialize pseudobulk vector
        if sparse.issparse(ref_adata.X):
            pseudobulk = np.zeros(ref_adata.shape[1])
        else:
            pseudobulk = np.zeros_like(ref_adata.X[0])
        
        # Add cells according to proportions
        for j, cell_type in enumerate(selected_types):
            # Get cells of this type
            type_cells = ref_adata[ref_adata.obs[annotation_key] == cell_type]
            
            # Skip if no cells
            if type_cells.shape[0] == 0:
                continue
                
            # Sample cells
            n_cells = max(1, int(props[j] * 100))  # At least 1 cell, scale by 100
            sampled_indices = np.random.choice(type_cells.shape[0], n_cells)
            
            # Add to pseudobulk
            if sparse.issparse(type_cells.X):
                cells_subset = type_cells.X[sampled_indices].toarray()
                pseudobulk += cells_subset.sum(axis=0) * props[j]
            else:
                cells_subset = type_cells.X[sampled_indices]
                pseudobulk += cells_subset.sum(axis=0) * props[j]
        
        # Add to storage
        pseudobulk_X.append(pseudobulk)
        
        # Create true proportions vector for all cell types
        true_prop_vec = np.zeros(len(cell_types))
        for j, cell_type in enumerate(selected_types):
            idx = cell_types.index(cell_type)
            true_prop_vec[idx] = props[j]
        
        true_props.append(true_prop_vec)
    
    # Create AnnData object
    pseudobulk_adata = ad.AnnData(
        X=np.vstack(pseudobulk_X),
        var=ref_adata.var.copy()
    )
    
    # Set sample names
    pseudobulk_adata.obs_names = [f"pseudobulk_{i}" for i in range(n_samples)]
    
    # Create true proportions DataFrame
    true_props_df = pd.DataFrame(
        np.vstack(true_props),
        index=pseudobulk_adata.obs_names,
        columns=cell_types
    )
    
    logger.info("Pseudobulk validation samples created")
    
    return pseudobulk_adata, true_props_df


def validate_deconvolution(ref_adata, annotation_key, signature_matrix, shared_genes, n_samples=20):
    """
    Validate deconvolution approach using synthetic pseudobulk samples.
    """
    logger.info("Validating deconvolution approach...")
    
    # Create pseudobulk samples
    pseudobulk, true_props = create_pseudobulk_validation(
        ref_adata, annotation_key, n_samples=n_samples
    )
    
    # Run deconvolution
    deconv_results = deconvolve_samples(pseudobulk, signature_matrix, shared_genes)
    
    # Calculate metrics
    correlations = []
    rmse_values = []
    
    for sample in pseudobulk.obs_names:
        try:
            true = true_props.loc[sample].values  # Convert to numpy array
            pred = deconv_results.loc[sample].values  # Convert to numpy array
            
            # Check for valid arrays
            if len(true) == 0 or len(pred) == 0:
                logger.warning(f"Empty arrays for sample {sample}, skipping")
                continue
                
            # Ensure we have arrays, not scalars
            true = np.array(true, dtype=float)
            pred = np.array(pred, dtype=float)
            
            # Now we can safely check for NaN values
            if np.isnan(true).any() or np.isnan(pred).any():
                logger.warning(f"NaN values found for sample {sample}, skipping")
                continue
                
            # Calculate correlation only if we have variation in both arrays
            if np.std(true) > 0 and np.std(pred) > 0:
                corr = np.corrcoef(true, pred)[0, 1]
                if not np.isnan(corr):
                    correlations.append(corr)
            else:
                logger.warning(f"No variation in data for sample {sample}, skipping correlation")
            
            # Calculate RMSE
            rmse = np.sqrt(np.mean((true - pred) ** 2))
            if not np.isnan(rmse):
                rmse_values.append(rmse)
                
        except Exception as e:
            logger.warning(f"Error processing sample {sample}: {e}")
            continue
    
    # Overall metrics
    if len(correlations) > 0:
        mean_corr = np.mean(correlations)
    else:
        mean_corr = 0
        logger.warning("No valid correlations calculated")
    
    if len(rmse_values) > 0:
        mean_rmse = np.mean(rmse_values)
    else:
        mean_rmse = 0
        logger.warning("No valid RMSE values calculated")
    
    logger.info(f"Validation results: Mean correlation = {mean_corr:.3f}, Mean RMSE = {mean_rmse:.3f}")
    
    # Create validation plot
    try:
        plt.figure(figsize=(10, 10))
        
        # Stack true and predicted proportions for plotting
        true_flat = []
        pred_flat = []
        
        for sample in pseudobulk.obs_names:
            try:
                for cell_type in true_props.columns:
                    true_val = float(true_props.at[sample, cell_type])
                    pred_val = float(deconv_results.at[sample, cell_type])
                    
                    # Skip NaN values
                    if not (np.isnan(true_val) or np.isnan(pred_val)):
                        true_flat.append(true_val)
                        pred_flat.append(pred_val)
            except:
                continue
        
        if len(true_flat) > 0 and len(pred_flat) > 0:
            plt.scatter(true_flat, pred_flat, alpha=0.6)
            plt.plot([0, 1], [0, 1], 'r--')
            plt.xlabel('True Proportion')
            plt.ylabel('Predicted Proportion')
            plt.title(f'Deconvolution Validation\nPearson r = {mean_corr:.3f}, RMSE = {mean_rmse:.3f}')
            plt.tight_layout()
            plt.savefig(os.path.join(sc.settings.figdir, 'deconvolution_validation.pdf'))
            plt.close()
        else:
            logger.warning("Not enough valid data points to create validation plot")
    except Exception as e:
        logger.warning(f"Error creating validation plot: {e}")
    
    return mean_corr

def plot_deconvolution_results(deconv_results, lower_ci=None, upper_ci=None, palette=None):
    """
    Create visualizations of the deconvolution results.
    
    Args:
        deconv_results: DataFrame with deconvolution results
        lower_ci: Lower confidence interval DataFrame
        upper_ci: Upper confidence interval DataFrame
        palette: Color palette for cell types
    """
    logger.info("Creating deconvolution visualizations...")
    
    # Make sure data is numeric
    try:
        deconv_results = deconv_results.astype(float)
        if lower_ci is not None:
            lower_ci = lower_ci.astype(float)
        if upper_ci is not None:
            upper_ci = upper_ci.astype(float)
    except Exception as e:
        logger.warning(f"Error converting results to float: {e}")
    
    # 1. Create heatmap of cell type proportions
    try:
        plt.figure(figsize=(14, 10))
        
        # Sort columns by average proportion
        sorted_cols = deconv_results.mean().sort_values(ascending=False).index
        
        # Create heatmap
        ax = sns.heatmap(
            deconv_results[sorted_cols],
            cmap="viridis",
            linewidths=0.5,
            vmin=0,
            vmax=deconv_results.values.max(),
            cbar_kws={"label": "Proportion"}
        )
        
        plt.title("Cell Type Proportions in Bulk Samples")
        plt.ylabel("Bulk Samples")
        plt.xlabel("Cell Types")
        plt.tight_layout()
        plt.savefig(os.path.join(sc.settings.figdir, 'deconvolution_heatmap.pdf'))
        plt.close()
        logger.info("Created deconvolution heatmap")
    except Exception as e:
        logger.warning(f"Error creating heatmap: {e}")
    
    # 2. Create stacked bar chart of cell type proportions
    try:
        plt.figure(figsize=(14, 10))
        
        # Sort cell types by average proportion
        sorted_cols = deconv_results.mean().sort_values(ascending=False).index.tolist()
        
        # Only include top 15 cell types for readability
        if len(sorted_cols) > 15:
            top_cols = sorted_cols[:14]
            # Create a copy to avoid SettingWithCopyWarning
            plot_df = deconv_results.copy()
            # Group remaining cell types as "Other"
            plot_df['Other'] = plot_df[sorted_cols[14:]].sum(axis=1)
            sorted_cols = top_cols + ['Other']
        else:
            plot_df = deconv_results
        
        # Plot stacked bars
        plot_df[sorted_cols].plot(
            kind='bar',
            stacked=True,
            figsize=(14, 10),
            colormap='tab20' if palette is None else palette
        )
        
        plt.title("Cell Type Composition of Bulk Samples")
        plt.xlabel("Bulk Samples")
        plt.ylabel("Proportion")
        plt.legend(title="Cell Types", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(sc.settings.figdir, 'deconvolution_stacked_bars.pdf'))
        plt.close()
        logger.info("Created stacked bar chart")
    except Exception as e:
        logger.warning(f"Error creating stacked bar chart: {e}")
    
    # 3. Create composition plot for each sample with confidence intervals
    if lower_ci is not None and upper_ci is not None:
        try:
            for sample in deconv_results.index:
                # Skip if sample doesn't exist in all DataFrames
                if sample not in lower_ci.index or sample not in upper_ci.index:
                    logger.warning(f"Sample {sample} missing from confidence intervals, skipping")
                    continue
                    
                plt.figure(figsize=(14, 8))
                
                # Get proportions and CIs for this sample
                props = deconv_results.loc[sample]
                lower = lower_ci.loc[sample]
                upper = upper_ci.loc[sample]
                
                # Sort by proportion
                sorted_idx = np.argsort(-props.values)
                sorted_types = props.index[sorted_idx]
                
                # Only plot top 15 cell types
                if len(sorted_types) > 15:
                    plot_types = sorted_types[:15]
                else:
                    plot_types = sorted_types
                
                # Plot proportions with error bars
                y_pos = np.arange(len(plot_types))
                
                # Convert to float arrays and ensure valid errors
                props_values = props[plot_types].astype(float).values
                lower_values = lower[plot_types].astype(float).values
                upper_values = upper[plot_types].astype(float).values
                
                # Calculate error bars
                xerr_low = np.maximum(0, props_values - lower_values)  # Can't have negative error bars
                xerr_high = np.maximum(0, upper_values - props_values)
                
                plt.barh(
                    y_pos,
                    props_values,
                    xerr=np.vstack([xerr_low, xerr_high]),
                    capsize=5,
                    alpha=0.7,
                    color='skyblue'
                )
                
                plt.yticks(y_pos, plot_types)
                plt.xlabel('Proportion')
                plt.title(f'Cell Type Composition: {sample}')
                plt.grid(axis='x', linestyle='--', alpha=0.7)
                plt.tight_layout()
                plt.savefig(os.path.join(sc.settings.figdir, f'sample_{sample}_composition.pdf'))
                plt.close()
                
            logger.info("Created sample composition plots")
        except Exception as e:
            logger.warning(f"Error creating sample composition plots: {e}")
    
    # 4. Create hierarchical clustering of samples based on cell type composition
    try:
        # Cluster samples
        g = sns.clustermap(
            deconv_results,
            cmap="viridis",
            standard_scale=None,  # Don't standardize
            figsize=(14, 10),
            linewidths=0.5,
            col_cluster=True,
            row_cluster=True,
            vmin=0,
            vmax=deconv_results.values.max(),
            cbar_kws={"label": "Proportion"}
        )
        
        g.fig.suptitle("Hierarchical Clustering of Samples by Cell Type Composition", 
                     fontsize=16, y=1.02)
        plt.savefig(os.path.join(sc.settings.figdir, 'deconvolution_clustering.pdf'))
        plt.close()
        logger.info("Created hierarchical clustering plot")
    except Exception as e:
        logger.warning(f"Error creating hierarchical clustering: {e}")
    
    logger.info("Visualization creation completed")


# %%
# def main():
"""Main function to run the deconvolution pipeline."""
# Log start of processing
logger.info("Starting cell type deconvolution pipeline")
logger.info(f"Bulk data: {bulk_path}")
logger.info(f"Reference data: {ref_path}")
logger.info(f"Output directory: {output_dir}")

# try:
# Option 1: Use the load_and_validate_data function if it exists
# bulk_adata, ref_adata = load_and_validate_data(bulk_data, ref_path)

# Option 2: Load the files directly
logger.info("Loading data files...")
bulk_adata = sc.read_h5ad(bulk_path)
logger.info(f"Bulk dataset loaded: {bulk_adata.shape} (cells × genes)")

ref_adata = sc.read_h5ad(ref_path)
logger.info(f"Reference dataset loaded: {ref_adata.shape} (cells × genes)")

# Validate data
for name, adata in [("Bulk", bulk_adata), ("Reference", ref_adata)]:
    if adata.shape[0] == 0 or adata.shape[1] == 0:
        raise ValueError(f"{name} dataset is empty: {adata.shape}")
    
    # Make sure var_names and obs_names are unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

# Find shared genes
shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names).tolist()
if len(shared_genes) == 0:
    raise ValueError("No shared genes between bulk and reference datasets!")
logger.info(f"Using {len(shared_genes)} shared genes")

# Identify cell type markers
markers_dict = identify_cell_type_markers(
    ref_adata, 
    annotation_key, 
    n_markers=n_markers
)

# Create signature matrix
signature_matrix = create_signature_matrix(
    ref_adata, 
    markers_dict, 
    annotation_key, 
    shared_genes
)

# Save signature matrix
signature_matrix.to_csv(os.path.join(output_dir, 'signature_matrix.csv'))

# Validate deconvolution approach
validation_score = validate_deconvolution(
    ref_adata,
    annotation_key,
    signature_matrix,
    shared_genes,
    n_samples=20
)


# Run deconvolution with confidence intervals
deconv_results, lower_ci, upper_ci = bootstrap_confidence_intervals(
    bulk_adata,
    signature_matrix,
    shared_genes,
    n_bootstrap=n_bootstrap
)


# %%

# Create visualizations
plot_deconvolution_results(deconv_results, lower_ci, upper_ci, custom_palette)

# Save results
deconv_results.to_csv(os.path.join(output_dir, 'deconvolution_results.csv'))
lower_ci.to_csv(os.path.join(output_dir, 'deconvolution_lower_ci.csv'))
upper_ci.to_csv(os.path.join(output_dir, 'deconvolution_upper_ci.csv'))

# Create summary report
with open(os.path.join(output_dir, 'deconvolution_summary.txt'), 'w') as f:
    f.write("Cell Type Deconvolution Summary\n")
    f.write("===============================\n\n")
    f.write(f"Bulk dataset: {bulk_data}\n")
    f.write(f"Reference dataset: {ref_path}\n")
    f.write(f"Number of bulk samples: {bulk_adata.shape[0]}\n")
    f.write(f"Number of reference cells: {ref_adata.shape[0]}\n")
    f.write(f"Number of shared genes: {len(shared_genes)}\n")
    f.write(f"Number of cell types: {len(signature_matrix.columns)}\n\n")
    f.write(f"Validation correlation score: {validation_score:.3f}\n\n")
    
    f.write("Top cell types by average proportion:\n")
    for cell_type, prop in deconv_results.mean().sort_values(ascending=False).items():
        f.write(f"  {cell_type}: {prop:.3f}\n")

logger.info("Deconvolution pipeline completed successfully")
logger.info(f"Results saved to {output_dir}")

# except Exception as e:
#     logger.error(f"Deconvolution pipeline failed: {e}")
#     import traceback
#     logger.error(traceback.format_exc())
#     return 1

# return 0


