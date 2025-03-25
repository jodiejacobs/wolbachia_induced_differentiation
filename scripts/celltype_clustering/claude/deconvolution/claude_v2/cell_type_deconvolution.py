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
from scipy.stats import zscore, percentileofscore, rankdata
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.preprocessing import normalize
import random
import warnings
from sklearn.utils import resample
from matplotlib.colors import LinearSegmentedColormap

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
    
    # Validate data
    for name, adata in [("Bulk", bulk_adata), ("Reference", ref_adata)]:
        if adata.shape[0] == 0 or adata.shape[1] == 0:
            logger.error(f"{name} dataset is empty: {adata.shape}")
            raise ValueError(f"{name} dataset is empty")
        
        # Make sure var_names and obs_names are unique
        adata.var_names_make_unique()
        adata.obs_names_make_unique()
    
    # Check for shared genes
    shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names)
    if len(shared_genes) == 0:
        logger.error("No shared genes between bulk and reference datasets!")
        raise ValueError("No shared genes between bulk and reference datasets")
    else:
        logger.info(f"Number of shared genes: {len(shared_genes)}")
        
    return bulk_adata, ref_adata

def scale_bulk_to_match_scrna(bulk_adata, ref_adata, shared_genes):
    """Scale bulk data gene-by-gene to match scRNA-seq expression ranges."""
    logger.info("Scaling bulk expression to match scRNA-seq levels (gene-specific scaling)...")
    
    # Extract data for shared genes
    bulk_subset = bulk_adata[:, shared_genes].copy()
    ref_subset = ref_adata[:, shared_genes].copy()
    
    # Calculate statistics for each gene in reference data
    ref_means = {}
    ref_stds = {}
    
    for i, gene in enumerate(shared_genes):
        if sparse.issparse(ref_subset.X):
            gene_expr = ref_subset[:, gene].X.toarray().flatten()
        else:
            gene_expr = ref_subset[:, gene].X.flatten()
        
        ref_means[gene] = np.mean(gene_expr)
        ref_stds[gene] = np.std(gene_expr) if np.std(gene_expr) > 0 else 1.0
    
    # Scale bulk data gene by gene
    if sparse.issparse(bulk_subset.X):
        bulk_matrix = bulk_subset.X.toarray()
    else:
        bulk_matrix = np.array(bulk_subset.X)
    
    # For each gene, scale the bulk expression
    for i, gene in enumerate(shared_genes):
        try:
            # Get bulk expression for this gene
            bulk_expr = bulk_matrix[:, i]
            
            # Calculate bulk statistics
            bulk_mean = np.mean(bulk_expr)
            bulk_std = np.std(bulk_expr) if np.std(bulk_expr) > 0 else 1.0
            
            # Z-score normalize bulk data
            bulk_expr_normalized = (bulk_expr - bulk_mean) / bulk_std
            
            # Scale to reference distribution
            bulk_matrix[:, i] = bulk_expr_normalized * ref_stds[gene] + ref_means[gene]
        except Exception as e:
            logger.warning(f"Failed to scale gene {gene}: {e}")
    
    # Update bulk AnnData object
    if sparse.issparse(bulk_subset.X):
        bulk_subset.X = sparse.csr_matrix(bulk_matrix)
    else:
        bulk_subset.X = bulk_matrix
    
    logger.info("Bulk expression scaling completed")
    return bulk_subset

def quantile_normalize_datasets(bulk_adata, ref_adata, shared_genes):
    """Apply quantile normalization to make bulk data distribution match reference."""
    logger.info("Performing quantile normalization between datasets...")
    
    # Extract expression matrices for shared genes
    bulk_subset = bulk_adata[:, shared_genes].copy()
    ref_subset = ref_adata[:, shared_genes].copy()
    
    # Get reference gene expression values
    if sparse.issparse(ref_subset.X):
        ref_expr = ref_subset.X.mean(axis=0).A1
    else:
        ref_expr = np.array(ref_subset.X.mean(axis=0))
        if len(ref_expr.shape) > 1:
            ref_expr = ref_expr.flatten()
    
    # For each bulk sample, perform quantile normalization to match reference
    if sparse.issparse(bulk_subset.X):
        bulk_expr = bulk_subset.X.toarray()
    else:
        bulk_expr = bulk_subset.X
    
    # Normalize each bulk sample separately
    for i in range(bulk_expr.shape[0]):
        sample = bulk_expr[i, :]
        
        # Get ranks of bulk sample values
        ranks = rankdata(sample)
        
        # Sort reference values
        sorted_ref = np.sort(ref_expr)
        
        # Create a mapping from ranks to reference values
        rank_to_value = {r: v for r, v in zip(range(1, len(sorted_ref) + 1), sorted_ref)}
        
        # Replace each value with the corresponding reference value based on rank
        normalized = np.array([rank_to_value.get(r, 0) for r in ranks])
        
        # Update the original data
        bulk_expr[i, :] = normalized
    
    # Update the bulk AnnData object
    if sparse.issparse(bulk_subset.X):
        bulk_subset.X = sparse.csr_matrix(bulk_expr)
    else:
        bulk_subset.X = bulk_expr
    
    logger.info("Quantile normalization completed")
    return bulk_subset

def harmonize_datasets(bulk_adata, ref_adata, shared_genes):
    """
    Apply a simplified harmonization approach to align bulk and single-cell datasets.
    """
    logger.info("Applying simplified harmonization between datasets...")
    
    # Extract data for shared genes
    bulk_subset = bulk_adata[:, shared_genes].copy()
    ref_subset = ref_adata[:, shared_genes].copy()
    
    # Normalize both datasets first
    # For reference data - normalize to 10,000 counts per cell
    sc.pp.normalize_total(ref_subset, target_sum=1e4)
    sc.pp.log1p(ref_subset)
    
    # For bulk data - normalize to match the same scale
    sc.pp.normalize_total(bulk_subset, target_sum=1e4)
    sc.pp.log1p(bulk_subset)
    
    # Calculate highly variable genes in reference dataset
    sc.pp.highly_variable_genes(ref_subset, n_top_genes=2000, flavor="seurat")
    hvg_mask = ref_subset.var['highly_variable']
    hvg_genes = shared_genes[hvg_mask.values] if hasattr(hvg_mask, 'values') else shared_genes[hvg_mask]
    
    if len(hvg_genes) < 50:
        logger.warning("Too few highly variable genes found, using all genes")
        hvg_genes = shared_genes
    
    logger.info(f"Using {len(hvg_genes)} highly variable genes for harmonization")
    
    # Subset data to highly variable genes
    bulk_hvg = bulk_subset[:, hvg_genes].copy()
    ref_hvg = ref_subset[:, hvg_genes].copy()
    
    # Scale data for PCA
    sc.pp.scale(bulk_hvg, max_value=10)
    sc.pp.scale(ref_hvg, max_value=10)
    
    # Run PCA on both datasets independently
    sc.pp.pca(bulk_hvg, n_comps=min(50, len(hvg_genes) - 1))
    sc.pp.pca(ref_hvg, n_comps=min(50, len(hvg_genes) - 1))
    
    # Use a simple approach to align the PCA spaces
    # Center the reference PCA space
    ref_mean = np.mean(ref_hvg.obsm['X_pca'], axis=0)
    ref_hvg.obsm['X_pca'] = ref_hvg.obsm['X_pca'] - ref_mean
    
    # Center and scale the bulk PCA space to match the reference variance
    bulk_mean = np.mean(bulk_hvg.obsm['X_pca'], axis=0)
    bulk_centered = bulk_hvg.obsm['X_pca'] - bulk_mean
    
    # Calculate the standard deviation for each PC in both datasets
    ref_std = np.std(ref_hvg.obsm['X_pca'], axis=0)
    bulk_std = np.std(bulk_centered, axis=0)
    
    # Scale bulk PCs to match reference standard deviation
    scaling_factors = ref_std / bulk_std
    scaling_factors[~np.isfinite(scaling_factors)] = 1.0  # Handle any zero divisions
    
    # Apply scaling to align the distributions
    for i in range(bulk_centered.shape[1]):
        bulk_centered[:, i] = bulk_centered[:, i] * scaling_factors[i]
    
    # Update the PCA embeddings
    bulk_hvg.obsm['X_pca'] = bulk_centered
    
    # Transfer the corrected embeddings back to the full dataset
    bulk_subset.obsm['X_pca'] = bulk_hvg.obsm['X_pca']
    ref_subset.obsm['X_pca'] = ref_hvg.obsm['X_pca']
    
    # Use PCA to project back to gene space (simplified approach)
    # This is a simple version of batch correction in the original space
    bulk_subset.X = bulk_hvg.X  # Keep the normalized and scaled expression for HVGs
    
    logger.info("Harmonization completed")
    return bulk_subset, ref_subset

def plot_distribution_comparison(bulk_adata, ref_adata, shared_genes, output_dir, stage="before"):
    """Plot gene expression distributions to compare bulk and reference data."""
    try:
        plt.figure(figsize=(12, 6))
        
        # Get bulk expression
        if sparse.issparse(bulk_adata.X):
            bulk_expr = bulk_adata[:, shared_genes].X.toarray().flatten()
        else:
            bulk_expr = bulk_adata[:, shared_genes].X.flatten()
        
        # Get reference expression
        if sparse.issparse(ref_adata.X):
            ref_expr = ref_adata[:, shared_genes].X.toarray().flatten()
        else:
            ref_expr = ref_adata[:, shared_genes].X.flatten()
        
        # Remove zeros for log scale
        bulk_nonzero = bulk_expr[bulk_expr > 0]
        ref_nonzero = ref_expr[ref_expr > 0]
        
        # Create density plots
        sns.kdeplot(bulk_nonzero, label='Bulk')
        sns.kdeplot(ref_nonzero, label='scRNA-seq')
        
        plt.title(f'Expression Distribution Comparison ({stage} normalization)')
        plt.xlabel('Expression Value')
        plt.ylabel('Density')
        plt.legend()
        
        plt.savefig(os.path.join(output_dir, f'expression_distribution_{stage}.pdf'))
        plt.close()
        
        # Also create a log-scale version
        plt.figure(figsize=(12, 6))
        sns.kdeplot(bulk_nonzero, label='Bulk', log_scale=True)
        sns.kdeplot(ref_nonzero, label='scRNA-seq', log_scale=True)
        
        plt.title(f'Expression Distribution Comparison ({stage} normalization) - Log Scale')
        plt.xlabel('Expression Value (log scale)')
        plt.ylabel('Density')
        plt.legend()
        
        plt.savefig(os.path.join(output_dir, f'expression_distribution_{stage}_log.pdf'))
        plt.close()
        
        # Plot per-gene correlation
        plt.figure(figsize=(10, 10))
        
        # Calculate mean expression per gene for both datasets
        if sparse.issparse(bulk_adata.X):
            bulk_means = bulk_adata[:, shared_genes].X.mean(axis=0).A1
        else:
            bulk_means = bulk_adata[:, shared_genes].X.mean(axis=0)
            if len(bulk_means.shape) > 1:
                bulk_means = bulk_means.flatten()
        
        if sparse.issparse(ref_adata.X):
            ref_means = ref_adata[:, shared_genes].X.mean(axis=0).A1
        else:
            ref_means = ref_adata[:, shared_genes].X.mean(axis=0)
            if len(ref_means.shape) > 1:
                ref_means = ref_means.flatten()
        
        # Add small pseudocount to avoid log(0)
        bulk_means = bulk_means + 1e-5
        ref_means = ref_means + 1e-5
        
        plt.scatter(bulk_means, ref_means, alpha=0.5, s=10)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Bulk Mean Expression (log scale)')
        plt.ylabel('scRNA-seq Mean Expression (log scale)')
        plt.title(f'Gene Expression Correlation ({stage} normalization)')
        
        # Add correlation line
        try:
            m, b = np.polyfit(np.log10(bulk_means), np.log10(ref_means), 1)
            x_line = np.logspace(-5, np.log10(max(bulk_means)), 100)
            y_line = 10**(m * np.log10(x_line) + b)
            plt.plot(x_line, y_line, 'r--')
            
            # Calculate correlation
            corr = np.corrcoef(np.log10(bulk_means), np.log10(ref_means))[0, 1]
            plt.text(0.05, 0.95, f'log-log correlation: {corr:.3f}', transform=plt.gca().transAxes)
        except Exception as e:
            logger.warning(f"Could not calculate correlation line: {e}")
        
        plt.savefig(os.path.join(output_dir, f'gene_correlation_{stage}.pdf'))
        plt.close()
        
        logger.info(f"Distribution comparison plots for {stage} normalization saved")
    except Exception as e:
        logger.warning(f"Failed to create distribution comparison plots: {e}")

def run_pca_comparison(bulk_adata, ref_adata, output_dir, stage="before"):
    """Run PCA on combined data and visualize batch effects."""
    try:
        # Create a concatenated dataset for PCA
        combined = ad.concat([bulk_adata, ref_adata], join='inner')
        combined.obs['dataset'] = ['bulk'] * bulk_adata.shape[0] + ['reference'] * ref_adata.shape[0]
        
        # Run PCA
        sc.pp.pca(combined)
        
        # Plot PCA
        plt.figure(figsize=(10, 8))
        bulk_idx = combined.obs['dataset'] == 'bulk'
        ref_idx = combined.obs['dataset'] == 'reference'
        
        plt.scatter(combined.obsm['X_pca'][ref_idx, 0], combined.obsm['X_pca'][ref_idx, 1], 
                  s=5, alpha=0.5, label='scRNA-seq')
        plt.scatter(combined.obsm['X_pca'][bulk_idx, 0], combined.obsm['X_pca'][bulk_idx, 1], 
                  s=50, alpha=0.8, label='Bulk')
        
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title(f'PCA Visualization ({stage} normalization)')
        plt.legend()
        
        plt.savefig(os.path.join(output_dir, f'pca_visualization_{stage}.pdf'))
        plt.close()
        
        # Calculate batch effect strength
        variance_ratio = combined.uns['pca']['variance_ratio']
        plt.figure(figsize=(10, 6))
        plt.plot(np.cumsum(variance_ratio), marker='o')
        plt.xlabel('Number of PCs')
        plt.ylabel('Cumulative Variance Explained')
        plt.title(f'PCA Variance Explained ({stage} normalization)')
        plt.grid(True, alpha=0.3)
        
        plt.savefig(os.path.join(output_dir, f'pca_variance_{stage}.pdf'))
        plt.close()
        
        logger.info(f"PCA comparison for {stage} normalization saved")
    except Exception as e:
        logger.warning(f"Failed to create PCA comparison plots: {e}")

def identify_cell_type_markers(adata, groupby, n_markers=100):
    """
    Identify marker genes for each cell type using a modified TF-IDF approach.
    
    Args:
        adata: AnnData object
        groupby: Column name in adata.obs for cell type annotations
        n_markers: Number of marker genes to select per cell type
    
    Returns:
        dict: Dictionary mapping cell types to marker genes with weights
    """
    logger.info(f"Identifying marker genes for each cell type using {groupby}...")
    
    # Run standard Scanpy differential expression to get initial markers
    sc.tl.rank_genes_groups(adata, groupby, method='wilcoxon')
    
    # Get the list of all cell types
    cell_types = adata.obs[groupby].cat.categories.tolist()
    logger.info(f"Found {len(cell_types)} cell types")
    
    # Create a "document" for each cell type consisting of gene expression
    # Convert adata to dense format for cell type aggregation if needed
    if sparse.issparse(adata.X):
        adata_dense = adata.X.toarray()
    else:
        adata_dense = adata.X
    
    # Create cell type expression profiles (mean expression per cell type)
    cell_type_profiles = {}
    gene_names = adata.var_names.tolist()
    
    for cell_type in cell_types:
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
    for i, cell_type in enumerate(cell_type_profiles):
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
    plot_marker_heatmap(adata, markers_dict, groupby)
    
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
            if len(mean_expr.shape) > 1:
                mean_expr = mean_expr.flatten()
        
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
            if len(b.shape) > 1:
                b = b.flatten()
        
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


def bootstrap_confidence_intervals(bulk_adata, signature_matrix, shared_genes, n_bootstrap=50):
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
        # Resample genes with replacement
        bootstrap_genes = resample(shared_genes, replace=True, n_samples=len(shared_genes))
        
        # Create bootstrapped signature matrix (only including resampled genes)
        bootstrap_sig = signature_matrix.loc[bootstrap_genes]
        
        # Run deconvolution
        bootstrap_results = deconvolve_samples(bulk_subset[:, bootstrap_genes], bootstrap_sig, bootstrap_genes)
        all_results.append(bootstrap_results)
        
        # Log progress
        if (i + 1) % 10 == 0:
            logger.info(f"Completed {i+1}/{n_bootstrap} bootstrap iterations")
    
    # Calculate confidence intervals (2.5th and 97.5th percentiles)
    stacked_results = np.stack([df.values for df in all_results], axis=0)
    lower_bound = np.percentile(stacked_results, 2.5, axis=0)
    upper_bound = np.percentile(stacked_results, 97.5, axis=0)
    
    # Convert to DataFrames
    lower_df = pd.DataFrame(
        lower_bound, 
        index=original_results.index, 
        columns=original_results.columns
    )
    
    upper_df = pd.DataFrame(
        upper_bound, 
        index=original_results.index, 
        columns=original_results.columns
    )
    
    logger.info("Confidence interval calculation completed")
    
    return original_results, lower_df, upper_df


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
    
    # 1. Create heatmap of cell type proportions
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
    
    # 2. Create stacked bar chart of cell type proportions
    plt.figure(figsize=(14, 10))
    
    # Sort cell types by average proportion
    sorted_cols = deconv_results.mean().sort_values(ascending=False).index.tolist()
    
    # Only include top 15 cell types for readability
    if len(sorted_cols) > 15:
        top_cols = sorted_cols[:14]
        # Group remaining cell types as "Other"
        deconv_results['Other'] = deconv_results[sorted_cols[14:]].sum(axis=1)
        sorted_cols = top_cols + ['Other']
    
    # Plot stacked bars
    deconv_results[sorted_cols].plot(
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
    
    # 3. Create composition plot for each sample with confidence intervals
    if lower_ci is not None and upper_ci is not None:
        for sample in deconv_results.index:
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
            
            plt.barh(
                y_pos,
                props[plot_types].values,
                xerr=[props[plot_types].values - lower[plot_types].values,
                      upper[plot_types].values - props[plot_types].values],
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
    
    # 4. Create hierarchical clustering of samples based on cell type composition
    plt.figure(figsize=(14, 10))
    
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
    
    logger.info("Visualizations created")