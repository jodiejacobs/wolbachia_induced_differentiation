#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Cell Type Deconvolution for Bulk RNA-seq Using Single-Cell References
=====================================================================

This script implements a regression-based deconvolution approach to estimate
cell type proportions in bulk RNA-seq data using a single-cell reference atlas.
It uses a modified TF-IDF weighting scheme to identify cell type-specific marker
genes and employs constrained optimization for deconvolution.

Author: Claude
Date: March 2025
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
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

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Cell Type Deconvolution for Bulk RNA-seq")
    
    parser.add_argument("--bulk_path", type=str, required=True,
                        help="Path to the bulk RNA-seq AnnData file (h5ad)")
    
    parser.add_argument("--ref_path", type=str, required=True,
                        help="Path to the reference single-cell AnnData file (h5ad)")
    
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save output files")
    
    parser.add_argument("--annotation_key", type=str, default="annotation",
                        help="Key in the reference AnnData.obs containing cell type annotations")
    
    parser.add_argument("--n_markers", type=int, default=100,
                        help="Number of marker genes to use per cell type")
    
    parser.add_argument("--n_bootstrap", type=int, default=500,
                        help="Number of bootstrap iterations for confidence intervals")
    
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    
    return parser.parse_args()


def setup_environment(args):
    """Set up output directory and plotting parameters."""
    np.random.seed(args.seed)
    random.seed(args.seed)
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    plots_dir = os.path.join(args.output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Set up log file
    log_file = os.path.join(args.output_dir, 'deconvolution_log.txt')
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    # Set scanpy settings
    sc.settings.figdir = plots_dir
    sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(10, 8), facecolor='white')
    
    # Define custom color palette for cell types
    custom_palette = sns.color_palette("husl", 100)  # Generate a large color palette
    
    return plots_dir, custom_palette


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
        # Resample genes with replacement
        bootstrap_genes = resample(shared_genes, replace=True, n_samples=len(shared_genes))
        
        # Create bootstrapped signature matrix (only including resampled genes)
        bootstrap_sig = signature_matrix.loc[bootstrap_genes]
        
        # Run deconvolution
        bootstrap_results = deconvolve_samples(bulk_subset[:, bootstrap_genes], bootstrap_sig, bootstrap_genes)
        all_results.append(bootstrap_results)
        
        # Log progress
        if (i + 1) % 50 == 0:
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
    
    Args:
        ref_adata: Reference AnnData
        annotation_key: Column name for cell type annotations
        signature_matrix: Signature matrix
        shared_genes: List of shared genes
        n_samples: Number of validation samples
    
    Returns:
        float: Overall correlation score
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
        true = true_props.loc[sample]
        pred = deconv_results.loc[sample]
        
        # Pearson correlation
        corr = np.corrcoef(true, pred)[0, 1]
        if not np.isnan(corr):
            correlations.append(corr)
        
        # RMSE
        rmse = np.sqrt(np.mean((true - pred) ** 2))
        rmse_values.append(rmse)
    
    # Overall metrics
    mean_corr = np.mean(correlations)
    mean_rmse = np.mean(rmse_values)
    
    logger.info(f"Validation results: Mean correlation = {mean_corr:.3f}, Mean RMSE = {mean_rmse:.3f}")
    
    # Create validation plot
    plt.figure(figsize=(10, 10))
    
    # Stack true and predicted proportions for plotting
    true_flat = []
    pred_flat = []
    
    for sample in pseudobulk.obs_names:
        for cell_type in true_props.columns:
            true_flat.append(true_props.at[sample, cell_type])
            pred_flat.append(deconv_results.at[sample, cell_type])
    
    plt.scatter(true_flat, pred_flat, alpha=0.6)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlabel('True Proportion')
    plt.ylabel('Predicted Proportion')
    plt.title(f'Deconvolution Validation\nPearson r = {mean_corr:.3f}, RMSE = {mean_rmse:.3f}')
    plt.tight_layout()
    plt.savefig(os.path.join(sc.settings.figdir, 'deconvolution_validation.pdf'))
    plt.close()
    
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


def main():
    """Main function to run the deconvolution pipeline."""
    # Parse arguments
    args = parse_arguments()
    
    # Set up environment
    plots_dir, custom_palette = setup_environment(args)
    
    # Log start of processing
    logger.info("Starting cell type deconvolution pipeline")
    logger.info(f"Bulk data: {args.bulk_path}")
    logger.info(f"Reference data: {args.ref_path}")
    logger.info(f"Output directory: {args.output_dir}")
    
    try:
        # Load and validate data
        bulk_adata, ref_adata = load_and_validate_data(args.bulk_path, args.ref_path)
        
        # Find shared genes
        shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names).tolist()
        logger.info(f"Using {len(shared_genes)} shared genes")
        
        # Identify cell type markers
        markers_dict = identify_cell_type_markers(
            ref_adata, 
            args.annotation_key, 
            n_markers=args.n_markers
        )
        
        # Create signature matrix
        signature_matrix = create_signature_matrix(
            ref_adata, 
            markers_dict, 
            args.annotation_key, 
            shared_genes
        )
        
        # Save signature matrix
        signature_matrix.to_csv(os.path.join(args.output_dir, 'signature_matrix.csv'))
        
        # Validate deconvolution approach
        validation_score = validate_deconvolution(
            ref_adata,
            args.annotation_key,
            signature_matrix,
            shared_genes,
            n_samples=20
        )
        
        # Run deconvolution with confidence intervals
        deconv_results, lower_ci, upper_ci = bootstrap_confidence_intervals(
            bulk_adata,
            signature_matrix,
            shared_genes,
            n_bootstrap=args.n_bootstrap
        )
        
        # Create visualizations
        plot_deconvolution_results(deconv_results, lower_ci, upper_ci, custom_palette)
        
        # Save results
        deconv_results.to_csv(os.path.join(args.output_dir, 'deconvolution_results.csv'))
        lower_ci.to_csv(os.path.join(args.output_dir, 'deconvolution_lower_ci.csv'))
        upper_ci.to_csv(os.path.join(args.output_dir, 'deconvolution_upper_ci.csv'))
        
        # Create summary report
        with open(os.path.join(args.output_dir, 'deconvolution_summary.txt'), 'w') as f:
            f.write("Cell Type Deconvolution Summary\n")
            f.write("===============================\n\n")
            f.write(f"Bulk dataset: {args.bulk_path}\n")
            f.write(f"Reference dataset: {args.ref_path}\n")
            f.write(f"Number of bulk samples: {bulk_adata.shape[0]}\n")
            f.write(f"Number of reference cells: {ref_adata.shape[0]}\n")
            f.write(f"Number of shared genes: {len(shared_genes)}\n")
            f.write(f"Number of cell types: {len(signature_matrix.columns)}\n\n")
            f.write(f"Validation correlation score: {validation_score:.3f}\n\n")
            
            f.write("Top cell types by average proportion:\n")
            for cell_type, prop in deconv_results.mean().sort_values(ascending=False).items():
                f.write(f"  {cell_type}: {prop:.3f}\n")
        
        logger.info("Deconvolution pipeline completed successfully")
        logger.info(f"Results saved to {args.output_dir}")
        
    except Exception as e:
        logger.error(f"Deconvolution pipeline failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1
    
    return 0


if __name__ == "__main__":
    main()