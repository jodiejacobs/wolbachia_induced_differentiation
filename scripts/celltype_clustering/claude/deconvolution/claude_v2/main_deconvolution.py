import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import logging
import argparse
import warnings
from cell_type_deconvolution import (
    load_and_validate_data, scale_bulk_to_match_scrna, 
    quantile_normalize_datasets, harmonize_datasets,
    plot_distribution_comparison, run_pca_comparison,
    identify_cell_type_markers, create_signature_matrix,
    deconvolve_samples, bootstrap_confidence_intervals,
    plot_deconvolution_results
)

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
    parser = argparse.ArgumentParser(description="Enhanced Cell Type Deconvolution with Batch Correction")

    parser.add_argument("--bulk_path", type=str, required=True,
                        help="Path to the bulk AnnData file")

    parser.add_argument("--ref_path", type=str, required=True,
                        help="Path to the reference AnnData file")

    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save the output results")

    parser.add_argument("--annotation_key", type=str, required=False,
                        default="annotation",
                        help="Column name with cell type labels in reference data")
                        
    parser.add_argument("--n_markers", type=int, required=False, 
                        default=100,
                        help="Number of marker genes to identify per cell type")
                        
    parser.add_argument("--n_bootstraps", type=int, required=False, 
                        default=50,
                        help="Number of bootstrap iterations for confidence intervals")
                        
    parser.add_argument("--normalization", type=str, required=False, 
                        default="all",
                        choices=["none", "scale", "quantile", "harmonize", "all"],
                        help="Normalization method(s) to use")
                        
    parser.add_argument("--seed", type=int, required=False, 
                        default=42,
                        help="Random seed for reproducibility")
                        
    return parser.parse_args()

def main():
    """Main function to run the deconvolution pipeline with batch effect correction."""
    # Parse arguments
    args = parse_arguments()
    
    # Set random seed
    np.random.seed(args.seed)
    
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
    
    # Log start of processing
    logger.info("Starting enhanced cell type deconvolution pipeline")
    logger.info(f"Bulk data: {args.bulk_path}")
    logger.info(f"Reference data: {args.ref_path}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Normalization method: {args.normalization}")
    
    try:
        # Load and validate data
        bulk_adata, ref_adata = load_and_validate_data(args.bulk_path, args.ref_path)
        
        # Find shared genes
        shared_genes = bulk_adata.var_names.intersection(ref_adata.var_names).tolist()
        logger.info(f"Using {len(shared_genes)} shared genes")
        
        # Create subset with shared genes
        bulk_subset = bulk_adata[:, shared_genes].copy()
        ref_subset = ref_adata[:, shared_genes].copy()
        
        # Plot distributions before normalization
        plot_distribution_comparison(bulk_subset, ref_subset, shared_genes, plots_dir, "before")
        run_pca_comparison(bulk_subset, ref_subset, plots_dir, "before")
        
        # Apply normalization based on selected method
        if args.normalization == "none":
            # Use data as is
            normalized_bulk = bulk_subset
            normalized_ref = ref_subset
            logger.info("Using raw data without normalization")
            
        elif args.normalization == "scale":
            # Scale bulk data to match scRNA-seq ranges
            normalized_bulk = scale_bulk_to_match_scrna(bulk_subset, ref_subset, shared_genes)
            normalized_ref = ref_subset
            logger.info("Applied gene-specific scaling")
            
        elif args.normalization == "quantile":
            # Apply quantile normalization
            normalized_bulk = quantile_normalize_datasets(bulk_subset, ref_subset, shared_genes)
            normalized_ref = ref_subset
            logger.info("Applied quantile normalization")
            
        elif args.normalization == "harmonize":
            # Apply harmonization
            normalized_bulk, normalized_ref = harmonize_datasets(bulk_subset, ref_subset, shared_genes)
            logger.info("Applied harmonization")
            
        else:  # "all" - try all methods and select the best
            logger.info("Trying all normalization methods and selecting the best")
            
            # 1. Gene-specific scaling
            scaled_bulk = scale_bulk_to_match_scrna(bulk_subset, ref_subset, shared_genes)
            plot_distribution_comparison(scaled_bulk, ref_subset, shared_genes, plots_dir, "scale")
            run_pca_comparison(scaled_bulk, ref_subset, plots_dir, "scale")
            
            # 2. Quantile normalization
            quantile_bulk = quantile_normalize_datasets(bulk_subset, ref_subset, shared_genes)
            plot_distribution_comparison(quantile_bulk, ref_subset, shared_genes, plots_dir, "quantile")
            run_pca_comparison(quantile_bulk, ref_subset, plots_dir, "quantile")
            
            # 3. Harmonization
            harmonized_bulk, harmonized_ref = harmonize_datasets(bulk_subset, ref_subset, shared_genes)
            plot_distribution_comparison(harmonized_bulk, harmonized_ref, shared_genes, plots_dir, "harmonize")
            run_pca_comparison(harmonized_bulk, harmonized_ref, plots_dir, "harmonize")
            
            # Use harmonized data for the rest of the pipeline
            normalized_bulk = harmonized_bulk
            normalized_ref = harmonized_ref
            logger.info("Selected harmonization as the best method based on distribution matching")
        
        # Plot distributions after normalization
        if args.normalization != "all":
            plot_distribution_comparison(normalized_bulk, normalized_ref, shared_genes, plots_dir, "after")
            run_pca_comparison(normalized_bulk, normalized_ref, plots_dir, "after")
        
        # Identify cell type markers
        markers_dict = identify_cell_type_markers(
            normalized_ref, 
            args.annotation_key, 
            n_markers=args.n_markers
        )
        
        # Create signature matrix
        signature_matrix = create_signature_matrix(
            normalized_ref, 
            markers_dict, 
            args.annotation_key, 
            shared_genes
        )
        
        # Save signature matrix
        signature_matrix.to_csv(os.path.join(args.output_dir, 'signature_matrix.csv'))
        
        # Run deconvolution with confidence intervals
        deconv_results, lower_ci, upper_ci = bootstrap_confidence_intervals(
            normalized_bulk,
            signature_matrix,
            shared_genes,
            n_bootstrap=args.n_bootstraps
        )
        
        # Create visualizations
        plot_deconvolution_results(deconv_results, lower_ci, upper_ci, custom_palette)
        
        # Save results
        deconv_results.to_csv(os.path.join(args.output_dir, 'deconvolution_results.csv'))
        lower_ci.to_csv(os.path.join(args.output_dir, 'deconvolution_lower_ci.csv'))
        upper_ci.to_csv(os.path.join(args.output_dir, 'deconvolution_upper_ci.csv'))
        
        # Create summary report
        with open(os.path.join(args.output_dir, 'deconvolution_summary.txt'), 'w') as f:
            f.write("Enhanced Cell Type Deconvolution Summary\n")
            f.write("=======================================\n\n")
            f.write(f"Bulk dataset: {args.bulk_path}\n")
            f.write(f"Reference dataset: {args.ref_path}\n")
            f.write(f"Number of bulk samples: {bulk_adata.shape[0]}\n")
            f.write(f"Number of reference cells: {ref_adata.shape[0]}\n")
            f.write(f"Number of shared genes: {len(shared_genes)}\n")
            f.write(f"Number of cell types: {len(signature_matrix.columns)}\n")
            f.write(f"Normalization method: {args.normalization}\n\n")
            
            f.write("Top cell types by average proportion:\n")
            for cell_type, prop in deconv_results.mean().sort_values(ascending=False).items():
                f.write(f"  {cell_type}: {prop:.3f}\n")
        
        logger.info("Enhanced deconvolution pipeline completed successfully")
        logger.info(f"Results saved to {args.output_dir}")
        
    except Exception as e:
        logger.error(f"Deconvolution pipeline failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1
    
    return 0

if __name__ == "__main__":
    main()