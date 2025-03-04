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

# Argument parsing
parser = argparse.ArgumentParser(description="Scanpy Analysis - Single Bootstrap Iteration")

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

parser.add_argument("--atlas", "-b", type=str, required=False, 
                    default="fly_cell_atlas", 
                    help="Atlas type")

args = parser.parse_args()

ref_adata_path = args.ref_adata_path
output_folder = args.output_folder
annotation = args.annotation  # Column with cell type labels
mem = args.mem
atlas = args.atlas 

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

## Color map 
color_dict={
    'JW18DOX':'#C7E576',
    'JW18wMel':'#241E4E',
    'S2DOX':'#F4E04D',
    'S2wMel':'#587792'

}

def preprocess(adata):
    """Preprocess the AnnData object."""
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=6)

def subsample_celltypes(adata, annotation):
    """Subsample cell types to ensure even representation."""
    
    # Reduce memory footprint <- these can be dropped
    adata.uns = {}
    adata.obsm = {}
    adata.varm = {}
    adata.obsp = {}

    celltype_counts = adata.obs[annotation].value_counts()
    min_count = celltype_counts.min()
    subsampled_adatas = []

    for celltype in celltype_counts.index:
        celltype_indices = adata.obs[annotation] == celltype
        subsampled_adata = adata[celltype_indices]  # No copy yet

        # Fix: Use obs_names for efficient subsetting
        subsampled_indices = np.random.choice(subsampled_adata.obs_names, min_count, replace=False)
        subsampled_adata = subsampled_adata[subsampled_indices, :].copy()  # Shallow copy
        
        subsampled_adatas.append(subsampled_adata)

    # Concatenate efficiently
    subsampled_adata = anndata.concat(
        subsampled_adatas, label="subsampled", keys=[f"sample_{i}" for i in range(len(subsampled_adatas))], index_unique="-"
    )

    return subsampled_adata

# Load and process
ref_adata = sc.read_h5ad(ref_adata_path)

# Remove raw data if not needed
ref_adata.raw = None

# Only keep highly variable genes if applicable
sc.pp.highly_variable_genes(ref_adata, n_top_genes=4000)
ref_adata = ref_adata[:, ref_adata.var.highly_variable]

preprocess(ref_adata)
ref_adata = subsample_celltypes(ref_adata, annotation)

# Fix reserved column name issue
ref_adata.var.reset_index(drop=True, inplace=True)

# Save output
ref_adata.write_h5ad(f"{output_folder}/{atlas}_subsampled.h5ad")