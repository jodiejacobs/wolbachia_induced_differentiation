#!/bin/bash
#SBATCH --job-name=scanpy_analysis
#SBATCH --partition=medium
#SBATCH --mem=1024G
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=64
#SBATCH --output=scanpy_analysis_%j.out
#SBATCH --error=scanpy_analysis_%j.err

# Load your modules or activate your conda environment
source ~/.bashrc
conda activate scanpy

# Set environment variables for memory management
export PYTHONUNBUFFERED=1  # Ensure Python output is not buffered
export NUMBA_CACHE_DIR=/tmp/numba_cache_$SLURM_JOB_ID  # Numba cache in temp directory
export MPLCONFIGDIR=/tmp/matplotlib_$SLURM_JOB_ID  # Matplotlib config in temp
mkdir -p $NUMBA_CACHE_DIR $MPLCONFIGDIR

# Print job info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Available memory: $(free -h)"

# Run the Python script
python /private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/20250303_scanpy_cKDtree_integration_3.py

# Clean up temp directories
rm -rf $NUMBA_CACHE_DIR $MPLCONFIGDIR

echo "Job finished at: $(date)"