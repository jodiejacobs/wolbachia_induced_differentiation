#!/bin/bash
#SBATCH --job-name=cell_deconv
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/output/logs/cell_deconv_%J.out
#SBATCH --error=/path/to/output/logs/cell_deconv_%J.err
#SBATCH --mail-type=END,FAIL

# Create logs directory if it doesn't exist
mkdir -p /private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/logs

# Define paths - EDIT THESE
SCRIPT_PATH="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/cell_type_deconvolution.py"
BULK_ADATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad"
REF_ADATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/blood_adata_all.h5ad"
OUTPUT_DIR="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/myeloid_atlas"
ANNOTATION_KEY="subclustering"  # Column in reference data with cell type labels

# # Activate conda/mamba environment
# eval "$(conda shell.bash hook)"
# conda activate scanpy  # Replace with your environment name

# # Copy the script to the output folder for reproducibility
# mkdir -p "$OUTPUT_DIR"
# cp "$SCRIPT_PATH" "$OUTPUT_DIR/$(basename $SCRIPT_PATH)"

# echo "Starting cell type deconvolution at $(date)"
# echo "===================================================="
# echo "Bulk dataset: $BULK_ADATA"
# echo "Reference dataset: $REF_ADATA"
# echo "Output directory: $OUTPUT_DIR"
# echo "Annotation key: $ANNOTATION_KEY"
# echo "===================================================="

# Run the deconvolution script with standard parameters
python "$SCRIPT_PATH" \
  --bulk_path "$BULK_ADATA" \
  --ref_path "$REF_ADATA" \
  --output_dir "$OUTPUT_DIR" \
  --annotation_key "$ANNOTATION_KEY" \
  --n_markers 100 \
  --n_bootstrap 500 \
  --seed 42 > "$OUTPUT_DIR/deconvolution.out" 2> "$OUTPUT_DIR/deconvolution.err"    

# Check if the script executed successfully
if [ $? -eq 0 ]; then
  echo "Cell type deconvolution completed successfully at $(date)"
else
  echo "Deconvolution failed at $(date)"
  echo "Attempting fallback with reduced parameters..."
  
  # Try with fewer markers and bootstrap iterations as fallback
  python "$SCRIPT_PATH" \
    --bulk_path "$BULK_ADATA" \
    --ref_path "$REF_ADATA" \
    --output_dir "${OUTPUT_DIR}_fallback" \
    --annotation_key "$ANNOTATION_KEY" \
    --n_markers 50 \
    --n_bootstrap 100
    
  if [ $? -eq 0 ]; then
    echo "Fallback deconvolution completed successfully at $(date)"
  else
    echo "Fallback deconvolution also failed at $(date)"
  fi
fi

# Create a simple report with the results summary
if [ -f "$OUTPUT_DIR/deconvolution_summary.txt" ]; then
  echo "=== DECONVOLUTION RESULTS SUMMARY ===" > "$OUTPUT_DIR/slurm_report.txt"
  echo "Date: $(date)" >> "$OUTPUT_DIR/slurm_report.txt"
  echo "" >> "$OUTPUT_DIR/slurm_report.txt"
  cat "$OUTPUT_DIR/deconvolution_summary.txt" >> "$OUTPUT_DIR/slurm_report.txt"
  
  echo "" >> "$OUTPUT_DIR/slurm_report.txt"
  echo "Plot files generated:" >> "$OUTPUT_DIR/slurm_report.txt"
  ls -l "$OUTPUT_DIR/plots/" >> "$OUTPUT_DIR/slurm_report.txt"
fi

# Check for memory usage summary
MAX_MEM_USED=$(grep "MaxRSS" /proc/self/status | awk '{print $2}')
if [ ! -z "$MAX_MEM_USED" ]; then
  MAX_MEM_GB=$(echo "scale=2; $MAX_MEM_USED / 1024 / 1024" | bc)
  echo "" >> "$OUTPUT_DIR/slurm_report.txt"
  echo "Maximum memory used: ${MAX_MEM_GB}GB" >> "$OUTPUT_DIR/slurm_report.txt"
fi

echo "Job completed. See logs and report in $OUTPUT_DIR"