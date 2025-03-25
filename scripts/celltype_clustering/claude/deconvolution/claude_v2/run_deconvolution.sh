#!/bin/bash
#SBATCH --job-name=deconv_enhanced
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=1024gb
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/logs/deconv_enhanced_%J.out
#SBATCH --error=/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/logs/deconv_enhanced_%J.err
#SBATCH --mail-user=jomojaco@ucsc.edu
#SBATCH --mail-type=END,FAIL


# Create logs directory if it doesn't exist
mkdir -p /private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/logs

# Define paths
SCRIPT_DIR="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/claude_v2"
BULK_DATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad"
REF_DATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/fca_subset.h5ad"
OUTPUT_DIR="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/claude_v2/fca_subset"
ANNOTATION_KEY="annotation"
N_MARKERS=100
N_BOOTSTRAPS=50
NORMALIZATION="all"

# # # Activate conda environment
# eval "$(conda shell.bash hook)"
# conda activate scanpy

echo "Starting enhanced cell type deconvolution at $(date)"
echo "===================================================="
echo "Bulk dataset: $BULK_DATA"
echo "Reference dataset: $REF_DATA"
echo "Output directory: $OUTPUT_DIR"
echo "Annotation column: $ANNOTATION_KEY"
echo "Normalization method: $NORMALIZATION"
echo "===================================================="

# Run the deconvolution script
python $SCRIPT_DIR/main_deconvolution.py \
  --bulk_path "$BULK_DATA" \
  --ref_path "$REF_DATA" \
  --output_dir "$OUTPUT_DIR" \
  --annotation_key "$ANNOTATION_KEY" \
  --n_markers "$N_MARKERS" \
  --n_bootstraps "$N_BOOTSTRAPS" \
  --normalization "$NORMALIZATION" \
  --seed 42

# Check if the script executed successfully
if [ $? -eq 0 ]; then
  echo "Enhanced cell type deconvolution completed successfully at $(date)"
else
  echo "Enhanced cell type deconvolution failed at $(date)"
  echo "Attempting fallback with simpler normalization method..."
  
  # Try with a simpler normalization method as fallback
  python $SCRIPT_DIR/main_deconvolution.py \
    --bulk_path "$BULK_DATA" \
    --ref_path "$REF_DATA" \
    --output_dir "${OUTPUT_DIR}_fallback" \
    --annotation_key "$ANNOTATION_KEY" \
    --n_markers "$N_MARKERS" \
    --n_bootstraps 20 \
    --normalization "quantile" \
    --seed 42
  
  if [ $? -eq 0 ]; then
    echo "Fallback deconvolution completed successfully at $(date)"
  else
    echo "Fallback deconvolution also failed at $(date)"
  fi
fi

# Create summary of results
echo "======= Cell Type Deconvolution Summary =======" > "$OUTPUT_DIR/deconvolution_report.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Bulk dataset: $BULK_DATA" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Reference dataset: $REF_DATA" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Output directory: $OUTPUT_DIR" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Annotation column: $ANNOTATION_KEY" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Normalization method: $NORMALIZATION" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "=======================================" >> "$OUTPUT_DIR/deconvolution_report.txt"

# List created files
echo "" >> "$OUTPUT_DIR/deconvolution_report.txt"
echo "Generated files:" >> "$OUTPUT_DIR/deconvolution_report.txt"
ls -l "$OUTPUT_DIR/" >> "$OUTPUT_DIR/deconvolution_report.txt"

echo "Job completed. See logs and report in $OUTPUT_DIR"