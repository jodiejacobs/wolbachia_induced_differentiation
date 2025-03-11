#!/bin/bash
#SBATCH --job-name=bulk_sc_integration
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/output/logs/bulk_sc_integration_%J.out
#SBATCH --error=/path/to/output/logs/bulk_sc_integration_%J.err
#SBATCH --mail-type=END,FAIL

# Create logs directory if it doesn't exist
mkdir -p /private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/MNN_kNN_tree/logs

# Define paths - EDIT THESE
SCRIPT_PATH="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/deconvolution/cell_type_deconvolution.py"
BULK_ADATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/bulk_adata.h5ad"
REF_ADATA="/private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/scanpy_objects/blood_adata_all.h5ad"
OUTPUT_DIR="/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/celltype_clustering/claude/MNN_kNN_tree/myeloid_atlas"
ANNOTATION_KEY="subclustering"  # Column in reference data with cell type labels

# Activate conda/mamba environment
eval "$(conda shell.bash hook)"
conda activate scanpy  # Replace with your environment name

# Copy the script to the output folder for reproducibility
mkdir -p "$OUTPUT_DIR"
cp "$SCRIPT_PATH" "$OUTPUT_DIR/$(basename $SCRIPT_PATH)"

echo "Starting bulk to single-cell integration at $(date)"
echo "===================================================="
echo "Bulk dataset: $BULK_ADATA"
echo "Reference dataset: $REF_ADATA"
echo "Output directory: $OUTPUT_DIR"
echo "Annotation key: $ANNOTATION_KEY"
echo "===================================================="

# Run the integration script
python "$SCRIPT_PATH" \
  --bulk_path "$BULK_ADATA" \
  --ref_path "$REF_ADATA" \
  --output_dir "$OUTPUT_DIR" \
  --annotation_key "$ANNOTATION_KEY" \
  --num_permutations 1000 \
  --seed 42 > "$OUTPUT_DIR/integration.out" 2> "$OUTPUT_DIR/integration.err"

# Check if the script executed successfully
if [ $? -eq 0 ]; then
  echo "Bulk to single-cell integration completed successfully at $(date)"
else
  echo "Integration failed at $(date)"
  echo "Attempting fallback with fewer neighbors and permutations..."
  
  # Try with explicit k value as fallback
  python "$SCRIPT_PATH" \
    --bulk_path "$BULK_ADATA" \
    --ref_path "$REF_ADATA" \
    --output_dir "${OUTPUT_DIR}_fallback" \
    --annotation_key "$ANNOTATION_KEY" \
    --k_neighbors 20 \
    --num_permutations 500
    
  if [ $? -eq 0 ]; then
    echo "Fallback integration completed successfully at $(date)"
  else
    echo "Fallback integration also failed at $(date)"
  fi
fi

# Create a simple report with the results
echo "======= Bulk to Single-cell Integration Report =======" > "$OUTPUT_DIR/integration_report.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/integration_report.txt"
echo "Bulk dataset: $BULK_ADATA" >> "$OUTPUT_DIR/integration_report.txt"
echo "Reference dataset: $REF_ADATA" >> "$OUTPUT_DIR/integration_report.txt"
echo "Output directory: $OUTPUT_DIR" >> "$OUTPUT_DIR/integration_report.txt"
echo "Annotation key: $ANNOTATION_KEY" >> "$OUTPUT_DIR/integration_report.txt"
echo "=======================================" >> "$OUTPUT_DIR/integration_report.txt"

# List created files
echo "" >> "$OUTPUT_DIR/integration_report.txt"
echo "Generated files:" >> "$OUTPUT_DIR/integration_report.txt"
ls -l "$OUTPUT_DIR/" | grep -v "logs" >> "$OUTPUT_DIR/integration_report.txt"

# Check for results file and count classifications
if [ -f "$OUTPUT_DIR/bulk_classification_results.csv" ]; then
  echo "" >> "$OUTPUT_DIR/integration_report.txt"
  echo "Classification summary:" >> "$OUTPUT_DIR/integration_report.txt"
  TOTAL_SAMPLES=$(wc -l < "$OUTPUT_DIR/bulk_classification_results.csv")
  echo "Total samples classified: $((TOTAL_SAMPLES - 1)) (excluding header)" >> "$OUTPUT_DIR/integration_report.txt"
  
  # Count samples by cell type
  echo "Cell type distribution:" >> "$OUTPUT_DIR/integration_report.txt"
  tail -n +2 "$OUTPUT_DIR/bulk_annotations.csv" | cut -d',' -f2 | sort | uniq -c | sort -nr >> "$OUTPUT_DIR/integration_report.txt"
fi

echo "Job completed. See logs and report in $OUTPUT_DIR"