# Run this script to convert a Seurat object saved as an .rds file to an .h5ad file using SeuratDisk
# Requires: mamba activate seurat 
# Example: Rscript seurat_to_h5ad.R -f seurat_obj.rds
# Example: find /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/data/atlas/peng_2024_fruitfly_organogenesis/ -type f -name '*.rds' | parallel Rscript /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/seurat_to_h5ad.R -f {}


library(Seurat)
library(SeuratDisk)
library(monocle3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] != "-f" || length(args) < 2) {
    stop("Please provide the file name as a command line argument using -f.")
}

input_file <- args[2]

file_name <- sub("\\.rds$", "", input_file)

#Load the data as a monocle3 object
cds <- readRDS(input_file)

# Extract count matrix
count_data <- counts(cds)
write.csv(as.matrix(count_data), paste0(file_name, "_count_matrix.csv"))

# Extract cell metadata
cell_metadata <- as.data.frame(colData(cds))
write.csv(cell_metadata, paste0(file_name, "_cell_metadata.csv"))

# Extract gene metadata
gene_metadata <- as.data.frame(rowData(cds))
write.csv(gene_metadata, paste0(file_name, "_gene_metadata.csv"))

# Extract cell type annotations 
subtypes <- as.data.frame(colData(cds)$subtypes)
write.csv(subtypes, paste0(file_name, "_subtypes.csv"))

# run python script to convert csv to h5ad
system(paste("python /private/groups/russelllab/jodie/wolbachia_induced_DE/scanpy_clustering/pub_scipts/csv_to_h5ad.py -f", file_name))