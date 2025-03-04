import argparse
import scanpy as sc
import pandas as pd
import numpy as np

# Set up argument parser
parser = argparse.ArgumentParser(description="Build a Scanpy AnnData object from CSV files and save as .h5ad")
parser.add_argument("-f", "--file_name", required=True, help="Base name of the CSV files (e.g., 'your_file_name')")
args = parser.parse_args()

# Retrieve the base file name from arguments
file_name = args.file_name

# Define file paths
count_matrix_file = f"{file_name}_count_matrix.csv"
cell_metadata_file = f"{file_name}_cell_metadata.csv"
gene_metadata_file = f"{file_name}_gene_metadata.csv"
subtypes_file = f"{file_name}_subtypes.csv"

# Load count matrix
print("Loading count matrix...")
count_matrix = pd.read_csv(count_matrix_file, sep=',', index_col=0)
count_matrix = count_matrix.T  # Transpose to have cells as rows and genes as columns
count_matrix.rename(index = {'Gene':'gene_symbol'})

# Extract Column and Row names 
Gene = count_matrix.columns.astype(str).to_list()
CellID = count_matrix.index

# Export Values to numpy dataframe as float
df = count_matrix.to_numpy().astype(int)

# Create a new AnnData object
adata = sc.AnnData(X=df, dtype=df.dtype) 

# Set the cell IDs as observation names
adata.obs_names = CellID
adata.var_names = Gene 

# Load cell metadata
print("Loading cell metadata...")
cell_metadata = pd.read_csv(cell_metadata_file, index_col=0)
adata.obs = cell_metadata.astype(str)

# Load gene metadata 
print("Loading gene metadata...")
gene_metadata = pd.read_csv(gene_metadata_file, index_col=0)
adata.var["gene_short_name"] = gene_metadata.astype(str)

# Load subtypes (cell types) #This works 
print("Loading subtypes...")
subtypes = pd.read_csv(subtypes_file, index_col=0)
subtypes.rename(columns={"colData(cds)$subtypes": "subtypes"}, inplace=True)
adata.obs["subtypes"] = subtypes["subtypes"].values.astype(str)

# Save the AnnData object as an .h5ad file
output_file = f"{file_name}.h5ad"
print(f"Saving AnnData object to {output_file}...")
adata.write(output_file)

print("Done!")
