import os
import pandas as pd
import argparse 

parser = argparse.ArgumentParser(description='Filter CSV files in a directory.')
parser.add_argument('--directory', type=str, default="/private/groups/russelllab/jodie/wolbachia_induced_DE/v04_deseq2/plots",help='Directory containing CSV files')

# Parse the arguments
args = parser.parse_args()


def filter_csv_files(directory):
    # Iterate over all files in the given directory
    for filename in os.listdir(directory):
        # Check if the file is a CSV
        if filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            # Read the CSV file
            df = pd.read_csv(file_path)

            # Apply the filtering
            pos_filtered_df = df[(df['padj'] < 0.001) & (df['log2FoldChange'] > 2)]
            neg_filtered_df = df[(df['padj'] < 0.001) & (df['log2FoldChange'] < -2)]

            pos_filtered = ['log2FC>2, ' + s for s in pos_filtered_df['Unnamed: 0'].tolist()]
            neg_filtered = ['log2FC<-2, '+ s for s in neg_filtered_df['Unnamed: 0'].tolist()]

            # Extract gene names
            filtered_genes = pos_filtered + neg_filtered
            filtered_genes = [gene.replace("Dmel_C", "C") for gene in filtered_genes]

            # Output the result
            output_file = os.path.join(directory, 'GO',f'GO_Search_{filename}')
            with open(output_file, 'w') as f:
                for gene in filtered_genes:
                    f.write(f"{gene}\n")
            
            print(f"Processed {filename} - {len(filtered_genes)} genes found")


# Run the filtering function
filter_csv_files(args.directory)
