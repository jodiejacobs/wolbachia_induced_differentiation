#Hard coded files:
fasta_file_path = '/scratch1/jodie/wolbachia/Micro-C/09Nov2023_Micro_C/map/final_enriched_genes/enriched_genes.fasta'
fasta_out_path = '/scratch1/jodie/wolbachia/Wolbachia_induced_differentiation/enriched_genes_flybase_id.fasta'
key_file_path = '/scratch1/jodie/wolbachia/Micro-C/reference_genomes/dmel_gene_id_key.uniq.tsv'


import pandas as pd

key = pd.read_csv(key_file_path, sep='\s+', header=0)
print(key)
# # Set a column as index
key = key.set_index('gene_symbol')

# Make the index unique (if necessary)
# This step will keep the first occurrence and drop the rest in case of duplicates
key = key[~key.index.duplicated(keep='first')]

# # Now you can search using the index
# # For example, to get the row for index 'CG17291'
# gene_ex_id = "CG17291"
# gene_ex_trx = key.loc[gene_ex_id]['transcript_id']
# gene_ex_sym = key.loc[gene_ex_id]['gene_id']

# print(f'Gene_ID: {gene_ex_id}, Transcript: {gene_ex_trx}, Gene_Symbol: {gene_ex_sym}')

# # Reading the FASTA file

def rename(gene_symbol, key):
    try:
        flybase_id = key.loc[gene_symbol]['transcript_id']
        return flybase_id
    except KeyError:
        return gene_symbol
        # print(gene_symbol)
        
fasta_renamed = []
with open(fasta_file_path, 'r') as file:
    for line in file:
        # filtered_df = df[(df['column_name'] > 2) & (df['column_name'] < 5)]
        gene = line.strip('\n')
        if gene[0][0] == ">":
            fasta_renamed.append(gene)
        else:
            try:
                gene_tx = rename(gene, key)
                fasta_renamed.append(gene_tx)
                # print(f'Gene_ID: {gene_ex_id}, Transcript: {gene_ex_trx}, Gene_Symbol: {gene_ex_sym}')
            except KeyError:
                # fasta_headers.append(gene)
                continue

print(fasta_renamed)
# Open the file in write mode
with open(fasta_out_path, 'w') as file:
    # Write each item on a new line
    for item in fasta_renamed:
        file.write(f"{item}\n")


#Convert output of scanpy clustering ranking to Flybase ID's
scanpy_genes = "/scratch1/jodie/wolbachia/celltype_clustering/12_20_23_gene_rankings.csv"
scanpy_genes_out = "/scratch1/jodie/wolbachia/Wolbachia_induced_differentiation/12_20_23_scanpy_gene_rankings.csv"
scanpy_genes_df = pd.read_csv(scanpy_genes, sep='\t', index_col=0)

scanpy_genes_df = scanpy_genes_df.applymap(lambda x: rename(x, key))

# Convert DataFrame to CSV
scanpy_genes_df.to_csv(scanpy_genes_out, index=True)


#Convert output of DESeq2 to Flybase IDs:

# Change the index of the key to gene_id
key = key.set_index('gene_id')


deseq2_out= '/scratch1/jodie/wolbachia/Wolbachia_induced_differentiation/Kallisto-DESeq2_WaldTest_infection_wMel_vs_DOX_flybase_id.csv'
       
deseq2_df=pd.read_csv(deseq2_results,sep=',')
deseq2_df['Unnamed: 0'] = deseq2_df['Unnamed: 0'].str.replace('Dmel_', '', regex=False)
deseq2_df['Unnamed: 0'] = deseq2_df['Unnamed: 0'].apply(lambda x: rename(x, key))


print(deseq2_df)
# print(deseq2_df)
deseq2_df.to_csv(deseq2_out, index=False)

