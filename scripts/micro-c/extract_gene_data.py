import pandas as pd
import pyBigWig
import pysam
import re
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Extract gene data with WGCNA modules')
parser.add_argument('--chrom', type=str, default='2L', help='Chromosome')
parser.add_argument('--start', type=int, default=6000000, help='Start position')
parser.add_argument('--end', type=int, default=9000000, help='End position')
parser.add_argument('--output', type=str, default='gene_data.tsv', help='Output file')
args = parser.parse_args()

# Define region of interest from args
chrom = args.chrom
start_bp = args.start
end_bp = args.end
output_file = args.output

# Paths to files
gtf_file = '/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/micro-c/WGCNA_bigwig/Dmel6wMel2978.6.clean.gtf'
wgcna_modules_file = '/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/micro-c/WGCNA_bigwig/significant_DE_WGCNA_modules.Dmel.tsv'

# BigWig files to extract data from
bw_wmel = "/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/mapped/tracks/JW18-wMel_1000.bw"
bw_dox = "/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/mapped/tracks/JW18-DOX_1000.bw"
rnaseq_cov = '/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/mapped/tracks/JW18wMel221117_vs_DOX.diff.coverage.bw'
microc_cov = '/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/mapped/tracks/JW18-wMel_vs_DOX.diff.coverage.bw'
mappability = '/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/mapped/tracks/dmel_6_genmap_e0_k150_final.bw'

# Map D. melanogaster chromosome names to RefSeq accessions
chrom_to_refseq = {
    '2L': 'NC_004354.4',  # Chromosome 2L
    '2R': 'NC_004355.4',  # Chromosome 2R
    '3L': 'NC_004356.4',  # Chromosome 3L
    '3R': 'NC_004357.4',  # Chromosome 3R
    '4': 'NC_004353.4',   # Chromosome 4
    'X': 'NC_004354.4',   # Chromosome X
    'Y': 'NC_024512.1'    # Chromosome Y
}

# Extract CG number from gene name or ID
def extract_cg_number(gene_string):
    match = re.search(r'(CG\d+)', str(gene_string))
    if match:
        return match.group(1)
    return None

# Load WGCNA module assignments
def load_wgcna_data(file_path):
    try:
        # Load the complete WGCNA data file
        wgcna_df = pd.read_csv(file_path, sep='\t')
        
        print("WGCNA file preview:")
        print(wgcna_df.head())
        
        # Identify gene ID column (usually first column)
        gene_id_col = wgcna_df.columns[0]
        
        # Create dictionaries mapping gene IDs to WGCNA data
        gene_to_wgcna_data = {}
        cg_to_wgcna_data = {}
        
        for _, row in wgcna_df.iterrows():
            gene_id = str(row[gene_id_col])
            row_dict = row.to_dict()
            
            # Store by original ID
            gene_to_wgcna_data[gene_id] = row_dict
            
            # Store by CG number directly - this is already a CG number in your case
            cg_to_wgcna_data[gene_id] = row_dict
        
        print(f"Loaded {len(gene_to_wgcna_data)} original gene IDs from WGCNA data")
        
        return gene_to_wgcna_data, cg_to_wgcna_data, list(wgcna_df.columns)
    except Exception as e:
        print(f"Error loading WGCNA data: {e}")
        return {}, {}, []

# Function to extract gene information from GTF
def extract_genes_from_gtf(gtf_file, chromosome, start, end):
    genes = []
    gene_count = 0
    
    # Convert from generic chromosome name to RefSeq accession
    refseq_chrom = chrom_to_refseq.get(chromosome)
    if not refseq_chrom:
        print(f"Warning: No RefSeq mapping for chromosome {chromosome}")
        refseq_chrom = chromosome  # Fall back to original name
    
    print(f"Looking for {refseq_chrom} in GTF file (mapped from {chromosome})")
    
    # Keep track of genes that have been processed
    processed_genes = set()
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                # Check if on the right chromosome
                if fields[0] != refseq_chrom:
                    continue
                
                # Parse gene_id
                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                if not gene_id_match:
                    continue
                
                gene_id = gene_id_match.group(1)
                
                # Only process gene features or collect unique genes
                if fields[2] == 'gene' or gene_id not in processed_genes:
                    feature_start = int(fields[3])
                    feature_end = int(fields[4])
                    
                    # Skip if the feature is outside our region of interest
                    if feature_end < start or feature_start > end:
                        continue
                    
                    # Extract gene_name from attributes if available
                    gene_name_match = re.search(r'gene_name "([^"]*)"', attributes)
                    gene_name = gene_name_match.group(1) if gene_name_match and gene_name_match.group(1) else gene_id
                    
                    # If gene_name is empty, try to extract CG number
                    if not gene_name:
                        cg_number = extract_cg_number(gene_id)
                        if cg_number:
                            gene_name = cg_number
                    
                    # Extract CG number from gene_id or gene_name
                    cg_number = extract_cg_number(gene_id)
                    if not cg_number:
                        cg_number = extract_cg_number(gene_name)
                    
                    # Only add if we haven't processed this gene yet
                    if gene_id not in processed_genes:
                        processed_genes.add(gene_id)
                        gene_count += 1
                        
                        strand = fields[6]
                        
                        # Use gene_id as the primary identifier
                        gene_identifier = gene_id
                        
                        # If gene_id is a CG number, use that
                        if re.match(r'^CG\d+$', gene_id):
                            gene_identifier = gene_id
                        # Otherwise if gene_name is a CG number, use that
                        elif re.match(r'^CG\d+$', gene_name):
                            gene_identifier = gene_name
                        # Finally, use CG number if found
                        elif cg_number:
                            gene_identifier = cg_number
                        
                        genes.append({
                            'chrom': chromosome,  # Use the original chromosome name
                            'start': feature_start,
                            'end': feature_end,
                            'gene_id': gene_identifier,
                            'strand': strand,
                            # Store the CG number for WGCNA matching
                            '_cg_number': cg_number
                        })
        
        print(f"GTF stats: Found {gene_count} genes on chromosome {chromosome} in the selected region {start}-{end}")
        return genes
    except Exception as e:
        print(f"Error reading GTF file: {e}")
        return []

# Function to get average value from a bigwig for a given region
def get_bigwig_mean(bw_file, chrom, start, end):
    try:
        bw = pyBigWig.open(bw_file)
        
        # Try with original chromosome name
        if chrom in bw.chroms():
            use_chrom = chrom
        # Try with "chr" prefix
        elif "chr" + chrom in bw.chroms():
            use_chrom = "chr" + chrom
        # Try without "chr" prefix
        elif chrom.startswith("chr") and chrom[3:] in bw.chroms():
            use_chrom = chrom[3:]
        # Try RefSeq accession
        elif chrom in chrom_to_refseq and chrom_to_refseq[chrom] in bw.chroms():
            use_chrom = chrom_to_refseq[chrom]
        else:
            print(f"Warning: Chromosome {chrom} not found in {bw_file}")
            return None
        
        # Get the mean value for the region
        mean_val = bw.stats(use_chrom, start, end, type="mean")[0]
        bw.close()
        return mean_val if mean_val is not None else 0
    except Exception as e:
        print(f"Error reading {bw_file}: {e}")
        return None

# Extract genes and add data from bigwigs and WGCNA
def create_gene_data_table(chromosome, start, end, gene_to_wgcna, cg_to_wgcna, wgcna_columns):
    # Get genes from the GTF file
    genes = extract_genes_from_gtf(gtf_file, chromosome, start, end)
    
    if not genes:
        print(f"No genes found in {gtf_file} for region {chromosome}:{start}-{end}")
        return []
    
    # Add data from bigwigs for each gene
    results = []
    wgcna_match_count = 0
    
    for gene in genes:
        gene_start = gene['start']
        gene_end = gene['end']
        gene_id = gene['gene_id']
        cg_number = gene['_cg_number']
        
        # Add bigwig data
        gene['contact_wmel'] = get_bigwig_mean(bw_wmel, chromosome, gene_start, gene_end)
        gene['contact_dox'] = get_bigwig_mean(bw_dox, chromosome, gene_start, gene_end)
        gene['rnaseq_cov'] = get_bigwig_mean(rnaseq_cov, chromosome, gene_start, gene_end)
        gene['microc_cov'] = get_bigwig_mean(microc_cov, chromosome, gene_start, gene_end)
        gene['mappability'] = get_bigwig_mean(mappability, chromosome, gene_start, gene_end)
        
        # Try to match gene to WGCNA data
        wgcna_data = None
        
        # Try gene_id first
        if gene_id in gene_to_wgcna:
            wgcna_data = gene_to_wgcna[gene_id]
        # Then try CG number
        elif cg_number and cg_number in cg_to_wgcna:
            wgcna_data = cg_to_wgcna[cg_number]
        
        # Add all WGCNA data if available
        if wgcna_data:
            wgcna_match_count += 1
            for col in wgcna_columns:
                if col != wgcna_columns[0]:  # Skip the gene ID column
                    gene[f'wgcna_{col}'] = wgcna_data[col]
        else:
            # Add NA for all WGCNA columns if no match
            for col in wgcna_columns:
                if col != wgcna_columns[0]:  # Skip the gene ID column
                    gene[f'wgcna_{col}'] = 'NA'
        
        # Remove temporary fields used for matching
        gene.pop('_cg_number', None)
        
        results.append(gene)
    
    print(f"Found WGCNA matches for {wgcna_match_count} out of {len(genes)} genes")
    return results

# Main execution
print(f"Processing region {chrom}:{start_bp}-{end_bp}")

# Load WGCNA module assignments
gene_to_wgcna, cg_to_wgcna, wgcna_columns = load_wgcna_data(wgcna_modules_file)
print(f"Loaded WGCNA data with columns: {wgcna_columns}")

# Create gene data table
gene_data = create_gene_data_table(chrom, start_bp, end_bp, gene_to_wgcna, cg_to_wgcna, wgcna_columns)

# Save results
if gene_data:
    df = pd.DataFrame(gene_data)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved gene data to {output_file}")
    print(f"Total genes extracted: {len(gene_data)}")
else:
    print("No data to save - check your GTF file and region specification")