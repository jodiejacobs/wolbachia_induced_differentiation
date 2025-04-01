import pandas as pd
import subprocess
import os

# Step 1: Extract gene positions from GTF
def parse_gtf(gtf_file):
    # This will store gene_id -> (chrom, start, end)
    gene_positions = {}
    gene_coords = {}  # Temporary store to track min start and max end for each gene
    chrom_lengths = {}  # Track chromosome lengths
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:  # Skip incomplete lines
                continue
                
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            
            # Update chromosome length
            if chrom in chrom_lengths:
                chrom_lengths[chrom] = max(chrom_lengths[chrom], end)
            else:
                chrom_lengths[chrom] = end
            
            # Extract gene_id from attribute field
            attr = fields[8]
            gene_id = None
            for attr_pair in attr.split(';'):
                if 'gene_id' in attr_pair:
                    try:
                        gene_id = attr_pair.split('"')[1]
                        break
                    except IndexError:
                        continue
            
            if gene_id:
                # Track the gene's min start and max end positions
                if gene_id in gene_coords:
                    current_min, current_max, current_chrom = gene_coords[gene_id]
                    # Only update if same chromosome
                    if current_chrom == chrom:
                        gene_coords[gene_id] = (min(start, current_min), max(end, current_max), chrom)
                else:
                    gene_coords[gene_id] = (start, end, chrom)
    
    # Convert the min start and max end to final gene positions
    for gene_id, (start, end, chrom) in gene_coords.items():
        gene_positions[gene_id] = (chrom, start, end)
    
    return gene_positions, chrom_lengths

# Step 2: Match positions with module data
def create_bedgraph(gene_positions, module_file, output_bedgraph, chrom_lengths=None):
    # Read module data
    module_df = pd.read_csv(module_file, sep='\t')
    
    # Create dictionary to store points
    # Key: (chrom, position) to ensure uniqueness
    # Value: module score
    points_dict = {}
    
    for _, row in module_df.iterrows():
        gene_id = row['flybaseCG']
        module = row['module']
        
        if gene_id in gene_positions:
            chrom, start, end = gene_positions[gene_id]
            
            # Get module value
            try:
                score = int(module)
            except ValueError:
                # If module is categorical, map it to numbers
                module_dict = {f"module{i}": i for i in range(1, 100)}
                score = module_dict.get(module, 0)
            
            # Use midpoint as position
            midpoint = (start + end) // 2
            
            # Ensure midpoint is within chromosome bounds if chrom_lengths is provided
            if chrom_lengths and chrom in chrom_lengths:
                midpoint = min(midpoint, chrom_lengths[chrom] - 1)
            
            # Store in dictionary - if duplicate positions, keep highest module value
            key = (chrom, midpoint)
            if key in points_dict:
                points_dict[key] = max(points_dict[key], score)
            else:
                points_dict[key] = score
    
    # Convert dictionary to list of bedGraph entries
    bedgraph_data = []
    for (chrom, pos), score in points_dict.items():
        bedgraph_data.append((chrom, pos, pos+1, score))
    
    # Custom sort for Drosophila chromosomes
    def chrom_key(entry):
        chrom = entry[0]
        # Define chromosome order for standard Drosophila chromosomes
        chrom_order = {'2L': 1, '2R': 2, '3L': 3, '3R': 4, '4': 5, 'X': 6, 'Y': 7, 
                      'MT': 8, 'wMel': 9, 'UN8': 10, 'rRNA': 11}
        return (chrom_order.get(chrom, 999), entry[1])
    
    bedgraph_data.sort(key=chrom_key)
    
    # Write to bedGraph file
    with open(output_bedgraph, 'w') as f:
        for entry in bedgraph_data:
            f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")
    
    return output_bedgraph

# Function to create chromosome sizes file
def create_chrom_sizes(chrom_lengths, output_file):
    # Custom sort for Drosophila chromosomes
    def chrom_key(chrom_tuple):
        chrom = chrom_tuple[0]
        chrom_order = {'2L': 1, '2R': 2, '3L': 3, '3R': 4, '4': 5, 'X': 6, 'Y': 7, 
                      'MT': 8, 'wMel': 9, 'UN8': 10, 'rRNA': 11}
        return chrom_order.get(chrom, 999)
    
    with open(output_file, 'w') as f:
        for chrom, length in sorted(chrom_lengths.items(), key=chrom_key):
            f.write(f"{chrom}\t{length}\n")
    
    return output_file

# Step 3: Convert bedGraph to bigWig
def bedgraph_to_bigwig(bedgraph_file, chrom_sizes_file, output_bigwig):
    cmd = f"bedGraphToBigWig {bedgraph_file} {chrom_sizes_file} {output_bigwig}"
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        if result.stdout:
            print(f"bedGraphToBigWig output: {result.stdout}")
        return output_bigwig
    except subprocess.CalledProcessError as e:
        print(f"Error running bedGraphToBigWig: {e}")
        if hasattr(e, 'stdout') and e.stdout:
            print(f"stdout: {e.stdout}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"stderr: {e.stderr}")
        
        # Try with bedtools sort as fallback
        try:
            print("Trying to fix bedGraph with bedtools sort...")
            sorted_file = bedgraph_file + ".sorted"
            sort_cmd = f"bedtools sort -i {bedgraph_file} > {sorted_file}"
            subprocess.run(sort_cmd, shell=True, check=True)
            
            # Try again with sorted file
            retry_cmd = f"bedGraphToBigWig {sorted_file} {chrom_sizes_file} {output_bigwig}"
            subprocess.run(retry_cmd, shell=True, check=True)
            print(f"Successfully created bigWig using sorted bedGraph!")
            
            # Replace original with sorted
            os.rename(sorted_file, bedgraph_file)
            return output_bigwig
        except Exception as sort_error:
            print(f"Error with bedtools sort approach: {sort_error}")
            raise e

# Main workflow
def main(gtf_file, module_file, output_bigwig):
    # Chromosome mapping dictionary
    chromosomes = {
        'NC_004353.4': '4',
        'NC_004354.4': 'X',
        'NC_024511.2': 'MT',
        'NT_033777.3': '3R',
        'NT_033778.4': '2R',
        'NT_033779.5': '2L',
        'NT_037436.4': '3L',
        'NW_001845431.1': 'Y',
        'NW_007931083.1': 'UN8',
        'NW_007931121.1': 'rRNA',
        'NC_002978.6': 'wMel',
    }
    
    # Parse GTF to get gene positions and chromosome lengths
    gene_positions, chrom_lengths = parse_gtf(gtf_file)
    print(f"Extracted positions for {len(gene_positions)} genes")
    
    # For each gene position, convert the chromosome IDs to symbols
    gene_positions_with_symbols = {}
    for gene_id, (chrom, start, end) in gene_positions.items():
        chrom_symbol = chromosomes.get(chrom, chrom)
        gene_positions_with_symbols[gene_id] = (chrom_symbol, start, end)
    
    # Create chromosome sizes file with symbols
    chrom_lengths_with_symbols = {}
    for chrom, length in chrom_lengths.items():
        chrom_symbol = chromosomes.get(chrom, chrom)
        chrom_lengths_with_symbols[chrom_symbol] = length
    
    # Create bedGraph
    bedgraph_file = "temp_modules.bedgraph"
    create_bedgraph(gene_positions_with_symbols, module_file, bedgraph_file, chrom_lengths_with_symbols)
    print(f"Created bedGraph file: {bedgraph_file}")
    
    # Create chromosome sizes file
    chrom_sizes_file = "chrom.sizes"
    create_chrom_sizes(chrom_lengths_with_symbols, chrom_sizes_file)
    print(f"Created chromosome sizes file: {chrom_sizes_file}")
    
    # Print the first few lines of the bedGraph file
    print("First 5 lines of bedGraph file:")
    with open(bedgraph_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            print(line.strip())
    
    # Convert to bigWig
    try:
        bedgraph_to_bigwig(bedgraph_file, chrom_sizes_file, output_bigwig)
        print(f"Created bigWig file: {output_bigwig}")
    except Exception as e:
        print(f"Error creating bigWig file: {e}")
        print("Please check your bedGraph and chrom.sizes files manually")
    
    # Keep temporary files for debugging
    print("Keeping temporary files for debugging")

# Example usage
if __name__ == "__main__":
    gtf_file = "/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/micro-c/WGCNA_bigwig/Dmel6wMel2978.6.clean.gtf"
    module_file = "/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/micro-c/WGCNA_bigwig/significant_DE_WGCNA_modules.Dmel.tsv"
    output_bigwig = "/private/groups/russelllab/jodie/wolbachia_induced_DE/wolbachia_induced_differentiation/scripts/micro-c/WGCNA_bigwig/WGCNA_modules.bw"
    
    main(gtf_file, module_file, output_bigwig)