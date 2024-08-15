import os
import subprocess
import pandas as pd
import concurrent.futures

os.chdir('/scratch1/jodie/wolbachia/Micro-C/09Nov2023_Micro_C')
#High sensitivity, restrivtive p-value to 0.05 and sparcity to 0.88 *recomended by mustache 
files={
    'JW18-DOX-1_JW18-wMel-1':	['JW18-DOX-1',	'JW18-wMel-1'],
    'JW18-DOX-2_JW18-wMel-1':	['JW18-DOX-2',	'JW18-wMel-1'],
    'JW18-DOX-1_JW18-wMel-2':	['JW18-DOX-1',	'JW18-wMel-2'],
    'JW18-DOX-2_JW18-wMel-2':	['JW18-DOX-2',	'JW18-wMel-2'],
}


gtf_path = '/scratch1/jodie/wolbachia/Dmelanogaster_wMel_RNAseq/reference_genomes/gtfs/Drosophila_melanogaster.BDGP6.32.109.filtered.gtf'
feature_path = '/scratch1/jodie/wolbachia/Micro-C/reference_genomes/FlyBase_Annotations_all.txt'
output=list(files.keys())
file_dict = {'1kb':[], '16kb':[], '128kb':[]}
feature_list = ['ANNOTATION_SYMBOL', 'FEATURE_TYPE', 'GENE_SNAPSHOT', 'GO_BIOLOGICAL_PROCESS', 'GO_CELLULAR_COMPONENT', 'GO_MOLECULAR_FUNCTION', 'LOCATION_ARM', 'LOCATION_MAX', 'LOCATION_MIN', 'LOCATION_STRAND', 'NAME', 'SPECIES_ABBREVIATION', 'SYMBOL', 'UNIPROT_PROTEIN_FAMILY']

# Open feature_path as a pandas DataFrame with gene_id as the index
feature_df = pd.read_csv(feature_path, sep='\t')
# Create a dictionary with gene_id as the index
feature_df = feature_df.set_index('FBID_KEY')
feature_df = feature_df[~feature_df.index.duplicated(keep='first')]
gene_dict = feature_df.to_dict(orient='index')

def search_feature_dict(chromosome, start, end):
    filtered_genes = []
    for gene_id, gene_data in gene_dict.items():
        if gene_data['LOCATION_ARM'] == chromosome and int(gene_data['LOCATION_MIN']) > int(start) and int(gene_data['LOCATION_MAX']) < int(end):
            filtered_genes.append(gene_id)
    return filtered_genes

def create_gene_matrix(gene_ids):
    # Create an empty matrix with columns as the keys
    matrix = pd.DataFrame(columns=feature_list)

    # Iterate over the gene_ids
    for gene_id in gene_ids:
        # Get the gene data from the gene_dict
        gene_data = gene_dict[gene_id]
        # Add the gene data as a new row in the matrix
        matrix.loc[gene_id] = list(gene_data.values())

    # Return the matrix
    return matrix


#Run mustache at 1kb, 16kb, 128kb resolutions:
for out in output:
    file_dict['1kb'].extend([f'map/mustache_v4/1kb/{out}.diffloop1',f'map/mustache_v4/1kb/{out}.diffloop2'])
    file_dict['16kb'].extend([f'map/mustache_v4/16kb/{out}.diffloop1',f'map/mustache_v4/16kb/{out}.diffloop2'])
    file_dict['128kb'].extend([f'map/mustache_v4/128kb/{out}.diffloop1',f'map/mustache_v4/128kb/{out}.diffloop2'])

    map1 = f'map/{files[out][0]}.matrix_1kb.mcool'
    map2 = f'map/{files[out][1]}.matrix_1kb.mcool'
    
    for file in [map1,map2]:
        if os.path.isfile(file) == False:
            print(f"{file} not found" )
            break
    if os.path.isfile(f'map/mustache_v4/1kb/{out}.diffloop1') == False:
        bashCommand = f'python3 /home/jodie/mustache/mustache/diff_mustache.py -f1 {map1} -f2 {map2} -pt 0.01 -pt2 0.05 -o map/mustache_v4/1kb/{out} -r 1000 -st 0.88 -p 16'
        os.system(bashCommand)
    if os.path.isfile(f'map/mustache_v4/16kb/{out}.diffloop1') == False:
        bashCommand = f'python3 /home/jodie/mustache/mustache/diff_mustache.py -f1 {map1} -f2 {map2} -pt 0.01 -pt2 0.05 -o map/mustache_v4/16kb/{out} -r 16000 -st 0.88 -p 16'
        os.system(bashCommand)
    if os.path.isfile(f'map/mustache_v4/128kb/{out}.diffloop1') == False:
        bashCommand = f'python3 /home/jodie/mustache/mustache/diff_mustache.py -f1 {map1} -f2 {map2} -pt 0.01 -pt2 0.05 -o map/mustache_v4/128kb/{out} -r 128000 -st 0.88 -p 16'
        os.system(bashCommand)
        # print(bashCommand)

#Function creates bedfiles with chromosome, start pos, end pos, fdr 
def diffloop_to_bed(filepath):
    if os.path.isfile(f'{filepath}.bed') == False:
        bashCommand = f"tail -n +2 {filepath} | cut -f1,2,6,7 > {filepath}.bed"
        os.system(bashCommand)
    else:
        return

bedlist=[]
#Reformat bedfiles from mustache for use with bedtools
for range in file_dict.keys():
    for file in file_dict[range]:
        diffloop_to_bed(file)
    if os.path.isfile(f'map/mustache_v4/{range}/diffloop1.bed') == False:
        bashCommand1 = f"bedtools intersect -a {file_dict[range][0]}.bed -b {file_dict[range][2]}.bed {file_dict[range][4]}.bed {file_dict[range][6]}.bed > map/mustache_v4/{range}/diffloop1.bed"
        os.system(bashCommand)
    if os.path.isfile(f'map/mustache_v4/{range}/diffloop2.bed') == False:
        bashCommand2 = f"bedtools intersect -a {file_dict[range][1]}.bed -b {file_dict[range][3]}.bed {file_dict[range][5]}.bed {file_dict[range][7]}.bed > map/mustache_v4/{range}/diffloop2.bed"
        os.system(bashCommand2)
    bedlist.extend([f"map/mustache_v4/{range}/diffloop1.bed",f"map/mustache_v4/{range}/diffloop2.bed"])

#Function extracts a list of id's to use in PANGEA
def extract_unique_gene_ids(gtf_file, output_file):
    # Correctly format the awk command with proper escaping
    if os.path.isfile(output_file) == True:
        return
    bash_command = f"""awk '{{match($0, /gene_id "([^"]+)"/, arr); if (arr[1] != "") print arr[1]}}' {gtf_file} | sort | uniq > {output_file}"""
    process = subprocess.run(bash_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(bash_command)
    if process.returncode == 0:
        print(f"Unique gene IDs have been exported to {output_file}")
    else:
        print(f"Command failed with error: {process.stderr.decode('utf-8')}")

def extract_annotation(bed_file):
    try:
        bed = pd.read_csv(bed_file, sep='\t', header=None)
        bed.columns = ['chromosome', 'start_pos', 'end_pos', 'fdr']
    except Exception as e:
        print(f'Error:{e}')
        return

    # Initialize an empty pandas DataFrame
    df = pd.DataFrame(columns=['chromosome', 'start_pos', 'end_pos', 'fdr', 'gene_id'] + feature_list)
    
    # Iterate over the bed file
    for index, row in bed.iterrows():
        # Access the values in each row
        chromosome = row['chromosome']
        start_pos = row['start_pos']
        end_pos = row['end_pos']
        fdr = row['fdr']

        # Filter the gene_dict for genes within the specified range
        gene_ids = search_feature_dict(chromosome, start_pos, end_pos)  
        # Create a new DataFrame with the gene_ids  
        gene_matrix = create_gene_matrix(gene_ids) 
        # Add the chromosome, start_pos, end_pos, and fdr to the DataFrame  
        gene_matrix['chromosome'] = chromosome  
        gene_matrix['start_pos'] = start_pos
        gene_matrix['end_pos'] = end_pos
        gene_matrix['fdr'] = fdr
        # Add the gene_id to the DataFrame
        gene_matrix['gene_id'] = gene_ids
        # Concatenate the new DataFrame with the existing DataFrame
        df = pd.concat([df, gene_matrix], axis=0, ignore_index=True)
        # print(df)
    # Export the DataFrame to a CSV file
    df.to_csv(f"{bed_file}.annotation", sep='\t', header=df.columns, index=False)
    print(f"Annotation has been exported to {bed_file}.annotation")


def process_bed(bed):
    bashCommand = f"bedtools intersect -a {gtf_path} -b {bed} > {bed}.gtf"
    os.system(bashCommand)
    print(f"Annotation has been exported to {bed}.gtf")
    extract_annotation(bed)
    extract_unique_gene_ids(f"{bed}.gtf", f"{bed}.gene_ids")

    

# Create a ThreadPoolExecutor with the desired number of threads
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    # Submit each bed file to the executor for processing
    futures = [executor.submit(process_bed, bed) for bed in bedlist]
    
    # Wait for all the futures to complete
    concurrent.futures.wait(futures)

#Concatenate the annotation files into one file for each diffloop (1 and 2)
loop1=[bedlist[0],bedlist[2],bedlist[4]]   
loop2=[bedlist[1],bedlist[3],bedlist[5]]


bashCommand = f"cat {' '.join([f'{bed}.annotation' for bed in loop1])} > map/mustache_v4/diffloop1.annotation.tsv"
os.system(bashCommand)
bashCommand = f"cat {' '.join([f'{bed}.annotation' for bed in loop2])} > map/mustache_v4/diffloop2.annotation.tsv"
os.system(bashCommand)
