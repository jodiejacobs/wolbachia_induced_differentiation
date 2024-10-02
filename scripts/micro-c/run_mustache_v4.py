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

# # Open feature_path as a pandas DataFrame with gene_id as the index
feature_df = pd.read_csv(feature_path, sep='\t')
feature_df.drop(columns=['LOCATION_ARM', 'LOCATION_MAX', 'LOCATION_MIN', 'LOCATION_STRAND', 'NAME'], inplace=True)
# feature_dict = feature_df.set_index('FBID_KEY').T.to_dict()

# def search_feature_dict(gene_id):
#     if gene_id in feature_dict:
#         return list(feature_dict[gene_id].values())
#     else:
#         return []

# gene_data = search_feature_dict('FBgn0267928')
# print(gene_data)
##['CR46209', 'SO0002127:lncRNA_gene', '-', '-', '-', '-', '2L', '19215967', '19215451', '1', 'long non-coding RNA:CR46209', 'Dmel', 'lncRNA:CR46209', '-']

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
    bashCommand = f"tail -n +2 {filepath} | cut -f1,2,6,7 > {filepath}.bed"
    os.system(bashCommand)

bedlist=[]
#Reformat bedfiles from mustache for use with bedtools
for range in file_dict.keys():
    for file in file_dict[range]:
        diffloop_to_bed(file)
    bashCommand1 = f"bedtools intersect -a {file_dict[range][0]}.bed -b {file_dict[range][2]}.bed {file_dict[range][4]}.bed {file_dict[range][6]}.bed | uniq > map/mustache_v4/{range}/diffloop1.bed"
    bashCommand2 = f"bedtools intersect -a {file_dict[range][1]}.bed -b {file_dict[range][3]}.bed {file_dict[range][5]}.bed {file_dict[range][7]}.bed | uniq > map/mustache_v4/{range}/diffloop2.bed"
    bedlist.extend([f"map/mustache_v4/{range}/diffloop1.bed",f"map/mustache_v4/{range}/diffloop2.bed"])
    os.system(bashCommand1)
    os.system(bashCommand2)

# Function to parse the attributes field
def parse_attributes(attributes):
    attr_dict = {}
    attributes = attributes.strip(';').split('; ')
    for attribute in attributes:
        key, value = attribute.split(' ')
        attr_dict[key] = value.strip('"')
    return attr_dict

#Function extracts a list of id's to use in PANGEA
def extract_unique_gene_ids(gtf_file, output_file):
    # Correctly format the awk command with proper escaping
    bash_command = f"""awk '{{match($0, /gene_id "([^"]+)"/, arr); if (arr[1] != "") print arr[1]}}' {gtf_file} | sort | uniq > {output_file}"""
    process = subprocess.run(bash_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(bash_command)
    if process.returncode == 0:
        print(f"Unique gene IDs have been exported to {output_file}")
    else:
        print(f"Command failed with error: {process.stderr.decode('utf-8')}")

def extract_annotation(gtf_file, bed_file):
    gtf = pd.read_csv(gtf_file, sep='\t', header=None)
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed.columns = ['chromosome', 'start_pos', 'end_pos', 'fdr']
    gtf.columns=['chromosome',  'source',   'feature',  'start_pos',    'end_pos',  'score',    'strand',   'frame',    'attribute']
                #3R             FlyBase     gene        7145880         7150968     .           -           .           gene_id "FBgn0250732"; gene_name "gfzf"; gene_source "FlyBase"; gene_biotype "protein_coding";
    
    # Parse the attributes field and expand it into separate columns
    attributes_df = gtf['attribute'].apply(parse_attributes).apply(pd.Series)

    # Concatenate the original DataFrame with the new attributes DataFrame
    gtf = pd.concat([gtf.drop(columns=['attribute']), attributes_df], axis=1)
    
    # #print(gtf.columns): 
    # # Index(['chromosome', 'source', 'feature', 'start_pos', 'end_pos', 'score',
    # #    'strand', 'frame', 'gene_id', 'gene_name', 'gene_source',
    # #    'gene_biotype', 'transcript_id', 'transcript_name', 'transcript_source',
    # #    'transcript_biotype', 'exon_number', 'exon_id', 'tag', 'protein_id'],
    # #   dtype='object')
    
    # Initialize an empty pandas DataFrame
    df = pd.DataFrame(columns=bed.columns.tolist() + gtf.columns.tolist())

    for index, row in bed.iterrows():
        # Access the values in each row
        chromosome = row['chromosome']
        start_pos = row['start_pos']
        end_pos = row['end_pos']
        fdr = row['fdr']
                
        gtf_filtered = gtf[(gtf['chromosome'] == chromosome) & (gtf['start_pos'] > start_pos) & (gtf['end_pos'] < end_pos)]
        
        for index, entry in gtf_filtered.iterrows():
            # gene = entry['gene_id']
            # features = feature_dict[gene]
            # Create a new row with the additional columns
            new_row = pd.Series([chromosome, start_pos, end_pos, fdr] + entry.tolist(), index=df.columns)
            # Append the new row to the DataFrame
            df.loc[len(df.index)] = new_row
            # Save the DataFrame to a new file

    df.to_csv(f"{bed_file}.annotation", sep='\t', header=None, index=False)
    print(f"Annotation has been exported to {bed_file}.annotation")

#Filter gtfs by bedfiles and 
# Filter gtfs by bedfiles

def process_bed(bed):
    bashCommand = f"bedtools intersect -a {gtf_path} -b {bed} > {bed}.gtf"
    os.system(bashCommand)
    extract_annotation(f"{bed}.gtf", bed)
    extract_unique_gene_ids(f"{bed}.gtf", f"{bed}.gene_ids")

# for bed in bedlist:
#     process_bed(bed)
print(bedlist)
# Create a ThreadPoolExecutor with the desired number of threads
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    # Submit each bed file to the executor for processing
    futures = [executor.submit(process_bed, bed) for bed in bedlist]
    
    # Wait for all the futures to complete
    concurrent.futures.wait(futures)
