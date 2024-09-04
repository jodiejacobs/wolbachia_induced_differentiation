

import os
import sys 
import h5py    
import numpy as np    
import argparse
import pandas as pd


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser()

parser.add_argument('-i','--inFile', 
                    default = '/scratch1/jodie/wolbachia/Dmelanogaster_wMel_RNAseq/01_09_23_Azenta/GeneExpressionStudy/kallisto_output/JW18wMel221117-1/abundance.h5',
                    type=str,
                    action = 'store',
                    help='path to input directory')
parser.add_argument('-o','--outDir',
                    default='/scratch1/jodie/wolbachia/Dmelanogaster_wMel_RNAseq/01_09_23_Azenta/GeneExpressionStudy/kallisto_output/JW18wMel221117-1/',
                    type=dir_path,
                    action='store',
                    help='output dir here')
parser.add_argument('-w','--Wolbachia',
                    default='H5_WP.h5',
                    type=str,
                    action='store',
                    help='output dir here')
parser.add_argument('-d','--Drosophila',
                    default='H5_DM.h5',
                    type=str,
                    action='store',
                    help='output dir here')
parser.add_argument('-s','--save',
                    default='True',
                    type=bool,
                    action='store',
                    help='Save Files Bool Here')

args = parser.parse_args()

inFile = args.inFile
wpOut = args.Wolbachia
dmOut = args.Drosophila 
outDir = args.outDir
saveBool = args.save

SAVE = saveBool  # whether to save the outputs of this notebook

if SAVE:
    print(f'SAVE={saveBool}')
    OUTPUT_PATH = outDir #Directory where exported files are saved
    os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist

def filter_counts(est_counts, id, sp):
    df = pd.DataFrame(
            data = est_counts,
            index = id
        )
    if sp == "w":
        df_sp = df[(df.index.str.contains(str("WD"))|df.index.str.contains(str("WP")))]
    else:
        df_sp = df[~(df.index.str.contains(str("WD"))|df.index.str.contains(str("WP")))]
    # print(f"Length of OG Dataframe:{len(df)},\n Number of Dmel: {len(df_dm)},\n Number of WD:{len(df_wd)}")
    
    return df_sp

def create_h5(oldFile, ids , newFile, species):
    df = filter_counts(f['est_counts'], ids, species)
    layer1 = list(oldFile.keys())
    for l in layer1:
        path = "/"+l
        if l == 'est_counts':
            newFile.create_dataset(path,data=df[0])
        else:
            layer2 = list(oldFile[l])
            for l2 in layer2:
                path = '/'+l+'/'+ l2 
                if l2 == "ids":
                    newFile.create_dataset(path,data=list(df.index), dtype = "|S50")
                elif l2 == "eff_lengths":
                    ef = filter_counts(f[l][l2], ids, species)
                    newFile.create_dataset(path,data=ef[0])
                elif l2 == "lengths":
                    ef = filter_counts(f[l][l2], ids, species)
                    newFile.create_dataset(path,data=ef[0], dtype = "<i4")
                elif l2[:2] == "bs":
                    bs = filter_counts(f[l][l2], ids, species)
                    newFile.create_dataset(path,data=bs[0])
                else:
                    newFile.create_dataset(path,data=oldFile[l][l2])

def check_layers(H5):
    for l1 in H5.keys():
        if l1 == 'est_counts':
            print(H5[l1])  
        else: 
            layer2 = list(H5[l1])
            for l2 in layer2:
                print(l1, H5[l1][l2])


with h5py.File(inFile,'r') as f:

    id = list(f["aux"]["ids"])
    id = [str(i) for i in id]

    print("H5:")
    check_layers(f)

    with h5py.File(outDir+wpOut, 'w') as WP:
        create_h5(oldFile = f, ids = id, newFile = WP, species ="w")  
        print("WP:")
        check_layers(WP)

    with h5py.File(outDir+dmOut, 'w') as DM:
        create_h5(oldFile = f,ids = id, newFile = DM, species = "d")
        print("DM:")
        check_layers(DM)


       

    



