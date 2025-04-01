import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
import argparse
import pandas as pd

'''
From gtf file, genome file, and coverage file obtainded using samtools depth function, plots 
coverage across genome.
'''

parser = argparse.ArgumentParser()

parser.add_argument('--outputDir','-o',
                    default='/scratch1/jodie/wolbachia/Micro-C/26June2023_Micro_C/coverage/test/',
                    type=str,
                    action='store',
                    help='output file goes here')
parser.add_argument('--coverage','-i',
                    default='/scratch1/jodie/wolbachia/Micro-C/26June2023_Micro_C/bam/test.coverage.gz',
                    type=str,
                    action='store',
                    help='input file goes here')
parser.add_argument('--gtfFile','-gtf',
                    default='/scratch1/jodie/wolbachia/Dmelanogaster_wMel_RNAseq/reference_genomes/gtfs/Dmel6wMel2978.6.clean.gtf',
                    type=str,action='store',
                    help='input file goes here')
parser.add_argument('--genomeFile','-g',
                    default='/scratch1/jodie/wolbachia/Micro-C/reference_genomes/merged_ref.genome',
                    type=str,
                    action='store',
                    help='input file goes here')

args = parser.parse_args()

outDir=args.outputDir
gtf=args.gtfFile
genome=args.genomeFile
coverage = args.coverage

'''
Plot Histogram (bottom)
'''
def plotHisto(edge1, edge2, yhisto, panel):
    iBlue=(88/255,85/255,120/255)

    xhisto = np.arange(edge1,edge2).tolist()
    
    for i in range(len(xhisto)):
        rectangle = mplpatches.Rectangle([xhisto[i],0],1,yhisto[i],
                        facecolor=iBlue,
                        edgecolor=iBlue,
                        linewidth=.05
                        )
        panel.add_patch(rectangle)
        panel.set_ylim(0,max(yhisto)*1.1)
    
    panel.set_ylabel("Coverage")
    yticks = list(np.arange(0,max(yhisto), step = 100))
    yticks.append(max(yhisto)*1.1)
    panel.set_yticks(yticks)

    panel.set_xlabel("Position")
    xticks = list(np.arange(0,len(yhisto), step = 1000))
    xticks.append(len(yhisto))
    panel.set_xticks(xticks)

'''
Plot GTF
'''
def plotGTF(id, gtfDict, start, stop, y, panel):
    height = 0.05
    rectangle = mplpatches.Rectangle([start,y - height/2],stop-start,height,
                                facecolor='grey',
                                edgecolor='black',
                                linewidth=.5)
    panel.add_patch(rectangle)

    blocks = gtfDict[id]
    for block in blocks:
        start, stop, feature = block
        if feature == "CDS":
            height = 0.5
            rectangle = mplpatches.Rectangle([start,y - height/2],stop-start,height,
                                facecolor='grey',
                                edgecolor='black',
                                linewidth=.5)
            panel.add_patch(rectangle)
        else:
            height = 0.25
            rectangle = mplpatches.Rectangle([start,y - height/2],stop-start,height,
                    facecolor='grey',
                    edgecolor='black',
                    linewidth=.5)
            panel.add_patch(rectangle) 

def makeGTFpannel(gtf, chr, edge1, edge2, panel):
    gtfDict = {}

    with open(gtf,'r') as f:
        line=f.readline().split('\t')
        line=f.readline().split('\t')
        line=f.readline().split('\t')
        line=f.readline().split('\t')
        line=f.readline().split('\t')
        line=f.readline().split('\t')
        line=f.readline().split('\t')

        for line in f:
            line=line.split('\t')
            chromosome=line[0]
            start=int(line[3])
            end=int(line[4])
            feature=line[2]
            if chromosome == chr and start <= edge2 and end >= edge1:
                if feature in ['exon','CDS']:
                    metaData=line[8].split('transcript_id "')[1].split('"')[0]
                    if metaData not in gtfDict:
                        gtfDict[metaData] = []
                    gtfDict[metaData].append([start, end, feature])

    genes= []
    for key in gtfDict.keys():
        genes.append([min(gtfDict[key])[0], max(gtfDict[key])[1], key])

    y = -4
    while len(genes) > 0:
        genes.sort()
        x = 0 
        i = 0
        while i <= len(genes)-1:
            if genes[i][0] >= x:
                gene = genes.pop(i)
                x = gene[1]
                plotGTF(gene[2], gtfDict, gene[0], gene[1], y, panel)
                # top.add_patch(rectangle)
            else:
                i += 1
        y += 1 #blockheight

'''
Make and Save Figure
'''

#open the coverage file with pandas
cov = pd.read_csv(coverage, sep = '\t', header=None,compression='gzip', engine='python')

#For each line of the genome file, make a new plot
with open(genome, 'r') as g:
    for line in g:
        chr,edge2=line.split('\t')
        yhisto = cov[cov[0]==chr][2].to_list()
        edge1 = int(0)
        edge2 = len(yhisto)

        if edge2 > 0:
            
            out = outDir+chr+"_coverage.pdf"
            
            '''
            Setup Figure
            '''
            figureWidth=10
            figureHeight=6

            plt.style.use('/home/jodie/pythonScripts/BME163.mplstyle')
            plt.figure(figsize=(figureWidth,figureHeight))

            panelWidth=9
            panelHeight=3.75

            coveragePanel = plt.axes([.6/figureWidth,.375/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
            gtfPanel = plt.axes([.6/figureWidth,4.25/figureHeight,panelWidth/figureWidth,.25*panelHeight/figureHeight])
                                #[x values],[y values]

            for panel in [gtfPanel,coveragePanel]:
                panel.set_xlim(edge1,edge2)
                panel.set_ylim(-5,5)

                panel.tick_params(bottom=False, labelbottom=False,\
                                left=False, labelleft=False, \
                                right=False, labelright=False,\
                                top=False, labeltop=False)
            coveragePanel.tick_params(bottom = True, labelbottom=True,\
                                      left = True, labelleft=True)



            plotHisto(edge1, edge2, yhisto, coveragePanel)
            makeGTFpannel(gtf, chr, edge1, edge2, gtfPanel)
            plt.savefig(out, format = 'pdf',dpi=600)
            


