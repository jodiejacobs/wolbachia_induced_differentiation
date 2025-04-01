#This one works! Generates a distance matrix plot using abundance.txt files from kallisto updated 8/28/23 Jodie 
#Make sure to use the correct conda environment or else these packages will not work 
# library(biomaRt)
#https://stackoverflow.com/questions/76470802/cannot-install-r-bioconductor-packages-via-conda
#https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html
#https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html
library(tximport)
library(DESeq2)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)

#library(beeswarm)

## kallisto --> tximport --> DESeq2

#################################################
####### 1) Tximport
## import transcript-level estimates (for DESeq2) with tximport
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
## http://127.0.0.1:23850/library/tximport/doc/tximport.html#kallisto_with_TSV_files
## https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#tximport
##

## First, we locate the directory containing the files.
## create sample ids
kallisto_dir <- "/private/groups/russelllab/jodie/wolbachia_induced_DE/v00_kallisto"
base_dir <- "/private/groups/russelllab/jodie/wolbachia_induced_DE/v01_deseq2"
sample_id <- dir(kallisto_dir)

## create file list of directories
sample_dirs <- file.path(kallisto_dir, sample_id)

## save sample metadata for modelling
## A --> JW18 unif
## B --> JW18 wMel
## C --> S2 uninf
## D --> S2 wMel
metadata_path <- file.path(base_dir,"samples.txt")
metadata <- read.table(metadata_path, header = TRUE)
##reassign "condition" from character to factor
metadata$condition <- as.factor(metadata$condition)
# metadata$infection <- as.factor(metadata$infection)
# metadata$assay <- as.factor(metadata$assay)

## append directories to table in a new column; label column "path"
metadata <- dplyr::mutate(metadata, path = sample_dirs)

########################################
## Transcripts need to be associated with gene IDs for gene-level summarization.
## If that information is present in the files, we can skip this step. For Salmon, Sailfish, and kallisto the files only provide the transcript ID.
#####################################################
## We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID.
## The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.
####

ref_db <- file.path(base_dir, "Dmel6_cds.rna_from_genomic_gene_ids.tsv")
t2g <- read.table(ref_db, header=FALSE)

###############################################################################

### create a named vector pointing to the quantification files
files <- file.path(kallisto_dir, metadata$sample, "abundance.tsv")
names(files) <- metadata$sample
all(file.exists(files))

## kallisto with TSV files
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE)

#######################
## use the gene-level estimated counts from the quantification tools, and additionally to use the transcript-level abundance estimates to calculate a gene-level offset that corrects for changes to the average transcript length across samples
## the function DESeqDataSetFromTximport takes care of creation of the offset for you

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~celltype + infection + celltype:infection, ignoreRank=TRUE)
# dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~genotype )
##using counts and average transcript lengths from tximport
##Warning message:
##In DESeqDataSet(se, design = design, ignoreRank) :
##  some variables in design formula are characters, converting to factors

## The ddsTxi object here can then be used as dds in the following analysis steps.
#######################
## prefilter low-count/coverage genes across samples
## keep rows that have at least 10 reads total
#keep <- rowSums(counts(dds)) >= 10
## ensure at least X samples with a count of 10 or more, where X can be chosen as the sample size of the smallest group of samples (==2):
print(dds)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]

## normalization happens automatically, this command is not necessary (unless making a gct file)
dds <- estimateSizeFactors(dds)


## Fit statistical model
dds_model <- DESeq(dds)

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(dds_model, normalized = TRUE)+1)

# # add the gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=F)

# # write normalized counts to text file
write.table(norm_counts,"/private/groups/russelllab/jodie/wolbachia_induced_DE/v02_deseq2/normalized_counts.txt", quote=F, col.names=T, row.names=T, sep="\t")
