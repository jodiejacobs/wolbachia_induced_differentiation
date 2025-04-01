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
kallisto_dir <- "/private/groups/russelllab/jodie/wolbachia_induced_DE/v04_deseq2/kallisto_output"
base_dir <- "/private/groups/russelllab/jodie/wolbachia_induced_DE/v04_deseq2"
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

#ref_db <- file.path(base_dir,e/groups/russelllab/jod"Dmel6_cds.rna_from_genomic_gene_ids.tsv")
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

dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData = metadata, design = ~celltype + infection + celltype*infection, ignoreRank=TRUE)

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
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

## normalization happens automatically, this command is not necessary (unless making a gct file)
# dds <- estimateSizeFactors(dds)


# Fit statistical model
dds_model <- DESeq(dds)


# ###############################
# ### Wald test ################
# ###############################
# ###############################

dds$group <- factor(paste0(dds$celltype, dds$infection))

# ################################################
# #####################################
# #################################################################
# #################################################################
# ###################### MODELING INTERACTIONS: ###################
# #################################################################
# #################################################################
# ########## WALD formatting: retain the terms to test ############
# #################################################################
# ## interactions modeled - Wald test
design(dds) <- ~ celltype + infection + celltype*infection 

dds <- DESeq(dds)
resultsNames(dds)


# > resultsNames(dds)
# [1] "Intercept"                "celltype_S2_vs_JW18"     
# [3] "infection_wMel_vs_DOX"    "celltypeS2.infectionwMel"

dds$celltype <- relevel(dds$celltype, "S2")
dds <- DESeq(dds)
resultsNames(dds)

# > resultsNames(dds)
# [1] "Intercept"                  "celltype_JW18_vs_S2"       
# [3] "infection_wMel_vs_DOX"      "celltypeJW18.infectionwMel"

S2 <- results(dds, (list(c("celltypeS2.infectionwMel"))))
JW18 <- results(dds, (list(c("celltypeJW18.infectionwMel"))))

# JW18_filtered_df <- as.data.frame(JW18) %>% filter(., padj<0.001)# %>% filter(.,log2FoldChange>1)
# > length(JW18_filtered_df$padj)
# # # [1] 132

S2_filtered_df <- as.data.frame(S2) %>% filter(., padj<0.001) #%>% filter(.,log2FoldChange>5)
# > length(S2_filtered_df$padj)
# # # [1] 166

# Create a function to test the effect of each interaction 
wald_test <- function(interaction){
  ## infection effect, controlling for interaction 

  dataset <- results(dds, contrast=list( c(interaction)))

  pdf_path  <- file.path(base_dir, "plots", paste0("Kallisto-DESeq2_WaldTest_",interaction,".pdf"))
  csv_path  <- file.path(base_dir, "plots", paste0("Kallisto-DESeq2_WaldTest_",interaction,".csv"))
  vol_path  <- file.path(base_dir, "plots", paste0("Kallisto-DESeq2_WaldTest_Volcano_",interaction,".pdf"))

  write.csv(dataset, file=(csv_path))

  ## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
  ## Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
  pdf(pdf_path)
  plotMA(dataset, ylim=c(-2,2))
  dev.off()
}

wald_plots =  c("Intercept")
for (i in resultsNames(dds)){
  if (!(i %in% wald_plots)) {
    wald_test(i)
    wald_plots = append(wald_plots, i)
  }
}

# dds$infection <- relevel(dds$infection, "DOX")
dds$celltype <- relevel(dds$celltype, "JW18")

dds <- DESeq(dds)
resultsNames(dds)

for (i in resultsNames(dds)){
  if (!(i %in% wald_plots)) {
    wald_test(i)
    wald_plots = append(wald_plots, i)
  }
}

######################################################
######## Volcano Plots from DESeq2 #######
######################################################

wald_path <- file.path(base_dir, "plots")

#[1] Kallisto-DESeq2_WaldTest_celltype_S2_vs_JW18.csv
csv <- "Kallisto-DESeq2_WaldTest_celltype_S2_vs_JW18.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

#[2] Kallisto-DESeq2_WaldTest_celltypeS2.infectionwMel.csv
csv <- "Kallisto-DESeq2_WaldTest_celltypeS2.infectionwMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

#[3] Kallisto-DESeq2_WaldTest_infection_wMel_vs_DOX.csv
csv <- "Kallisto-DESeq2_WaldTest_infection_wMel_vs_DOX.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()

# [4] Kallisto-DESeq2_WaldTest_celltypeJW18.infectionwMel.csv
csv <- "Kallisto-DESeq2_WaldTest_celltypeJW18.infectionwMel.csv"
pval_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-pval_volcano.pdf"))
padj_volcano <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-padj_volcano.pdf"))
histogram <- file.path(wald_path, paste0(tools::file_path_sans_ext(csv),"-histogram.pdf"))

csv_path <-file.path(wald_path, csv)
data <- read.csv(csv_path)

pdf(pval_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='pvalue', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(padj_volcano)
EnhancedVolcano(data, lab=data$X, x='log2FoldChange', y='padj', pCutoff = 1e-03, col = c("grey30", "cyan3", "skyblue", "orangered1"))
dev.off()

pdf(histogram)
hist(data$log2FoldChange)
dev.off()



# Specify the contrast to test for the effect of infection (wMel) in both cell types
contrast <- c("celltypeS2.infectionwMel", "celltypeJW18.infectionwMel")
res <- results(dds, contrast = c(JW18, S2))

# Filter genes that are significantly differentially expressed in both cell types
enriched_genes <- rownames(subset(res, padj < 0.001))


genes_to_plot <- c("Dmel_CG33946","Dmel_CG14688", "Dmel_CG15571", "Dmel_CG13972", "Dmel_CG5402","Dmel_CG13299","Dmel_CG33943","Dmel_CG1780","Dmel_CG11742","Dmel_CG30268","Dmel_CG5076","Dmel_CG13230","Dmel_CG17330,""Dmel_CG4356","Dmel_CG30275","Dmel_CG43895","Dmel_CR46075","Dmel_CG11037","Dmel_CG10062","Dmel_CG12885")

# ##################################################################################
## infection effect, controlling for genotype
#outputI <- results(dds, contrast=list( c("infection_wMel_vs_DOX", "genotypeOreR.infectionwMel")))

# outputI <- results(dds, contrast=list( c("infection_wMel_vs_DOX")))
# write.csv(outputI, file="/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-infection_minusGeno.1.csv")

# ## In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet.
# ## Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
# pdf("/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-infection_minusGeno_MAplot.pdf")
# plotMA(outputI, ylim=c(-2,2))
# dev.off()

# # Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes.
# outputIlfc <- lfcShrink(dds, coef="infection_wMel_vs_DOX", type="apeglm")

# # > summary(outputIlfc)

# # out of 12244 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 3497, 29%
# # LFC < 0 (down)     : 2740, 22%
# # outliers [1]       : 4, 0.033%
# # low counts [2]     : 0, 0%
# # (mean count < 0)
# # [1] see 'cooksCutoff' argument of ?results
# # [2] see 'independentFiltering' argument of ?results

# ###############################
# ## interactions modeled - Wald test
# ## genotype effect, controlling for infection
# #outputG <- results(dds, contrast=list( c("genotype_OreR_vs_meiP261", "genotypeOreR.infectionwMel")))

# outputG <- results(dds, contrast=list( c("celltype_S2_vs_JW18")))
# write.csv(outputG, file="/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-celltype_minusInf.1.csv")

# pdf("/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-celltype_minusInf_MAplot.pdf")
# plotMA(outputG, ylim=c(-2,2))
# dev.off()

# # > summary(outputG)

# # out of 11080 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 3949, 36%
# # LFC < 0 (down)     : 4408, 40%
# # outliers [1]       : 2, 0.018%
# # low counts [2]     : 0, 0%
# # (mean count < 1)
# # [1] see 'cooksCutoff' argument of ?results
# # [2] see 'independentFiltering' argument of ?results

# #################################
# ## interaction between genotype and infection
# outputGI <- results(dds, name="celltypeS2.infectionwMel")
# write.csv(outputGI, file="/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-GxIinteraction.csv")

# pdf("/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-GxIinteraction_MAplot.pdf")
# plotMA(outputGI, ylim=c(-2,2))
# dev.off()

# # > summary(outputGI)

# # out of 11080 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)       : 3430, 31%
# # LFC < 0 (down)     : 2885, 26%
# # outliers [1]       : 2, 0.018%
# # low counts [2]     : 0, 0%
# # (mean count < 1)
# # [1] see 'cooksCutoff' argument of ?results
# # [2] see 'independentFiltering' argument of ?results


# #################################
# ## interaction between genotype and infection
# outputA1 <- results(dds, name="assay_ribodepletRNAseq_vs_pseudorandom_hexamer")
# write.csv(outputGI, file="/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-assay_ribodepletRNAseq_vs_pseudorandom_hexamer.csv")

# pdf("/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-assay_ribodepletRNAseq_vs_pseudorandom_hexamer_MAplot.pdf")
# plotMA(outputA1, ylim=c(-2,2))
# dev.off()

# > summary(outputA1)

# out of 11080 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2250, 20%
# LFC < 0 (down)     : 4175, 38%
# outliers [1]       : 2, 0.018%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# #################################
# ## interaction between genotype and infection
# outputA2 <- results(dds, name="assay_rRNA.degredation_vs_pseudorandom_hexamer")
# write.csv(outputGI, file="/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-assay_rRNA.degredation_vs_pseudorandom_hexamer.csv")

# pdf("/private/groups/russelllab/jodie/rRNA_degredation/v02_deseq2/plots/Kallisto-DESeq2_WaldTest-assay_rRNA.degredation_vs_pseudorandom_hexamer.pdf")
# plotMA(outputA1, ylim=c(-2,2))
# dev.off()

# > summary(outputA2)

# out of 11080 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2414, 22%
# LFC < 0 (down)     : 3827, 35%
# outliers [1]       : 2, 0.018%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results




# ##############################################
# ###########################################
# ################ LRTs #######################
# ###############################################
# ##########################################
# ## The LRT test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero
# ## The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables
# ## in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model
# ## The likelihood ratio test can be performed by specifying test="LRT" when using the DESeq function, and providing a reduced design formula, e.g. one in which a number of terms from design(dds) are removed.

# #################################################################
# ########## LRT formatting: leave out the term to test ###########
# #################################################################
# ## interactions modeled - LRT
# ## full model = design
# ## only "reduced" is recognized

# design = ~ celltype + infection + assay + celltype*infection

# ## impact of infection + assay + celltype:infection
# dds <- DESeq(dds, test="LRT", full=design, reduced= ~celltype)
# LRT_IGxI <- results(dds)
# write.csv(LRT_A_IGxI, file="/Users/jodiejacobs/Documents/wolbachia/celltype_mapping/Kallisto-DESeq2_LRT-I+GxI.1.csv")

# ## impact of genotype + genotype:infection
# dds <- DESeq(dds, test="LRT", full=design, reduced= ~infection)
# LRT_GGxI <- results(dds)
# write.csv(LRT_GGxI, file="/Users/jodiejacobs/Documents/wolbachia/celltype_mapping/Kallisto-DESeq2_LRT-G+GxI.csv")

# ## interaction effect
# dds <- DESeq(dds, test="LRT", full=design, reduced= ~genotype + infection)
# LRT_GxI <- results(dds)
# write.csv(LRT_GxI, file="/Users/jodiejacobs/Documents/wolbachia/celltype_mapping/Kallisto-DESeq2_LRT-GxI.csv")

# ## genotype effect
# dds <- DESeq(dds, test="LRT", full=design, reduced= ~infection + genotype:infection)
# LRT_G <- results(dds)
# write.csv(LRT_G, file="/Users/jodiejacobs/Documents/wolbachia/celltype_mapping/Kallisto-DESeq2_LRT-G.csv")
# ###Error in nbinomLRT(object, full = full, reduced = reduced, quiet = quiet,  :
# ###  less than one degree of freedom, perhaps full and reduced models are not in the correct order

# ## infection effect
# dds <- DESeq(dds, test="LRT", full=design, reduced= ~genotype + genotype:infection)
# LRT_I <- results(dds)
# write.csv(LRT_I, file="/Users/jodiejacobs/Documents/wolbachia/celltype_mapping/Kallisto-DESeq2_LRT-I.csv")
# ## Error in nbinomLRT(object, full = full, reduced = reduced, quiet = quiet,  :
# ##  less than one degree of freedom, perhaps full and reduced models are not in the correct order
# # ==> see "deseq2_modelling_notes.txt for more information"

# ####################################
# ####################################
# ####################################

# pdf("Kallisto-DESeq2_Dmel_genes_LRT-geno-baseMean.histo.pdf")
# hist(LRT_I$baseMean)
# dev.off()

# pdf("Kallisto-DESeq2_Dmel_genes_LRT-geno-baseMean.log10histo.pdf")
# hist(log10(LRT_I$baseMean))
# dev.off()

# pdf("Kallisto-DESeq2_Dmel_genes_LRT-geno-log2FoldChange.histo.pdf")
# hist(LRT_I$log2FoldChange)
# dev.off()

# pdf("Kallisto-DESeq2_Dmel_genes_LRT-geno-log2FoldChange.log10histo.pdf")
# hist(log10(LRT_I$log2FoldChange))
# dev.off()

# ####################################
# ####################################
# ####### plot gene counts ###########
# ####################################
# ####################################

# pdf("KallistoDESeq2_min_SigGeno_gene.pdf")
# plotCounts(dds, gene=which.min(outputG$padj), intgroup="condition")
# dev.off()

# pdf("KallistoDESeq2_min_SigInfgene.pdf")
# plotCounts(dds, gene=which.min(outputI$padj), intgroup="condition")
# dev.off()

# pdf("KallistoDESeq2_min_SigIGInteract_gene.pdf")
# plotCounts(dds, gene=which.min(outputGI$padj), intgroup="condition")
# dev.off()


# ##########################################
# ##########################################
# pdf("KallistoDESeq2_Dmel_CG12218_mei-P26_counts.pdf")
# plotCounts(dds, gene="Dmel_CG12218", intgroup="condition")
# dev.off()

# ## returnData: should the function only return the data.frame of counts and covariates for custom plotting (default is FALSE)
# #######################
# ## geno (-infection) ##
# #######################
# ## assign "group"
# dds$group <- factor(paste0(dds$genotype, dds$infection))

# meiP26 <- plotCounts(dds, gene="Dmel_CG12218", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG12218_mei-P26_counts.pdf")
# ggplot(data=meiP26, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG12218_mei-P26_log2.counts.pdf")
# ggplot(data=meiP26, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# orb <- plotCounts(dds, gene="Dmel_CG10868", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG10868_orb_counts.pdf")
# ggplot(data=orb, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG10868_orb_log2.counts.pdf")
# ggplot(data=orb, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# wuho <- plotCounts(dds, gene="Dmel_CG15897", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG15897_wuho_counts.pdf")
# ggplot(data=wuho, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG15897_log2.wuho_counts.pdf")
# ggplot(data=wuho, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()



# ###############################
# ###### infection (-geno) ######
# ###############################
# Dmel_CG32834 <- plotCounts(dds, gene="Dmel_CG32834", intgroup="group", returnData=TRUE)

# OreRuninfected <- subset(Dmel_CG32834, group == "OreRuninfected")
# OreRwMel <- subset(Dmel_CG32834, group == "OreRwMel")
# meiP261uninfected <- subset(Dmel_CG32834, group == "meiP261uninfected")
# meiP261wMel <- subset(Dmel_CG32834, group == "meiP261wMel")

# wilcox.test(OreRuninfected$count, OreRwMel$count)
# wilcox.test(meiP261uninfected$count, meiP261wMel$count)
# wilcox.test(meiP261uninfected$count, OreRuninfected$count)
# wilcox.test(OreRwMel$count, meiP261wMel$count)
# wilcox.test(meiP261wMel$count, OreRuninfected$count)
# wilcox.test(OreRwMel$count, meiP261uninfected$count)

# > wilcox.test(OreRuninfected$count, OreRwMel$count)

# 	Wilcoxon rank sum test with continuity correction

# data:  OreRuninfected$count and OreRwMel$count
# W = 21, p-value = 0.4047
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(OreRuninfected$count, OreRwMel$count) :
#   cannot compute exact p-value with ties
# > wilcox.test(meiP261uninfected$count, meiP261wMel$count)

# 	Wilcoxon rank sum exact test

# data:  meiP261uninfected$count and meiP261wMel$count
# W = 28, p-value = 0.132
# alternative hypothesis: true location shift is not equal to 0

# > wilcox.test(meiP261uninfected$count, OreRuninfected$count)

# 	Wilcoxon rank sum test with continuity correction

# data:  meiP261uninfected$count and OreRuninfected$count
# W = 30, p-value = 0.0562
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(meiP261uninfected$count, OreRuninfected$count) :
#   cannot compute exact p-value with ties
# > wilcox.test(OreRwMel$count, meiP261wMel$count)

# 	Wilcoxon rank sum test with continuity correction

# data:  OreRwMel$count and meiP261wMel$count
# W = 3, p-value = 0.009622
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(OreRwMel$count, meiP261wMel$count) :
#   cannot compute exact p-value with ties
# > wilcox.test(meiP261wMel$count, OreRuninfected$count)

# 	Wilcoxon rank sum test with continuity correction

# data:  meiP261wMel$count and OreRuninfected$count
# W = 27.5, p-value = 0.124
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(meiP261wMel$count, OreRuninfected$count) :
#   cannot compute exact p-value with ties
# > wilcox.test(OreRwMel$count, meiP261uninfected$count)

# 	Wilcoxon rank sum test with continuity correction

# data:  OreRwMel$count and meiP261uninfected$count
# W = 0, p-value = 0.002778
# alternative hypothesis: true location shift is not equal to 0

# Warning message:
# In wilcox.test.default(OreRwMel$count, meiP261uninfected$count) :
#   cannot compute exact p-value with ties

# #########################
# Dmel_CG14606 <- plotCounts(dds, gene="Dmel_CG14606", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG14606_transmemTransport_counts.pdf")
# ggplot(data=Dmel_CG14606, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG14606_transmemTransport_log2.counts.pdf")
# ggplot(data=Dmel_CG14606, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# aEst2 <- plotCounts(dds, gene="Dmel_CG2505", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG2505_α-Est2_counts.pdf")
# ggplot(data=aEst2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG2505_α-Est2_log2.counts.pdf")
# ggplot(data=aEst2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CG17167 <- plotCounts(dds, gene="Dmel_CG17167", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG17167_KNaAntiporter_counts.pdf")
# ggplot(data=Dmel_CG17167, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG17167_KNaAntiporter_log2.counts.pdf")
# ggplot(data=Dmel_CG17167, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# stlk <- plotCounts(dds, gene="Dmel_CG40293", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG40293_stlk_counts.pdf")
# ggplot(data=stlk, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG40293_stlk_log2.counts.pdf")
# ggplot(data=stlk, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CG13871 <- plotCounts(dds, gene="Dmel_CG13871", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG13871_hypProt_counts.pdf")
# ggplot(data=Dmel_CG13871, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG13871_hypProt_log2.counts.pdf")
# ggplot(data=Dmel_CG13871, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CG42335 <- plotCounts(dds, gene="Dmel_CG42335", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG42335_metalloaminopeptidase_counts.pdf")
# ggplot(data=Dmel_CG42335, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG42335_metalloaminopeptidase_log2.counts.pdf")
# ggplot(data=Dmel_CG42335, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# gstD9 <- plotCounts(dds, gene="Dmel_CG10091", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG10091_gstD9_counts.pdf")
# ggplot(data=gstD9, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG10091_gstD9_log2.counts.pdf")
# ggplot(data=gstD9, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# ##############################################
# ########## G+I
# ############################################
# ### sordd2 is significant for both G and GxI ###
# sordd2 <- plotCounts(dds, gene="Dmel_CG32581", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG32581_sordd2_counts.pdf")
# ggplot(data=sordd2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG32581_sordd2_log2.counts.pdf")
# ggplot(data=sordd2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# Dmel_CG32640 <- plotCounts(dds, gene="Dmel_CG32640", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG32640_hypProt_counts.pdf")
# ggplot(data=Dmel_CG32640, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG32640_hypProt_log2.counts.pdf")
# ggplot(data=Dmel_CG32640, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# mthl8 <- plotCounts(dds, gene="Dmel_CG32475", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG32475_mthl8_counts.pdf")
# ggplot(data=mthl8, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG32475_mthl8_log2.counts.pdf")
# ggplot(data=mthl8, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# dysf <- plotCounts(dds, gene="Dmel_CG32474", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG32474_dysf_counts.pdf")
# ggplot(data=dysf, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG32474_dysf_log2.counts.pdf")
# ggplot(data=dysf, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# antisenseRNACR45330 <- plotCounts(dds, gene="Dmel_CR45330", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CR45330_antisenseRNACR45330_counts.pdf")
# ggplot(data=antisenseRNACR45330, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CR45330_antisenseRNACR45330_log2.counts.pdf")
# ggplot(data=antisenseRNACR45330, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# tRNA_m6t6A37_MT <- plotCounts(dds, gene="Dmel_CG18853", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG18853_tRNA_m6t6A37_MT_counts.pdf")
# ggplot(data=tRNA_m6t6A37_MT, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG18853_tRNA_m6t6A37_MT_log2.counts.pdf")
# ggplot(data=tRNA_m6t6A37_MT, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CG31997 <- plotCounts(dds, gene="Dmel_CG31997", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG31997_hypProt_counts.pdf")
# ggplot(data=Dmel_CG31997, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG31997_hypProt_log2.counts.pdf")
# ggplot(data=Dmel_CG31997, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# H_ACA_box_snoRNA_gene <- plotCounts(dds, gene="Dmel_CR34601", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CR34601_H_ACA_box_snoRNA_gene_counts.pdf")
# ggplot(data=H_ACA_box_snoRNA_gene, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CR34601_H_ACA_box_snoRNA_gene_log2.counts.pdf")
# ggplot(data=H_ACA_box_snoRNA_gene, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Tep2 <- plotCounts(dds, gene="Dmel_CG7052", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG7052_Tep2_counts.pdf")
# ggplot(data=Tep2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG7052_Tep2_log2.counts.pdf")
# ggplot(data=Tep2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# RpL22like <- plotCounts(dds, gene="Dmel_CG9871", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG9871_RpL22-like_counts.pdf")
# ggplot(data=RpL22like, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG9871_RpL22-like_log2.counts.pdf")
# ggplot(data=RpL22like, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# his1 <- plotCounts(dds, gene="Dmel_CG33801", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG33801_his1_counts.pdf")
# ggplot(data=his1, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG33801_his1_log2.counts.pdf")
# ggplot(data=his1, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# arc2 <- plotCounts(dds, gene="Dmel_CG13941", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG13941_arc2_counts.pdf")
# ggplot(data=arc2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG13941_Arc2_log2.counts.pdf")
# ggplot(data=arc2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# GALNT1 <- plotCounts(dds, gene="Dmel_CG31776", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG31776_GALNT1_counts.pdf")
# ggplot(data=GALNT1, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG31776_GALNT1_log2.counts.pdf")
# ggplot(data=GALNT1, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# HIG1 <- plotCounts(dds, gene="Dmel_CG11825", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG11825_HIG1_counts.pdf")
# ggplot(data=HIG1, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG11825_HIG1_log2.counts.pdf")
# ggplot(data=HIG1, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Cyp6d2 <- plotCounts(dds, gene="Dmel_CG4373", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG4373_Cyp6d2_counts.pdf")
# ggplot(data=Cyp6d2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG4373_Cyp6d2_log2.counts.pdf")
# ggplot(data=Cyp6d2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# hml <- plotCounts(dds, gene="Dmel_CG7002", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG7002_hml_counts.pdf")
# ggplot(data=hml, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG7002_hml_log2.counts.pdf")
# ggplot(data=hml, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# damm <- plotCounts(dds, gene="Dmel_CG18188", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG18188_damm_counts.pdf")
# ggplot(data=damm, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG18188_damm_log2.counts.pdf")
# ggplot(data=damm, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# hist1 <- plotCounts(dds, gene="Dmel_CG33816", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG33816_hist1_counts.pdf")
# ggplot(data=hist1, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG33816_hist1_log2.counts.pdf")
# ggplot(data=hist1, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CR45631 <- plotCounts(dds, gene="Dmel_CR45631", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CR45631_lncRNA_counts.pdf")
# ggplot(data=Dmel_CR45631, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CR45631_lncRNA_log2.counts.pdf")
# ggplot(data=Dmel_CR45631, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CG10013 <- plotCounts(dds, gene="Dmel_CG10013", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG10013_hypProt_counts.pdf")
# ggplot(data=Dmel_CG10013, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG10013_hypProt_log2.counts.pdf")
# ggplot(data=Dmel_CG10013, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# Dmel_CR46123 <- plotCounts(dds, gene="Dmel_CR46123", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CR46123_lncRNA_counts.pdf")
# ggplot(data=Dmel_CR46123, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CR46123_lncRNA_log2.counts.pdf")
# ggplot(data=Dmel_CR46123, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# miple2 <- plotCounts(dds, gene="Dmel_CG18321", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG18321_miple2_counts.pdf")
# ggplot(data=miple2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG18321_miple2_log2.counts.pdf")
# ggplot(data=miple2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# ##############

# sordd2 <- plotCounts(dds, gene="Dmel_CG32581", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG32581_sordd2_counts.pdf")
# ggplot(data=sordd2, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG32581_sordd2_log2.counts.pdf")
# ggplot(data=sordd2, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()

# #####



# bgcn <- plotCounts(dds, gene="Dmel_CG30170", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG30170_bgcn_counts.pdf")
# ggplot(data=bgcn, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG30170_bgcn_log2.counts.pdf")
# ggplot(data=bgcn, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# ####################################################
# ################ NS plots ###########################
# #####################################################
# YL1 <- plotCounts(dds, gene="Dmel_CG4621", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG4621_YL-1_counts.pdf")
# ggplot(data=YL1, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG4621_YL-1_log2.counts.pdf")
# ggplot(data=YL1, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# hts <- plotCounts(dds, gene="Dmel_CG43443", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG43443_hts_counts.pdf")
# ggplot(data=hts, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG43443_hts_log2.counts.pdf")
# ggplot(data=hts, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# bam <- plotCounts(dds, gene="Dmel_CG10422", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG10422_bam_counts.pdf")
# ggplot(data=bam, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG10422_bam_log2.counts.pdf")
# ggplot(data=bam, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()


# sxl <- plotCounts(dds, gene="Dmel_CG43770", intgroup="group", returnData=TRUE)
# pdf("KallistoDESeq2_Dmel_CG43770_sxl_counts.pdf")
# ggplot(data=sxl, aes(x=group, y=count)) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
# pdf("KallistoDESeq2_Dmel_CG43770_sxl_log2.counts.pdf")
# ggplot(data=sxl, aes(x=group, y=log2(count))) + geom_boxplot() + scale_x_discrete(guide = guide_axis(angle = 90)) + geom_point() + theme_minimal()
# dev.off()
