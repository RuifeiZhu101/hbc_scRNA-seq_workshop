# Nov 2023
# HBC single-cell RNA-seq workshop
# Single-cell RNA-seq analysis - Normalization and regressiong out unwanted variation

## Goals:
# To normalize the gene expression values to account for differences in sequencing depth and overdispersed count values.
# To identify the most variant genes likely to be indicative of the different cell types present.

## Challenges:
# Checking and removing unwanted variation to avoid artifact cell clusters in downstream analysis

## Recommendations:
# Have an expectation for the cell types to be present prior to performing the clustering.
# know whether you expect cell types of low complexity? higher mitochondrial content? cells are differentiating?

# Regress out num of UMIs(default using sctransform)


# Load libraries
library(cowplot)
library(Seurat)
library(tidyverse)
library(RCurl)
library(ggplot2)

# Normalization: make expression counts comparable across genes and/or samples
# Two main factors: sequencing depth, gene length(for full-length sequencing only)
# Two main steps: scaling(get all cells to have the same UMI counts), transformation(log(simple) trans/ Pearson residual trans)

### Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

### Evaluating effects of cell cycle 

# The most common biological data correction in scRNA-seq is the effects of the cell cycle on the transcriptome. 
# Assign each cell a score based on its expression of G2/M and S phase markers
# If not working with human data, additional materials detailing how to acquire cell cycle markers for other organisms of interest
# can be found in the workshop(cell_cycle_scoring.md).

# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)                                

# Using PCA to evaluate the effects of cell cycle

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by(split by) cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# no large differences due to cell cycle phase were found from this plot
# -> no need to regress out variations due to the cell cycle 
# regress out cell cycle strategies (https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html)


### Evaluating effects of mitochondrial expression

# oftentimes, it is useful to regress out variation due to mitochondrial expression
# if mito differences refresent a biological pheno that helps to distinguish cell cluster, do not regress out

# step 1. turn the mito ratio into a new catagorical variable based on quartiles

# check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into catagorical factor vector 
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks = c(-Inf, 0.01438, 0.01993,0.02669,Inf),
                                     labels = c("Low","Medium","Medium High", "High"))

# step 2. plot PCA
# Plot the PCA colored by cell mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")
# there is a pattern related to the mitoRatio in cell cluster 
# we need to regress out mitochondrial fraction as a source of unwanted variation


### Normalization and regressing out sources of unwanted variation using SCTransform (https://satijalab.org/seurat/articles/sctransform_vignette.html)

# SCTransform function: 
# 1. normalize data (constructs a generalized linear model (GLM) for each gene
# with UMI counts as the response and sequencing depth as the explanatory variable. )
# 2. performs a variance stabilization and regresses out additional covariates

# There are two samples in our dataset (from two conditions), 
# need to keep them as separate objects and transform them as that is what is required for integration.

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# SCT can generate large R objects/variables in terms of memory
# Adjust the limit for allowable object sizes within R(Default is 500 * 1024^2 = 500MB)
options(future.globals.maxSize = 4000 * 1024^2)
# perform the sct on all samples
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
}


# Check which assays are stored in objects
split_seurat$ctrl@assays
split_seurat$stim@assays
# A SCT component was added to the assays slot. 
# the most variable features will be only genes stored inside the SCT assay.
# move through the scRNA-seq analysis, the most appropriate assay needs to be chosen for different steps in the analysis

### Save the object!
# t can take a while to get back to this stage especially when working with large datasets, 
# it is best practice to save the object as an easily loadable file locally.

# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")


