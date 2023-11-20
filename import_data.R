# Nov 2023
# HBC single-cell RNA-seq workshop


# Single-cell RNA-seq analysis - QC set-up
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(readr)
library(stringr)
library(ggplot2)

#---load data--#

##--- Read in data using readMM() from Matrix package to create a sparse matrix to reduce memory(RAM), processing capacity(CPU) and storage
# Read in `matrix.mtx`
counts <- readMM("data/ctrl_raw_feature_bc_matrix/matrix.mtx.gz")

# Read in `genes.tsv`
genes <- read_tsv("data/ctrl_raw_feature_bc_matrix/features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1

# Read in `barcodes.tsv`
cell_ids <- read_tsv("data/ctrl_raw_feature_bc_matrix/barcodes.tsv.gz", col_names = FALSE)$X1

# Make the column names as the cell IDs and the row names as the gene IDs
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids


##--- read in data using Read10X() function from Seurat package and create a seurat object ----

### ---- Reading in a Single Sample -----

# How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")
# Turn count matric into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100) #This argument will filter out poor quality cells that likely just have random barcodes encapsulated without any cell present.
                                                #Usually, cells with less than 100 genes detected are not considered for analysis.

# Seurat **automatically creates some metadata for each of the cells** when you use the Read10X() function to read in data.
# This information is stored in the meta.data slot within the Seurat object.
# Explore the metadata
head(ctrl@meta.data)
# orig.ident -> sample identity if known, default -> SeuratProject
# nCount_RNA -> # of UMIs per cell
# nFeature_RNA -> # of genes detected per cell


### ---- Reading in multiple samples with a `for loop` -----

# In R, the for loop has the following structure
# for (variable in input){
#       command1
#       command2
#       command3
# }

# Use for loop to interate over the two sample folders and execute two commands for each sample
# as steps took in reading a single sample
# (1)read in the cound data
# (2) create the Seurat object

for (file in c("ctrl_raw_feature_bc_matrix","stim_raw_feature_bc_matrix")){
  seurat_data <- Read10X(data.dir = paste0("data/",file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj) # this command assigns the created seurat_obj to a new variable
}

# check the metadata
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)