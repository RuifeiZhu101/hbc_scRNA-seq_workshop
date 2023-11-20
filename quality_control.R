# Nov 2023
# HBC single-cell RNA-seq workshop

# Single-cell RNA-seq analysis - QC
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



##---- Merge objects together into a single Seurat Object -----

# This step will make it easier to run the QC steps for both sample groups
# and enable us to easily compare the data quality for all the samples
# we can use merge() function from the Seurat package to do this

merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix,
                       y = stim_raw_feature_bc_matrix, # y could be multiple, e.g. c("sample2","sample3"...)
                       add.cell.ids = c("ctrl", "stim")) # this add.cell.ids add sample-specific prefix to each of our cell IDs
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)


## Generating quality metrics

# Explore merged metadata
View(merged_seurat@meta.data)

# In order to create appropriate plots for the QC analysis
# we need to calculate some additional metrics
# 1. num of genes detected per UMI(more genes detected per UMI, more complex our data)
# 2. mitochondrial ratio: give us percentage of cell reads comes from mitochondrial genes

### Novelty score = log10(#genes)/log10(#UMI) per cell
# Add number of genes per UMI for each cell to metadata
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

### Mitochondrial Ratio
# Seurat has a convenient function `PercentageFeatureSet()`
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat,
                                                pattern = "^MT-") #  The pattern provided ("^MT-") works for human gene names.
                                                                  # You may need to adjust the pattern argument depending on your organism of interest. 
#it is advisable to manually compute this metric. see detail (https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md).

# update using ratio rather than percentage in this analysis
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100 
head(merged_seurat@meta.data)

### Additional metadata columns
# Aim: add additional information like cell IDs(currently the rownames of metadata) and conditions as columns to the metadata
# To do so, we would like to extract the metadata as a separate data frame to avoid taking risks to affect the seurat objects

# 1. Create metadata dataframe
metadata <- merged_seurat@meta.data

# 2. Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# 3. Add sample information to metadata
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells,"^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells,"^stim_"))] <- "stim"
View(metadata)

# 4. rename some of the existing columns to be more intuitive
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
View(metadata)


## Saving the updated metadata to our Seurat object

#Add metadata back to Seurat Object
merged_seurat@meta.data <- metadata
#Create .RData object to load at any time
save(merged_seurat, file = "data/merged_filtered_seurat.RData")


#### Assessing quality metrics

###cell counts
# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# We see over 15,000 cells per sample, which is quite a bit more than the 12-13,000 expected

### UMI counts
# The UMI counts per cell should generally be above 500(Why? due to cell ranger second calling), that is the low end of what we expect. 
# cells with less than 500 UMI would be called as a real cell
# However, this threshold might need to be adjusted due to the real condition
# For example, samples are highly heterogeneous that contains single cells with very high/low RNA contents
# such as neutrophils

#Visualize the number UMIs/transcripts per cell
metadata %>%
  ggplot(aes(color = sample, x = nUMI, fill = sample))+
  geom_density(alpha = 0.2) +
  scale_x_log10()+
  theme_classic()+
  ylab("cell density")+
  geom_vline(xintercept = 500) 

###Genes detected per cell
# similar expectations for gene detection as for UMI detection, slightly lower
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


### Complexity
# generally we expect the novelty score to be above 0.8(nGene/nUMIs ~ 6.7) for good quality cells
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


### Mitochondrial Counts Ratio
# to identify mitochondrial contamination from dead/dying cells
# usually we define poor quality samples with mitochondrial ratio surpass 0.2
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

### Reads per cell is another metric that can be useful to explore
# generally, with this metric you hope to see all of the samples with peaks in relatively the same location 
# between 10,000 and 100,000 reads per cell.


### Joint filtering effects
# A general rule of thumb when performing QC is to set thresholds for individual metrics to be as permissive as possible
# and always consider the joint effects to reduce the risk of filtering out  any viable cells

#  Visualize the correlation between genes detected and number of UMIs 
# and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs (upper right quadrant of the plot).
# Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. 
# With this plot we also evaluate the slope of the line, and any scatter of data points in the bottom right hand quadrant of the plot. 
# These cells have a high number of UMIs but only a few number of genes. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# Mitochondrial read fractions are only high in particularly low count cells with few detected genes (darker colored data points). 
# This could be indicative of damaged/dying cells whose cytoplasmic mRNA has leaked out through a broken membrane, and thus, only mRNA located in the mitochondria is still conserved. 
# We can see from the plot, that these cells are filtered out by our count and gene number thresholds.

## Filtering
### Cell-level filtering
# Based on the above metrics, we can decide the thresholds to remove low quality cells
# often the recommendations mentioned earlier are a rough guideline, and the specific experienment needs to inform the exact thresholds chosen
# For this analysis we will use:
# nUMI > 500
# nGene > 250
# log10GenesPerUMI >0.8
# mitoRatio < 0.2

dim(merged_seurat)

# Filter out low quality cells using selected thresholds
filtered_seurat <- subset(x = merged_seurat,
                          subset = (nUMI >= 500)&
                            (nGene >= 250)&
                            (log10GenesPerUMI > 0.80)&
                            (mitoRatio < 0.20))
filtered_seurat
dim(filtered_seurat)
head(rownames(filtered_seurat))
head(colnames(filtered_seurat))
names(filtered_seurat)

### Gene-level filtering
# genes with zero counts can dramatically reduce the average expression for a cell
# start by identifying which genes have a zero count in each cell

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
head(colnames(counts))
# create a logical matrix specifying for each gene on whether or not zero count per cell
nonezero <- counts >0

# keep only genes which are expressed in 10 or more cells
keep_genes <- Matrix::rowSums(nonezero) >= 0
filtered_counts <- counts[keep_genes,]

# reassign to filtered seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Re-assess QC metrics
metadata_filtered <- filtered_seurat@meta.data
###cell counts
metadata_filtered %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# There are just under 15K cells left for both the control and stim cells. 
# The number of cells removed is reasonably low.
# While it would be ideal to have 12K cells, we do not expect that due to the lower capture efficiency (i.e. the number of actual cells encapsulated within droplets containing barcodes) of these technologies. 
# If we still see higher than expected numbers of cells after filtering, this means we could afford to filter more stringently (but we don't necessarily have to).


### UMI counts
metadata_filtered %>%
  ggplot(aes(color = sample, x = nUMI, fill = sample))+
  geom_density(alpha = 0.2) +
  scale_x_log10()+
  theme_classic()+
  ylab("cell density")+
  geom_vline(xintercept = 500) 

###Genes detected per cell
metadata_filtered %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)
# The smaoll shoulder may represent the cell populations with higher complexity?
# exhibit more diversity in its transcriptome (with the larger number of genes detected).


### Complexity
metadata_filtered %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)


### Mitochondrial Counts Ratio
metadata_filtered %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

### Joint filtering effects
metadata_filtered %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)
