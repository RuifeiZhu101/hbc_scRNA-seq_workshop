# Dec 2023
# hbc_scRNA-seq training - culstering analysis


# Single-cell RNA-seq clustering analysis


## Learning Objectives:

# -   Describe methods for evaluating the number of PC used for clustering
# -   Perform clustering of cells baksed on significant PCs

#After performing integratation, we have high quality cells from different samples/batches integrated, we want to know the different cell types presented whithin our population of cells.

# Goals:
# 
# -   To generate cell type-specific clusters*
# -   use known cell type marker genes to determine the identities of the clusters.
# -   To determine whether clusters represent true cell types or cluster due to biological or technical variation**, such as clusters of cells in the S phase of the cell cycle, clusters of specific batches, or cells with high mitochondrial content.*
# 
# Challenges:
# 
# -   Identifying poor quality clusters that may be due to uninteresting biological or technical variation
# -   Identifying the cell types of each cluster
# 
# Recommendations:
# 
# -   Have a good idea of your expectations for the cell types to be present prior to performing the clustering. Know whether you expect cell types of low complexity or higher mitochondrial content AND whether the cells are differentiating*
# -   If you have more than one condition, it's often helpful to perform integration to align the cells
# -   Regress out number of UMIs (by default with sctransform), mitochondrial content, and cell cycle, if needed and appropriate for experiment, so not to drive clustering
# -   Identify any junk clusters for removal or re-visit QC filtering. Possible junk clusters could include those with high mitochondrial content and low UMIs/genes. If comprised of a lot of cells, then may be helpful to go back to QC to filter out, then re-integrate/cluster.
# -   If not detecting all cell types as separate clusters, try changing the resolution or the number of PCs used for clustering



# Clustering cells based on top PCs(metagenes)

## Set up
# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

## Identify significant PCs
# Seurat assigns cells to clusters based on their PCA scores derived from the expression of the integrated most variable genes,
# Determinint how many PCs to include in the clustering step is important to ensure that 
# we can capture the majority of the variation/cell types present in our database

# ways to explore PCs before clustering:
##(1)Heatmap
## The idea here is to look at the PCs and determine whether the genes driving them make sense for differentiating the different cell types

# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, # specifies the number of cells with the most negative or postive PCA scores to use for the plotting.
           balanced = TRUE)
# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

##(2) Elbow plot
#The elbow plot visualizes the standard deviation of each PC, and we are looking for where the standard deviations begins to plateau. 
#Essentially, where the elbow appears is usually the threshold for identifying the majority of the variation. 

# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
# Based on this plot, we could roughly determine the majority of the variation by where the elbow occurs around PC8 - PC10
# However, this method can be quite subjective.
# a quantitative method can be found in the script elbow_plot_metric.R

## Note: While the above 2 methods were used a lot more with older methods from Seurat for normalization and identification of variable genes, 
## they are no longer as important as they used to be. 

## The older methods incorporated some technical sources of variation into some of the higher PCs, so selection of PCs was more important. 
## SCTransform estimates the variance better and does not frequently include these sources of technical variation in the higher PCs.

## In theory, with SCTransform, the more PCs we choose the more variation is accounted for when performing the clustering, but it takes a lot longer to perform the clustering. 
## Therefore for this analysis, we will use the first 40 PCs to generate the clusters.


## Cluster the cells
# Seurat uses a graph-based clustering approach using a K-nearest neibor
# A nice in-depth description of clustering methods(https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/clustering-and-cell-annotation.html).

### Find neighbors
# The first step is to construct a K-nearest neighbor (KNN) graph based on the euclidean distance in PCA space.
# Using Seurat:
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)
### Find clusters
#Next, Seurat will iteratively group cells together with the goal of optimizing the standard modularity function.

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
#The resolution is an important argument that sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment.
#For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. 
# Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.


## Visualize clusters of cells
# The most popular dimensionality reduction methods: 
# t-distributed stochastic neighbor embedding (t-SNE) and Uniform Manifold Approximation and Projection (UMAP)
# Here will use UMAP

#We can only visualize the results of one resolution setting at a time. 
#Look at the metadata of our Seurat object, a separate column for each of the different resolutions calculated can be found.

# Explore resolutions
seurat_integrated@meta.data %>% 
  View()

# To choose a resolution to start with, we often pick something in the middle of the range like 0.6 or 0.8.
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

## Calculation of UMAP

seurat_integrated <- RunUMAP(seurat_integrated,
                 reduction = "pca",
                 dims = 1:40)

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# because my clusters looks different from the lesson
# Load in the object to your R session and overwrite the existing one:
load(bzfile("data/additional_data/seurat_integrated.RData.bz2"))

# repeat the steps
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
# now they look identical!
