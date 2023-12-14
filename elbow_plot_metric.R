# # Dec 2023
# hbc_scRNA-seq training -  Elbow plot: quantitative approach

# First, try top 40 dimensions:

# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)

# We can determin the elbow using a more quantitative approach
# We can calculate where the PCs start to elbow by assessing two metrics:
# 1. The point where the PCs only contribute 5% of SD and the PCs cumulatively contribute 90% of the SD
# 2. The point where the percent change in variation between the consecutive PCs < 0.1%

## Start by calculating the first metric:
# Determine percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
# The first metric returns PC42 as the PC matching these requirements. 

## Check the second metric:
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

#This second metric returns PC18.

# Usually,choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.
# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

