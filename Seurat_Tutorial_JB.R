library(dplyr)
library(Seurat)
library(patchwork)

# 1. Open Data and create Seurat Object
# Download data, decompress with Zip7 and store in your desired folder

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/Admin/Desktop/R-script/Seurat Tutorial JB/Raw_Data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 709591472
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# 29905192 bytes
sparse.size <- object.size(pbmc.data)
sparse.size

# 23.7 bytes
dense.size/sparse.size

# 2. Pre-Processing the Data

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# After this preliminary analysis, we use these data to filter cells
# We filter cells that have unique feature (=gene) counts over 2,500 (or less than 200) (more seems like duplicats, less like empty droplets)
# We filter cells that have >5% mitochondrial counts (dying cells have abberant high mitochondrial gene expression)

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# visualize subset now
pbmc

# 3. Normalization of Data (feature expression measurements for each cell by the total expression)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# To check out which commands you have already performed on your Seurat Object
str(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 4. Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# See the names of the 10 most variable genes
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 5. Scaling the Data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 6. Perform linear dimensional reduction (creates a list with genes with the most positive and negative loadings)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the dimensionality of the dataset
ElbowPlot(pbmc)

# 7. Clustering the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95965
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8723
## Number of communities: 9
## Elapsed time: 0 seconds

# Find clusters can also be performed in a different resolution (0.1 - 1), the lower the number the fewer clusters
# To look at the metadata and see the different resolutions in columns to choose which one works best for you
view(pbmc@meta.data)

# Visualize data reduction plot at different resolutions to identify which one fits best
# Avoid overlapping clusters in higher resolutions
DimPlot(pbmc, group.by = "RNA_snn_res.0.5", label = TRUE)

# To set identity of clusters (just do this, if you want to join resolution)
Idents(pbmc) <- "RNA_snn_res.0.5"

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# 8. Non-dimensional reduction (tSNE/UMAP)
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

# Save this object from here, that it can easily be opened from here
saveRDS(pbmc, file = "C:/Users/Admin/Desktop/R-script/Seurat Tutorial JB/Output/pbmc_tutorial.rds")

# 9. Finding differentially expressed features (cluster biomarkers)

# Activate presto package for faster DE (differential expression) analysis
library(presto)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Plot expression probability distributions of individual gene markers across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features.
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# 10. Assigning cell type identity to clusters
# Known markers for specific cell types
#   Cluster ID 	Markers 	    Cell Type
#   0 	        IL7R, CCR7 	  Naive CD4+ T
#   1 	        CD14, LYZ 	  CD14+ Mono
#   2 	        IL7R, S100A4 	Memory CD4+
#   3 	        MS4A1 	      B
#   4 	        CD8A 	CD8+    T
#   5 	        FCGR3A, MS4A7 FCGR3A+ Mono
#   6 	        GNLY, NKG7 	  NK
#   7 	        FCER1A, CST3 	DC
#   8 	        PPBP 	        Platelet

#renaming clusters just individually because you already know the names
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#Plotting and saving plot as image
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

#Saving final document
saveRDS(pbmc, file = "C:/Users/Admin/Desktop/R-script/Seurat Tutorial JB/Output/pbmc3k_final.rds")