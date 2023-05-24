
# load libraries
library(Seurat)
library(tidyverse)

# load the dataset (Non Small Cell Lung Cancer)
nsclc.sparse.m <- Read10X_h5(filename='C:/Users/angel/Desktop/Rthings/Bioinformagician/scRNAseq analysis/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
# check modalities
str(nsclc.sparse.m)
# get counts for gene expression
cts <- nsclc.sparse.m$`Gene Expression`
head(cts)

# initialize the Seurat object with the raw (non normalized) data
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj


# 1) QC ---------------------------------
# view the metadata of the seurat object
View(nsclc.seurat.obj@meta.data)
# % MT reads
# add another column to the metadata with the percent of mitochondrial genes in each cell
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj,feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


# 2) Filtering ---------------------------------
# filter low quality cells based on the number of reads and mt genes percentage
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                             percent.mt < 5)


# 3) Normalization ---------------------------------
# divide the gene expression measurements in each cell by the total expression, 
# multiply it by a scaling factor and the log-transform it
# we do all that in order to get the expression measurements in relative measure so we can compare 
# across different cells
# nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "logNormalize", scale.factor = 10000)
# the above are the default values so we can simple:
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# we can see the commands and their parameters 
# that we have run on our seurat object by using str()
# and checking the command slot
str(nsclc.seurat.obj@commands)


# 4) Identify highly variable features ---------------------------------
# we want to see only those features that exhibit high cell-to-cell variation
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj),10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5) Scaling ---------------------------------
# we need to account for unwanted sources of variation
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)


# 6) Perform linear dimensionality reduction ---------------------------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)
# we need to capture all the principal components that explain the higher percentage of the variance
# in this case we can use all the PCs up to PC 15


# 7) Clustering ---------------------------------
# we want to cluster together cells that have similar expression patterns
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
# resolution parameter: the lower the number the lower the clusters
# try different number and see what works best
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3,0.5,0.7,1))
View(nsclc.seurat.obj@meta.data)

# visualize how many clusters you get for each resolution and choose the best for you
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE) # good clustering
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE) # bad clustering (overlaps)

# setting identity of clusters
Idents(nsclc.seurat.obj) # it by default picked 19-cluster resolution
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1" # we set it to our preference
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction ---------------------------------
# install UMAP by reticulate::py_install(packages = 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# you can set 'label=TRUE' or use LabelClusters to help label individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")

sessionInfo()
