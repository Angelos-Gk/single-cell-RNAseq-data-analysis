#setwd("C:/Users/angel/Desktop/Rthings/Bioinformagician/scRNAseq analysis")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Input Data (10X CellRanger, .HDF5 format) ------------------
hdf5_obj <- Read10X_h5(filename="20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)

pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)

# QC and Filtering ------------------
# calculate % of mitochondrial genes
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 & # remove cells with less than 800 counts
         nFeature_RNA > 500 & # remove cells with less than 500 genes
         mitoPercent < 10) # remove cells with mitochondrial percentage greater than 10%

# check the difference
pbmc.seurat
pbmc.seurat.filtered

# Pre-process standard workflow ------------------ 
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# check the results
View(pbmc.seurat.filtered@meta.data)
DimPlot(pbmc.seurat.filtered, reduction = 'umap')

# Get reference data ------------------  
ref <- celldex::HumanPrimaryCellAtlasData()
ref
# check which celltypes are included in the reference 
# in order to see if it includes our data's celltypes
View(as.data.frame(colData(ref)))

# Run SingleR (default mode) ------------------  
# default means that SingleR performs annotation of each individual cell
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot ='counts')

pred <- SingleR(test = pbmc_counts,
        ref = ref,
        labels = ref$label.main)

pred

pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'singleR.labels')

# Annotation diagnostics ------------------  
# check how well SingleR performed the annotation

# Annotation diagnostics based on the scores within cells
pred
pred$scores

plotScoreHeatmap(pred)

# Annotation diagnostics based on deltas across cells
plotDeltaDistribution(pred)

# Annotation diagnostics via comparing to supervised clustering
tab <- table(Assigned=pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))

sessionInfo()
