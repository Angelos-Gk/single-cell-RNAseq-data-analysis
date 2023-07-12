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

# run SingleR with multiple reference datasets (default mode)
hpca <- celldex::HumanPrimaryCellAtlasData()
dice <- celldex::DatabaseImmuneCellExpressionData()

# ========== Strategy 1: Using reference-specific labels ==========
hpca$label.main
dice$label.main

# adding ref info to labels 
# because some labels are shared in both datasets
# we label them so we know form which dataset they originate
hpca$label.main <- paste0('HPCA.',hpca$label.main)
dice$label.main <- paste0('DICE.',dice$label.main)

# create a combined ref based on shared genes
shared <- intersect(rownames(hpca),rownames(dice))
combined <- cbind(hpca[shared,],dice[shared,])
combined$label.main

# run SingleR using combined ref
# saving counts into a separate object
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels)

pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data),rownames(com.res1)), 'labels']
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)

# ========== Strategy 2: Comparing scores across references ==========
# remove the labels "HPCA" and "DICE"
hpca$label.main <- gsub('HPCA\\.','',hpca$label.main )
dice$label.main <- gsub('DICE\\.','',dice$label.main )
hpca$label.main 
dice$label.main

com.res2 <- SingleR(test = pbmc_counts, 
        ref = list(HPCA = hpca, DICE = dice),
        labels = list(hpca$label.main, dice$label.main))

# Check the final label from the combined assignment
table(com.res2$labels)

# Which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.',(com.res2$reference))
best_ref <- as.data.frame(split(com.res2,grouping))

#get differentially expressed genes from each individual references 
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes

#combined diagnostics
plotScoreHeatmap(com.res2)

# ========== Strategy 3: Using harmonized labels ==========
hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same set of genes
shared <- intersect(rownames(hpca.ont),rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)),10)

# Using label.ont instead of label.main while running SingleR
com.res3 <- SingleR(test = pbmc_counts,
        ref = list(HPCA = hpca.ont, DICE = dice.ont),
        labels = list(hpca.ont$label.ont, dice.ont$label.ont))

table(com.res3$labels)

# How to map cell ontology terms
colData(hpca.ont)
colData(dice.ont)

hpca.file <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.file, header = F)

dice.file <- system.file("mapping","dice.tsv", package = "celldex")
dice.mapping <- read.delim(dice.file, header = F)

sessionInfo()

