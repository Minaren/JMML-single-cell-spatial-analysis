library(Seurat) 
library(SeuratDisk)
load("D:/aJMML/scRNA1/sce_har_anno.RData")
SaveH5Seurat(sce.harm, filename = "seurat_obj.h5Seurat")
Convert("seurat_obj.h5Seurat", dest = "h5ad")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("SeuratDisk")
