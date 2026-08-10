rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
levels(Idents(HSC))

HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.5)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
DimPlot(HSC, label = T,pt.size = 1)

#HSC <- NormalizeData(HSC)

####1.读取目标genneset文件####
genesets <- read.csv("function_mouse.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

####2.提取分组平均表达矩阵####
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(HSC, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

####3.GSVA富集分析####
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)

####gsva.df####
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row")
pheatmap::pheatmap(gsva.res, show_colnames = T)

exprTable_t <- as.data.frame(t(gsva.res))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
pheatmap::pheatmap(gsva.res, cluster_cols = hclust_1)


manual_order = c("HSC signature","MK","Ery","Gr","Ly","Quiescence","S_phase", "G2M_phase")

dend = reorder(as.dendrogram(hclust_1), wts=order(match(colnames(exprTable_t),manual_order)))

# 默认为mean，无效时使用其他函数尝试
# dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)

row_cluster <- as.hclust(dend)
pheatmap::pheatmap(gsva.res, cluster_rows = row_cluster)
