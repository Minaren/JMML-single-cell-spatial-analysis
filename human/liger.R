####数据预处理####
##==准备seurat对象列表==##
library(Seurat)
library(tidyverse)
rm(list=ls())
dir = c('D:/aJMML/JMML/JMMLID5',
        'D:/aJMML/JMML/PBM1')
samples_name= c('JMML','Normal')
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features = 200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 10) 
}   
saveRDS(scRNAlist, "scRNAlist.rds")

####LIGER整合数据####
##==多个scRNA-seq数据整合==##
library(Seurat)
#install.packages("liger")
library(liger)
library(tidyverse)
library(patchwork)
#BiocManager::install("SeuratWrappers")
#devtools::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
#为了方便展示效果，只取其中的2个样本演示
scRNA <- merge(scRNAlist[[1]], scRNAlist[[2]]) -> scRNA.orig
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, split.by="orig.ident", do.center=FALSE)
nFactors=20    #设置矩阵分解的因子数，一般取值20-40
##因式分解
scRNA <- RunOptimizeALS(scRNA, k=nFactors, split.by="orig.ident")
##多样本整合
scRNA <- RunQuantileNorm(scRNA, split.by="orig.ident")
#整理因子顺序
scRNA$clusters <- factor(scRNA$clusters, 
                         levels=1:length(levels(scRNA$clusters)))
##聚类
scRNA <- FindNeighbors(scRNA, reduction="iNMF", 
                       dims=1:nFactors) %>% FindClusters()
scRNA <- RunUMAP(scRNA, dims=1:nFactors, reduction="iNMF")
##不整合数据的降维聚类
scRNA.orig <- scRNA.orig %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
  RunPCA() %>% FindNeighbors(dims=1:20) %>% FindClusters() %>% RunUMAP(dims=1:20)
##可视化
p1 = DimPlot(scRNA, group.by="orig.ident", pt.size=0.05) + ggtitle("Integrated by liger")
p2 = DimPlot(scRNA.orig, group.by="orig.ident", pt.size=0.05) + ggtitle("No integrated")
p3 = DimPlot(scRNA, group.by="seurat_clusters", label=T, label.size=2) + ggtitle("Clustered by seurat")#seurat_clusters列是seurat聚类的结果，clusters列是liger聚类的结果，其聚类数量与RunOptimizeALS函数运行时k参数的值相同。
p4 = DimPlot(scRNA, group.by="clusters", label=T, label.size=2) + ggtitle("Clustered by liger")

plot1 = p1 + p2 + plot_layout(guides = 'collect')
plot2 = p3|p4

ggsave("integrated_liger.png", plot=plot1, width=8, height=3.6)
ggsave("clustered_liger.png", plot=plot2, width=10, height=6)
save(scRNA,file = "sce_liger.RData")
