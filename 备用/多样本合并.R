setwd("JMML")
library(Seurat)
library(tidyverse)
library(patchwork)
dir.create('cluster1')
dir.create('cluster2')
dir.create('cluster3')
set.seed(123) #设置随机数种子，使结果可重复
####合并数据集####
##使用目录向量合并
dir = c('D:/aJMML/JMML/Kras', 
        'D:/aJMML/JMML/WT')

names(dir) = c('Kras', 'WT')
dir
counts <- Read10X(data.dir = dir)
scRNA1 = CreateSeuratObject(counts, min.cells=1)
dim(scRNA1)#查看基因数和细胞总数
table(scRNA1@meta.data$orig.ident)  #查看每个样本的细胞数

scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))
plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA1, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)
#选取主成分
pc.num=1:30

##细胞聚类
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)
metadata <- scRNA1@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

##非线性降维
#tSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
embed_tsne <- Embeddings(scRNA1, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNA1, reduction = "tsne", label=T) 
ggsave("cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
embed_umap <- Embeddings(scRNA1, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("cluster1/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("cluster1/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)
setwd("D:/aJMML/JMML")
#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
ggsave()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, min.cells=1)
}
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##利用锚点整合数据，运行时间较长
scRNA3 <- IntegrateData(anchorset = scRNA.anchors)
saveRDS(scRNA3, file = "scRNA3.Rds")
saveRDS(scRNA1, file = "scRNA1.Rds")



scRNA3 <- NormalizeData(scRNA3)
scRNA3 <- FindVariableFeatures(scRNA3, selection.method = "vst")
scRNA3 <- ScaleData(scRNA3, features = VariableFeatures(scRNA3))

scRNA3<- readRDS('scRNA3.rds')
scRNA3 <- RunPCA(scRNA3, features = VariableFeatures(scRNA3))
plot1 <- DimPlot(scRNA3, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA3, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster2/pca.png", plot = plotc, width = 8, height = 4)
#选取主成分
pc.num=1:30

##细胞聚类
scRNA3 <- FindNeighbors(scRNA3, dims = pc.num) 
scRNA3 <- FindClusters(scRNA3, resolution = 0.5)
table(scRNA3@meta.data$seurat_clusters)
metadata <- scRNA3@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster2/cell_cluster.csv',row.names = F)

##非线性降维
#tSNE
scRNA3 = RunTSNE(scRNA3, dims = pc.num)
embed_tsne <- Embeddings(scRNA3, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster2/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNA3, reduction = "tsne", label=T) 
ggsave("cluster2/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA3, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster2/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster2/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
scRNA3 <- RunUMAP(scRNA3, dims = pc.num)
embed_umap <- Embeddings(scRNA3, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster2/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA3, reduction = "umap", label=T) 
ggsave("cluster2/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNA3, reduction = "umap", group.by='orig.ident')
ggsave("cluster2/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster2/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)
setwd("D:/aJMML/JMML")
#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster2/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
ggsave()
