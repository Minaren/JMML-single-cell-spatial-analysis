rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
####构建参考数据集####
rm(list = ls())
sc<-readRDS("D:/aJMML/scRNA_mouse_R0.6_6/GSE147729/notlabel.RDS")
str(sc)
slotNames(sc)
data_matrix <- read.table("metaInfo.txt", header = TRUE, row.names = 1, sep = "\t")
head(data_matrix)
# 提取 CellType 列
celltype_data <- data_matrix[rownames(sc@meta.data), "CellType..based.on.scmap."]
# 将 CellType 添加到 Seurat 对象的 meta.data 中
sc@meta.data$CellType <- celltype_data
head(sc@meta.data)
sc.updated = UpdateSeuratObject(object = sc)
saveRDS(sc, file = "healthy_mouse_参考数据集.rds")

####构建询问数据集####
load("D:/aJMML/scRNA1/sce_harm_tu.RData")
anchors <- FindTransferAnchors(reference = sc.updated, query = sce.harm, dims = 1:30)
refdata <- sc.updated$CellType
predictions <- TransferData(anchorset = anchors, refdata = refdata, dims = 1:30)
sc_q <- AddMetaData(sce.harm, metadata = predictions)
sc_q$celltype <- sc_q$predicted.id
sc_q$celltype <- as.factor(sc_q$celltype)
head(sc_q$celltype)
table(sc_q$celltype)
DimPlot(sc_q, reduction = "umap",group.by = "celltype", label = TRUE)
ggsave('umap_Kras_注释后.png',width = 6,height = 6,dpi = 500)
ggsave('umap_sc_q_注释后.pdf',width = 20,height = 20,dpi = 500)


sc_WT<- readRDS("sc_WT_注释前.rds")
anchors <- FindTransferAnchors(reference = sc.updated, query = sc_WT, dims = 1:30)
refdata <- sc.updated$CellType
predictions <- TransferData(anchorset = anchors, refdata = refdata, dims = 1:30)
sc_q <- AddMetaData(sc_WT, metadata = predictions)
sc_q$celltype <- sc_q$predicted.id
sc_q$celltype <- as.factor(sc_q$celltype)
head(sc_q$celltype)
table(sc_q$celltype)
DimPlot(sc_q, reduction = "umap",group.by = "celltype", label = TRUE)
ggsave('umap_Kras_注释后.png',width = 6,height = 6,dpi = 500)
ggsave('umap_sc_WT_注释后.pdf',width = 20,height = 20,dpi = 500)

