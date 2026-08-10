#单细胞教程|Seurat v5分析流程
####1_创建Seurat对象####
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the WT dataset
WT<- Read10X(data.dir = "D:/aJMML/JMML/WT")
sc_WT <- CreateSeuratObject(counts = WT, project = "sc_WT", min.cells = 3, min.features = 200)
sc_WT

####2_数据质控####
# Visualize QC metrics as a violin plot
sc_WT[["percent.mt"]] <- PercentageFeatureSet(sc_WT, pattern = "^mt-")
VlnPlot(sc_WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#这里和之前分析保持一致
sc_WT <- subset(sc_WT, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & percent.mt < 10)
dim(sc_WT)
####3_数据标准化####
sc_WT <- NormalizeData(sc_WT, normalization.method = "LogNormalize", scale.factor = 10000)
sc_WT <- NormalizeData(sc_WT)
sc_WT <- FindVariableFeatures(sc_WT, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_WT), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc_WT)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


####4_数据归一化####
all.genes <- rownames(sc_WT)
sc_WT <- ScaleData(sc_WT, features = all.genes)
#删除不需要的变异源
sc_WT <- ScaleData(sc_WT, vars.to.regress = "percent.mt")

####5_进行PCA线性降维####
sc_WT <- RunPCA(sc_WT, features = VariableFeatures(object = sc_WT))
# Examine and visualize PCA results a few different ways
print(sc_WT[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc_WT, dims = 1:2, reduction = "pca")
DimPlot(sc_WT, reduction = "pca") + NoLegend()
DimHeatmap(sc_WT, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc_WT, dims = 1:15, cells = 500, balanced = TRUE)

####6_确定维度####
ElbowPlot(sc_WT,ndims=50)

####7_数据聚类####
sc_WT <- FindNeighbors(sc_WT, dims = 1:30)
sc_WT <- FindClusters(sc_WT, resolution = 0.6)

####8_降维tsen/umap####
sc_WT <- RunUMAP(sc_WT, dims = 1:30)
DimPlot(sc_WT, reduction = "umap",label=TRUE)
ggsave('umap_注释前.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap_注释前.png',width = 6,height = 6,dpi = 500)
saveRDS(sc_WT, file = "sc_WT_注释前.rds")

####9_寻找标记基因####
sc_WT.markers <- FindAllMarkers(sc_WT, only.pos = TRUE)
write.csv(sc_WT.markers,'allmarkers_mouse.csv')
sc_WT.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
genes_to_check =c("Prss34","Cd200r3","Lmo4",
                  "Vpreb3","Cd79a","Cd79b","Cd74","Cd19",
                  "Prg2","Itga2b","Plek","Cebpe","Lmo4","Prss34",
                  "Blvrb","Rhd","Hmbs","Gata1","Car1",
                  "Hba-a1","Hbb-bt","Hbb-bs","Blvrb","Hba-a2","Hmbs",
                  "Cd68","Cd74","Mpeg1","Adgre1","Csf1r","C1qc",
                  "Cd74","Mpeg1","Cd68","Csf1r","Adgre1",
                  "Csf1r","Cd74","Mpeg1","Cd68","S100a4",
                  "Pf4","Ctla2a","Itga2b","Gata2","Plek","Cd9",
                  "Prtn3","Elane","Mpo","Ctsg",
                  "Elane","Mpo","Prtn3",
                  "Ms4a6c","F13a1","Irf8","Csf1r",
                  "Cd34","Flt3","Ctla2a",
                  "Camp","Ngp","Lcn2","Cebpe","Fcnb",
                  "Mpo","Cd177",
                  "Lcn2",
                  "Ngp","Lcn2","Cd177","Cebpe","Ltf",
                  "Ngp","Lcn2","Cd177","Cebpe",
                  "Mmp8","Lcn2","Ngp","Camp","Retnlg",
                  "Lmo4","Runx1",
                  "Slamf1","Kit","Cd34","Sca1","Cd48","Procr","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf","Cd201","Cd47",
                  "Ms4a1",
                  "Il7r,Ccr7,Cd8a")
DefaultAssay(sc_WT) = "RNA"
DotPlot(sc_WT,features = unique(genes_to_check)) + coord_flip()
ggsave("genes_to_check.pdf",width = 14,height = 12,dpi = 500)
ggsave("genes_to_check.png",width = 14,height = 12,dpi = 500)

surface_markers =c("Kit","Ly6a","Flt3","Cd34",
                   "Slamf1","Cd48","Cd16","Cd32","Il7","Tfr","Cd41","Cd11b","Ly6g","Eporr","Cd115")
DotPlot(sc_WT,features = unique(surface_markers)) + coord_flip()
ggsave("surface_markers.pdf",width = 14,height = 12,dpi = 500)
ggsave("surface_markers.png",width = 14,height = 12,dpi = 500)

####10_细胞注释####
sc_WT<- readRDS("sc_WT_注释前.rds")
sc_anno<- sc_WT
celltype = data.frame(ClusterID = 0:19,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,1,4,5,7,14,17),2] = "Granolucyte"
celltype[celltype$ClusterID %in% c(12,19),2] = "Macrophage"
celltype[celltype$ClusterID %in% c(19),2] = "Monocyte"
celltype[celltype$ClusterID %in% c(15),2] = 'CLP'
celltype[celltype$ClusterID %in% c(6),2] = 'MPP'
celltype[celltype$ClusterID %in% c(8,9,16),2] = 'Erythroid' 
celltype[celltype$ClusterID %in% c(13),2] = "HSC"
celltype[celltype$ClusterID %in% c(3),2] = "LMPP"
celltype[celltype$ClusterID %in% c(2),2] = "CMP"
celltype[celltype$ClusterID %in% c(10),2] = "T cell"
celltype[celltype$ClusterID %in% c(18),2] = "B cell"
celltype[celltype$ClusterID %in% c(11),2] = "Megakaryocyte"

head(celltype)
celltype 
table(celltype$celltype)
sc_anno$celltype = "NA"
for(i in 1:nrow(celltype)){
  sc_anno@meta.data[which(sc_anno@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sc_anno@meta.data$celltype)

Idents(sc_anno) = sc_anno$celltype
Idents(sc_anno) = factor(Idents(sc_anno),levels = c('HSC','MPP',"LMPP","CMP",'CLP','Erythroid',
                                                      'Granolucyte', "Macrophage","Monocyte","T cell","B cell",'Megakaryocyte'))
sc_anno$celltype <- Idents(sc_anno)
Idents(sc_anno) <- factor(Idents(sc_anno),levels = rev(levels(Idents(sc_anno))))

#检查注释后marker
DotPlot(sc_anno,features = unique(genes_to_check)) + coord_flip()
ggsave("genes_to_check_注释后.pdf",width = 14,height = 12,dpi = 500)
ggsave("genes_to_check_注释后.png",width = 14,height = 12,dpi = 500)
saveRDS(sc_anno, file = "sc_WT_注释后.rds")

####11_验证注释可行性####
rm(list = ls())
library(msigdbr)
library(GSVA)
#BiocManager::install("GSVA")
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
library(Seurat)
library(ggplot2)

sc_anno<- readRDS("sc_WT_注释后.rds")
#读取目标genneset文件
genesets <- read.csv("gsva_mouse_cluster.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

#提取分组平均表达矩阵#
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(sc_anno, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

#GSVA富集分析#
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "sc_WT_ssgsea.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res_mouse.csv", row.names = F)
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_row = FALSE,cluster_cols = FALSE,filename = "ssGSEA验证注释可行性.png")
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row",cluster_row = FALSE,cluster_cols = FALSE,filename = "ssGSEA验证注释可行性.pdf")

####12_构建询问数据集####
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Load the Kras dataset
Kras<- Read10X(data.dir = "D:/aJMML/JMML/Kras")
sc_Kras <- CreateSeuratObject(counts = Kras, project = "sc_Kras", min.cells = 3, min.features = 200)
sc_Kras

#数据质控
# Visualize QC metrics as a violin plot
sc_Kras[["percent.mt"]] <- PercentageFeatureSet(sc_Kras, pattern = "^mt-")
VlnPlot(sc_Kras, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#这里和之前分析保持一致
sc_Kras <- subset(sc_Kras, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & percent.mt < 10)
dim(sc_Kras)
#数据标准化#
sc_Kras <- NormalizeData(sc_Kras, normalization.method = "LogNormalize", scale.factor = 10000)
sc_Kras <- NormalizeData(sc_Kras)
sc_Kras <- FindVariableFeatures(sc_Kras, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_Kras), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc_Kras)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#数据归一化#
all.genes <- rownames(sc_Kras)
sc_Kras <- ScaleData(sc_Kras, features = all.genes)
#删除不需要的变异源
sc_Kras <- ScaleData(sc_Kras, vars.to.regress = "percent.mt")

#进行PCA线性降维#
sc_Kras <- RunPCA(sc_Kras, features = VariableFeatures(object = sc_Kras))
# Examine and visualize PCA results a few different ways
print(sc_Kras[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc_Kras, dims = 1:2, reduction = "pca")
DimPlot(sc_Kras, reduction = "pca") + NoLegend()
DimHeatmap(sc_Kras, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc_Kras, dims = 1:15, cells = 500, balanced = TRUE)

#确定维度#
ElbowPlot(sc_Kras,ndims=50)

#数据聚类#
sc_Kras <- FindNeighbors(sc_Kras, dims = 1:30)
sc_Kras <- FindClusters(sc_Kras, resolution = 0.6)

#降维tsen/umap#
sc_Kras <- RunUMAP(sc_Kras, dims = 1:30)
DimPlot(sc_Kras, reduction = "umap",label=TRUE)
ggsave('umap_Kras_注释前.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap_Kras_注释前.png',width = 6,height = 6,dpi = 500)
saveRDS(sc_Kras, file = "sc_Kras_注释前.rds")

####13_Transferdata数据映射####
rm(list = ls())
library(Seurat)
sc_Kras<- readRDS("sc_Kras_注释前.rds")
sc_anno<- readRDS("sc_WT_注释后.rds")
anchors <- FindTransferAnchors(reference = sc_anno, query = sc_Kras, dims = 1:30)
refdata <- sc_anno$celltype
predictions <- TransferData(anchorset = anchors,refdata = refdata, dims = 1:30)
sc_Kras <- AddMetaData(sc_Kras, metadata = predictions)
sc_Kras$celltype <- sc_Kras$predicted.id
sc_Kras$celltype <- as.factor(sc_Kras$celltype)
head(sc_Kras$celltype)
table(sc_Kras$celltype)
DimPlot(sc_Kras, reduction = "umap",group.by = "celltype", label = TRUE)
ggsave('umap_Kras_注释后.png',width = 6,height = 6,dpi = 500)
ggsave('umap_Kras_注释后.pdf',width = 20,height = 20,dpi = 500)
