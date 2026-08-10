rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
load("D:/aJMML/scRNA1/sc_seurat_integr.RData")

####1.读取目标genneset文件####
genesets <- read.csv("human_signature_cluster.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

####2.提取分组平均表达矩阵####
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(sc_integr, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

####3.GSVA富集分析####
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row")


#注释后gsva
celltype = data.frame(ClusterID = 0:13,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0),2] = "HSC"
celltype[celltype$ClusterID %in% c(1),2] = "MPP"
celltype[celltype$ClusterID %in% c(10),2] = "LMPP"
celltype[celltype$ClusterID %in% c(5),2] = "CMP"
celltype[celltype$ClusterID %in% c(8,9),2] = 'GMP'
celltype[celltype$ClusterID %in% c(11),2] = 'CLP'
celltype[celltype$ClusterID %in% c(3),2] = 'MEP'
celltype[celltype$ClusterID %in% c(4),2] = 'Ery_prog' 
celltype[celltype$ClusterID %in% c(12),2] = 'MK' 
celltype[celltype$ClusterID %in% c(13),2] = "B" 
celltype[celltype$ClusterID %in% c(2,6,7),2] = 'Granulocyte'

head(celltype)
celltype 
table(celltype$celltype)
sc_integr$celltype = "NA"
for(i in 1:nrow(celltype)){
  sc_integr@meta.data[which(sc_integr@meta.data$integrated_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sc_integr@meta.data$celltype)

Idents(sc_integr) = sc_integr$celltype
Idents(sc_integr) = factor(Idents(sc_integr),levels = c('HSC','MPP',"LMPP",'CLP','CMP',"GMP",'MEP','Ery_prog',"MK",'B','Granulocyte'))
sc_integr$celltype <- Idents(sc_integr)
Idents(sc_integr) <- factor(Idents(sc_integr),levels = rev(levels(Idents(sc_integr))))

genesets <- read.csv("human_signature_cluster1.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

expr <- AverageExpression(sc_integr, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row")
