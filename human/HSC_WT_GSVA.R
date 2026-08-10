rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

load("D:/aJMML/scRNA1/sc_seurat_integr.RData")

celltype = data.frame(ClusterID = 0:13,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0),2] = "HSC"
celltype[celltype$ClusterID %in% c(10,1,5),2] = "MPP"
celltype[celltype$ClusterID %in% c(8,9),2] = 'CMP'
celltype[celltype$ClusterID %in% c(11),2] = 'CLP'
celltype[celltype$ClusterID %in% c(3),2] = 'MEP'
celltype[celltype$ClusterID %in% c(4,12),2] = 'Ery' 
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
Idents(sc_integr) = factor(Idents(sc_integr),levels = c('HSC','MPP','CLP','CMP','MEP','Ery','B','Granulocyte'))
sc_integr$celltype <- Idents(sc_integr)
Idents(sc_integr) <- factor(Idents(sc_integr),levels = rev(levels(Idents(sc_integr))))

HSC <- subset(sc_integr, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.4)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
DimPlot(HSC, label = T,pt.size = 1)
DefaultAssay(HSC) <- "RNA"
HSC <- NormalizeData(HSC)

HSC_N <-subset(x=HSC,subset=orig.ident=="Normal")
HSC_J <-subset(x=HSC,subset=orig.ident=="JMML")
####1.读取目标genneset文件####
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

####2.提取分组平均表达矩阵####
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr_N <- AverageExpression(HSC_N, assays = "RNA", slot = "data")[[1]]
expr_N <- expr[rowSums(expr_N)>0,]  #选取非零基因
expr_N <- as.matrix(expr_N)

expr_J <- AverageExpression(HSC_J, assays = "RNA", slot = "data")[[1]]
expr_J <- expr[rowSums(expr_J)>0,]  #选取非零基因
expr_J <- as.matrix(expr_J)

####3.GSVA富集分析####
# gsva默认开启全部线程计算
gsva_res_N <- gsva(expr_N, genesets, method="ssgsea") 
saveRDS(gsva_res_N, "gsva_res_N.rds")
gsva_df_N <- data.frame(Genesets=rownames(gsva_res_N), gsva_res_N, check.names = F)
write.csv(gsva_df_N, "gsva_res_N.csv", row.names = F)

gsva_res_J <- gsva(expr_J, genesets, method="ssgsea") 
saveRDS(gsva_res_J, "gsva_res_J.rds")
gsva_df_J <- data.frame(Genesets=rownames(gsva_res_J), gsva_res_J, check.names = F)
write.csv(gsva_df_J, "gsva_res_j.csv", row.names = F)
####gsva.df####
pheatmap::pheatmap(gsva_res_N, show_colnames = T, scale = "row")
pheatmap::pheatmap(gsva_res_J, show_colnames = T, scale = "row")

####组间差异####
tHSC_WT0<-subset(HSC_WT, idents = 0)

#准备geneset
genesets <- msigdbr(species = "Mus musculus", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

#准备表达矩阵
expr <- as.data.frame(tHSC_WT0@assays$RNA@data)
expr <- as.matrix(expr)
expr[1:3,1:3]
gsva1<-gsva(expr,genesets,kcdf="Gaussian",method="gsva", parallel.sz=12)

mydata <- t(as.data.frame(gsva1))%>%as.data.frame() 
mydata1 <- data.frame(names = row.names(mydata), mydata)
mydata1$group<-"WT"
a=str_which(mydata1$names,'Kras')
mydata1$group[a]<-"Kras"
b=str_which(mydata1$names,'WT')
mydata1$group[b]<-"WT"
mydata2<-mydata1[,-1]

group_list <- data.frame(sample = rownames(mydata2), group =mydata2$group )
group_list<- factor(mydata2$group)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- rownames(mydata)
contrast.matrix <- makeContrasts('Kras-WT',levels = design)
fit <- lmFit(t(mydata[,-1]),design = design)#不要把mydata改错了
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
alldiff=topTable(fit2,coef = 1,n = Inf)

alldiff <- alldiff[order(abs(alldiff$logFC),decreasing = T),]
plotdata <- t(mydata2[,rownames(alldiff[1:20,])])

library(ComplexHeatmap)
library(circlize)
Group <- c('#F8766c','#008892')
names(Group) <- c('Kras','WT')
Top <- HeatmapAnnotation(Group= factor(mydata2$group),
                         annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                      title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                      ncol=1),
                         border = T,
                         col=list(Group=Group),
                         show_annotation_name = TRUE,
                         annotation_name_side="left",
                         annotation_name_gp = gpar(fontsize = 10))

Heatmap(plotdata,##对基因表达量进行限定,也可以进行转置t后scale进行过滤
        name='Score',#热图值图例名
        top_annotation = Top,#顶部注释  
        cluster_rows = FALSE,#不对行进行聚类
        col=colorRamp2(c(-1,0,1),c("navy", "white", "firebrick3")),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        #show_row_names = FALSE,  
        row_names_gp = gpar(fontsize = 8),
        column_split = c(rep(1,table(mydata2$group)[1]),rep(2,table(mydata2$group)[2])),
        gap = unit(4, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize=10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), 
                                  title_gp = gpar(fontsize = 10, fontface = "bold")))
sample_class<-mydata1$group
sample_class<-as.data.frame(sample_class)
rownames(sample_class)<-mydata1$names
colnames(sample_class)<-"group"
ann_colors=list(group=c(Kras='#F8766c',WT='#008892'))
pheatmap::pheatmap(plotdata, show_colnames = F, scale = "row",
                   annotation_col=sample_class,annotation_colors=ann_colors)

#MPP
#MPP <- readRDS("MPP.classified.rds")
load("D:/aJMML/scRNA1/sce_har_anno.RData")
MPP <- subset(sce.harm, celltype=="MPP")
MPP <- ScaleData(MPP, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
MPP <- FindVariableFeatures(MPP, nfeatures = 4000)
MPP <- RunPCA(MPP, npcs = 50, verbose = FALSE)
MPP <- FindNeighbors(MPP, reduction = "pca", dims = 1:50)
MPP <- FindClusters(MPP, resolution =0.7 )
MPP <- RunUMAP(MPP, reduction = "pca", dims = 1:50)
MPP$seurat_clusters <- MPP@active.ident
DimPlot(MPP, label = T,pt.size = 1)

DefaultAssay(MPP) <- "RNA"
MPP <- NormalizeData(MPP)

####基因集准备####
setwd("Function")
dir.create("GSVA/avg_hallmark")
setwd("./Function/GSVA/avg_hallmark")
# 选择基因集
genesets <- msigdbr(species = "Mus musculus", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

####提取分组平均表达矩阵####
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）

expr <- AverageExpression(MPP, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

####GSVA富集分析####
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "MPP_gsva_res.csv", row.names = F)

####gsva.df####
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row")

expr <- AverageExpression(tMPPo, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

####GSVA富集分析####
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "MPP_gsva_res.csv", row.names = F)




tMPP0<-subset(MPP, idents = 0)
expr <- as.data.frame(tMPP0@assays$RNA@data)

genesets <- msigdbr(species = "Mus musculus", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

expr <- as.matrix(expr)
expr[1:3,1:3]
gsva1<-gsva(expr,genesets,kcdf="Gaussian",method="gsva", parallel.sz=12)
gsva.res <- gsva(expr, genesets, method="gsva") 
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)

mydata <- t(as.data.frame(gsva1))%>%as.data.frame() 
mydata1 <- data.frame(names = row.names(mydata), mydata)
mydata1$group<-"WT"
a=str_which(mydata1$names,'Kras')
mydata1$group[a]<-"Kras"
b=str_which(mydata1$names,'WT')
mydata1$group[b]<-"WT"
mydata2<-mydata1[,-1]

group_list <- data.frame(sample = rownames(mydata2), group =mydata2$group )
group_list<- factor(mydata2$group)
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- rownames(mydata)
contrast.matrix <- makeContrasts('Kras-WT',levels = design)
fit <- lmFit(t(mydata[,-1]),design = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
alldiff=topTable(fit2,coef = 1,n = Inf)

alldiff <- alldiff[order(abs(alldiff$logFC),decreasing = T),]
plotdata <- t(mydata2[,rownames(alldiff[1:20,])])

library(ComplexHeatmap)
library(circlize)
Group <- c('#D2691E','#DDA0DD')
names(Group) <- c('Kras','WT')
Top <- HeatmapAnnotation(Group= factor(mydata2$group),
                         annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                      title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                      ncol=1),
                         border = T,
                         col=list(Group=Group),
                         show_annotation_name = TRUE,
                         annotation_name_side="left",
                         annotation_name_gp = gpar(fontsize = 10))

Heatmap(plotdata,##对基因表达量进行限定,也可以进行转置t后scale进行过滤
        name='Score',#热图值图例名
        top_annotation = Top,#顶部注释  
        cluster_rows = FALSE,#不对行进行聚类
        col=colorRamp2(c(-1,0,1),c('#008B8B','#F5F5F5','#DC143C')),#颜色'steelblue','white','firebrick3'蓝白红
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        #show_row_names = FALSE,  
        row_names_gp = gpar(fontsize = 8),
        column_split = c(rep(1,table(mydata2$group)[1]),rep(2,table(mydata2$group)[2])),
        gap = unit(4, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize=10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), 
                                  title_gp = gpar(fontsize = 10, fontface = "bold")))
