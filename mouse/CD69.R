####素有类型细胞比例####
library(Seurat) 
library(tidyverse) 
load("D:/aJMML/scRNA1/sce_har_anno.RData")
table(sce.harm$orig.ident)#查看各组细胞数
prop.table(table(Idents(sce.harm)))
table(Idents(sce.harm), sce.harm$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(sce.harm), sce.harm$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)


colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage.png',width = 7,height = 6,dpi = 500)
ggsave('percentage.pdf',width = 7,height = 6,dpi = 500)


library(gplots)
tab.1=table(sce.harm$orig.ident,sce.harm$celltype) 
balloonplot(tab.1)
ggsave('balloonplot.png',width = 7,height = 6,dpi = 500)



####HSC亚群再分析####
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.5)
#HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
library(clustree)
clustree(HSC)
Idents(HSC) <- "RNA_snn_res.0.5"
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)

celltype = data.frame(ClusterID = 0:6,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,2,3,5),2] = "De_HSC"
celltype[celltype$ClusterID %in% c(1,4,6),2] = "In_HSC"
celltype[celltype$ClusterID %in% c(6),2] = "Un_HSC"
head(celltype)
celltype 
table(celltype$celltype)
HSC$celltype = "NA"
for(i in 1:nrow(celltype)){
  HSC@meta.data[which(HSC@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(HSC@meta.data$celltype)

Idents(HSC) = HSC$celltype
Idents(HSC) = factor(Idents(HSC),levels = c("De_HSC","In_HSC","Un_HSC"))
HSC$celltype <- Idents(HSC)
Idents(HSC) <- factor(Idents(HSC),levels = rev(levels(Idents(HSC))))

DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Cd69","Kras","Cd63","Cd52"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Kras","Mapk1"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69","Cd34"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69","Cd34","Kras","Cd59a"),pt.size=0,ncol = 2)
VlnPlot(HSC, features = c("Cd74","Cd34"),pt.size=0)
VlnPlot(HSC, features = c("Cd69","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a",
                          "Cd34","Cd63",'Cd9',"Cd117","Kras"),
        stacked=T,pt.size=0)
VlnPlot(HSC, features = c("Cd69"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a",
                          "Cd34","Cd63",'Cd9',"Cd48","Cd74","Kras"),pt.size=0,
        split.by = "orig.ident")
ggsave("3.pdf",width = 40,height = 40,units = "cm")


HSC_Kras <- subset(HSC, orig.ident=="Kras")
HSC_WT <- subset(HSC, orig.ident=="WT")
VlnPlot(HSC_Kras, features = c("Cd69"),pt.size=0)+geom_boxplot(width=.2,col="black",fill="white")
VlnPlot(HSC_WT, features = c("Cd69"),pt.size=0)+geom_boxplot(width=.2,col="black",fill="white")
