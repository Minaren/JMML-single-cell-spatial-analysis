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
celltype[celltype$ClusterID %in% c(0,1,3,4,5,6),2] = "Bright_HSC"
celltype[celltype$ClusterID %in% c(2),2] = "Dark_HSC"
head(celltype)
celltype 
table(celltype$celltype)
HSC$celltype = "NA"
for(i in 1:nrow(celltype)){
  HSC@meta.data[which(HSC@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(HSC@meta.data$celltype)

Idents(HSC) = HSC$celltype
Idents(HSC) = factor(Idents(HSC),levels = c('Bright_HSC','Dark_HSC'))
HSC$celltype <- Idents(HSC)
Idents(HSC) <- factor(Idents(HSC),levels = rev(levels(Idents(HSC))))

DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Cd69","Cd34","Trib2"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69","Cd34"),pt.size=0)
VlnPlot(HSC, features = c("Cd69","Cd34","Cd55","Cd59a"),pt.size=0,ncol = 2)
VlnPlot(HSC, features = c("Cd74","Cd34"),pt.size=0) 
VlnPlot(HSC, features = c("Cd69","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                          "Cd34","Cd63",'Cd9',"Cd117","Cd48","Cd74"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Gfi1b","Erg1","Gata2","Gata1","Meisl",
                          "Hox9"),pt.size=0,split.by = "orig.ident",ncol = 3)#PMID:34724565

ggsave("Cd69_Cd34human.pdf",width = 40,height = 20,units = "cm")
VlnPlot(HSC, features = c("Cd69"),pt.size=0,split.by = "orig.ident")
VlnPlot(sce.harm, features = c("Cd69"),pt.size=0)
FeaturePlot(sce.harm,features = c("Cd69"))
FeaturePlot(sc_integr,features = c("CD69"), min.cutoff =0.5, max.cutoff = 10)
FeaturePlot(HSC,features = c("CD69"),min.cutoff =3, max.cutoff = 10)
