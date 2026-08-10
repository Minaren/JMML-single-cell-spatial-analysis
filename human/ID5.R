
####准备工作####
rm(list = ls())
library(Seurat) 
library(tidyverse) 
library(patchwork)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))

counts1 <- Read10X(data.dir = "D:/aJMML/JMML/JMMLID5") 
ID5<-CreateSeuratObject(counts = counts1, project = "CB35", min.cells = 3,min.features=200
)
ID5[["percent.MT"]]<-PercentageFeatureSet(ID5,pattern = "MT")
ID5<-subset(ID5,subset = nFeature_RNA>500&nFeature_RNA<6000&percent.MT<10)
VlnPlot(ID5,features = c("nFeature_RNA","nCount_RNA","percent.MT"),ncol = 3)
ID5<-NormalizeData(ID5)
ID5<- FindVariableFeatures(ID5, nfeatures = 2000)
ID5<- ScaleData(ID5)
ID5 <- RunPCA(ID5, npcs = 50, verbose = FALSE)
ElbowPlot(ID5, ndims = 50)
pc.num=1:30
ID5 <- ID5 %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
ID5 <- FindNeighbors(ID5, dims=pc.num) %>% FindClusters()
DimPlot(ID5, label = T)

####9.UMAP/tSNE非线性降维####
#(1)一般认为，umap适合大型数据
ID5<-RunUMAP(ID5,dims = 1:30)
DimPlot(ID5,label = T,pt.size = 1)
#(2)tsne适合小型数据，因为计算起来比较慢，但分离效果优于umap
ID5<-RunTSNE(ID5,dims = 1:30)
DimPlot(ID5,reduction = "tsne")

genes_to_check =c("CD2","ITGAM","ITGAL","ANPEP","CD14","FUT4","CD19","CD33","IL6R","IL7R","MME","CD7","CD127",
                  "TFRC","EPOR","ITGA2B","GP1BA","PTPRC","CD34","CD38","THY1","ITGA6","KIT","FLT3","IL3RA")
DotPlot(ID5,features = unique(genes_to_check),cols = c("blue","red")) + theme_bw(base_line_size = 0) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid = element_blank()) + labs(x='',y='')+ coord_flip()


HSC = ID5[,ID5@meta.data$seurat_clusters %in% c(1)]
HSC <- subset(sce.h, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
library(clustree)
clustree(HSC)
Idents(HSC) <- "RNA_snn_res.0.5"
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)
allmarkers <- FindAllMarkers(HSC,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_HSC_human_ID5.csv') 


VlnPlot(ID5, features = c("CD69"),pt.size=0)


VlnPlot(HSC, features = c("CD69"),pt.size=0,split.by = "orig.ident")
ggsave("Cd69_huamn_split.pdf",width = 40,height = 20,units = "cm")
VlnPlot(HSC, features = c("CD69","CD36"),pt.size=0)
ggsave("Cd69_huamn.pdf",width = 40,height = 40,units = "cm")
FeaturePlot(HSC,features = c("CD69","CD36"),min.cutoff =1, max.cutoff = 10)

VlnPlot(HSC, features = c("CD69","EGR1","GATA2","JUNB",'SOX4',"FTL
"),pt.size=0)
ggsave("Cd69_CD36_ID5.pdf",width = 40,height = 20,units = "cm")
