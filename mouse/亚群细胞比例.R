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

allmarkers <- FindAllMarkers(HSC,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_HSC_mouse.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_HSC_mouse.csv')

p<-DoHeatmap(HSC, features = top10$gene) + NoLegend()
ggsave("HSC_top10_mouse.png", plot = p, width = 17, height = 17)
ggsave("HSC_top10_mouse.pdf", plot = p, width = 17, height = 17)

DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Hlf","Ifitm1","Flt3","Dntt","Gata1","Car1","Mki67","Elane"
                          ,"Mpo","Pf4","Vwf","Itga2b",'Cebpe'),stacked=T,pt.size=0)
VlnPlot(HSC, features = c("Cd74","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                          "Cd34","Cd63",'Cd9',"Cd117"),
        stacked=T,pt.size=0)


VlnPlot(HSC, 
        features = c("Cd74","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                     "Cd34","Cd63",'Cd9',"Cd117","Cd48"),
        pt.size = 0,
        ncol = 2,
        split.by = "orig.ident")
ggsave("vln9.pdf",width = 20,height = 40,units = "cm")

table(HSC$orig.ident)#查看各组细胞数
prop.table(table(Idents(HSC)))
table(Idents(HSC), HSC$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(HSC), HSC$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'percentage')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage_HSC_r0.5.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_HSC_r0.5.pdf',width = 7,height = 6,dpi = 500)


colourCount = length(unique(Arthritis$Var1))
Arthritis<-table(Idents(HSC), HSC$orig.ident)
Arthritis<-as.data.frame(Arthritis)
ggplot(Arthritis) + 
  geom_bar(aes(x =Var2, y= Freq,fill=Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Count')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('count_HSC_r0.5.png',width = 7,height = 6,dpi = 500)
ggsave('count_HSC_r0.5.pdf',width = 7,height = 6,dpi = 500)


####MPP亚群再分析####
MPP <- subset(sce.harm, celltype=="MPP")
MPP <- ScaleData(MPP, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
MPP <- FindVariableFeatures(MPP, nfeatures = 4000)
MPP <- RunPCA(MPP, npcs = 50, verbose = FALSE)
ElbowPlot(MPP, ndims = 50)
MPP <- FindNeighbors(MPP, reduction = "pca", dims = 1:50)
MPP <- FindClusters(MPP, 
                    resolution =0.3)
MPP <- RunUMAP(MPP, reduction = "pca", dims = 1:50)
library(clustree)
clustree(MPP)
Idents(MPP) <- "RNA_snn_res.0.7"
MPP$seurat_clusters <- MPP@active.ident
DimPlot(MPP, label = T,pt.size = 1)

table(MPP$orig.ident)#查看各组细胞数
prop.table(table(Idents(MPP)))
table(Idents(MPP), MPP$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(MPP), MPP$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'percentage')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage_MPP.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_MPP.pdf',width = 7,height = 6,dpi = 500)


Arthritis<-table(Idents(MPP), MPP$orig.ident)
Arthritis<-as.data.frame(Arthritis)
colourCount = length(unique(Arthritis$Var1))
ggplot(Arthritis) + 
  geom_bar(aes(x =Var2, y= Freq,fill=Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Count')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('count_MPP.pdf',width = 7,height = 6,dpi = 500)


allmarkers <- FindAllMarkers(MPP,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_MPP.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10.csv')

p<-DoHeatmap(MPP, features = top10$gene) + NoLegend()
ggsave("markers_MPP.png", plot = p, width = 17, height = 17)
ggsave("markers_MPP.pdf", plot = p, width = 17, height = 17)

VlnPlot(MPP, 
        features = c("Pf4","Pbx1","Sox4","Junb","Klf1","Aqp1"),
        pt.size = 0,
        ncol = 2)
VlnPlot(MPP, 
        features = c("Lgals3","Cd74","Dusp1","Atf3","Runx2","Bst2","Hmgb2",
                     "Pcna"),
        pt.size = 0,
        ncol = 2,
        split.by = "orig.ident")


