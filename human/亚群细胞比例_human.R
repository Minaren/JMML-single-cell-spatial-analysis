####素有类型细胞比例####
library(Seurat) 
load("D:/aJMML/scRNA1/sc_seurat_integr.RData")

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
Idents(sc_integr) = factor(Idents(sc_integr),levels = c('HSC','MPP',"LMPP",'CLP','CMP','MEP','Ery','B','Granulocyte'))
sc_integr$celltype <- Idents(sc_integr)
Idents(sc_integr) <- factor(Idents(sc_integr),levels = rev(levels(Idents(sc_integr))))


table(sc_integr$orig.ident)#查看各组细胞数
prop.table(table(Idents(sc_integr)))
table(Idents(sc_integr), sc_integr$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(sc_integr), sc_integr$orig.ident), margin = 2)#计算各组样本不同细胞群比例
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
tab.1=table(sc_integr$orig.ident,sc_integr$celltype) 
balloonplot(tab.1)
ggsave('balloonplot.png',width = 7,height = 6,dpi = 500)

####HSC亚群再分析####

HSC <- subset(sc_integr, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
#HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0,by = 0.1))
HSC <- FindClusters(HSC, resolution = 0.4)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
#library(clustree)
#clustree(HSC)
HSC$seurat_clusters<- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)

allmarkers <- FindAllMarkers(HSC,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_HSC_human_CD69.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_HSC_human.csv')

p<-DoHeatmap(HSC, features = top10$gene) + NoLegend()
ggsave("markers_HSC_human.heatmap.png", plot = p, width = 17, height = 17)
ggsave("markers_HSC_human.heatmap.pdf", plot = p, width = 17, height = 17)

VlnPlot(HSC, 
        features = c("PTPN","PTPN11","NRAS","KRAS","NF1","CBL","NF1","JUNB"),
        pt.size = 0,
        ncol = 2)
VlnPlot(HSC, 
        features = c("Pf4","Pbx1","Sox4","Junb","Klf1","Aqp1"),
        pt.size = 0,
        ncol = 2,
        split.by = "orig.ident")
VlnPlot(HSC, features = c("CD69"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("CD69"),pt.size=0)

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
ggsave('percentage_HSC_human.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_HSC_human.pdf',width = 7,height = 6,dpi = 500)


colourCount = length(unique(Arthritis$Var1))
Arthritis<-table(Idents(HSC), HSC$orig.ident)
Arthritis<-as.data.frame(Arthritis)
ggplot(Arthritis) + 
  geom_bar(aes(x =Var2, y= Freq,fill=Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Count')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))



####MPP亚群再分析####
MPP <- subset(sc_integr, celltype=="MPP")
MPP <- ScaleData(MPP, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
MPP <- FindVariableFeatures(MPP, nfeatures = 4000)
MPP <- RunPCA(MPP, npcs = 50, verbose = FALSE)
ElbowPlot(MPP, ndims = 50)
MPP <- FindNeighbors(MPP, reduction = "pca", dims = 1:50)
#MPP_r <- FindClusters(MPP, resolution = seq(from = 0.1, to = 1.0,by = 0.1))
MPP <- FindClusters(MPP, 
                    resolution =0.2)
MPP <- RunUMAP(MPP, reduction = "pca", dims = 1:50)
#library(clustree)
#clustree(MPP_r)

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
ggsave('percentage_MPP_human.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_MPP_human.pdf',width = 7,height = 6,dpi = 500)


Arthritis<-table(Idents(MPP), MPP$orig.ident)
Arthritis<-as.data.frame(Arthritis)
colourCount = length(unique(Arthritis$Var1))
ggplot(Arthritis) + 
  geom_bar(aes(x =Var2, y= Freq,fill=Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Count')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('count_MPP_human.pdf',width = 7,height = 6,dpi = 500)


allmarkers <- FindAllMarkers(MPP,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_MPP_human.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_MPP_human.csv')

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


