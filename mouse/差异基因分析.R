rm(list=ls())
library(Seurat)
library(dplyr)
load("D:/aJMML/scRNA1/sce_har_anno.RData")
sce.harm
DimPlot(sce.harm)
ggsave("sce.harm.pdf",width = 10, height = 8)
names(sce.harm@meta.data)
unique(sce.harm$orig.ident)
DimPlot(sce.harm,split.by = 'orig.ident')
sce.harm$celltype.group <- paste(sce.harm$celltype, sce.harm$orig.ident, sep = "_")
sce.harm$celltype <- Idents(sce.harm)
Idents(sce.harm) <- "orig.ident"
mydeg <- FindMarkers(sce.harm,ident.1 = 'Kras',ident.2 = 'WT', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)

cellfordeg<-levels(sce.harm$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(sce.harm,ident.1 = "Kras", ident.2 ="WT", verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()
library(dplyr)
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10
sce.harm <- ScaleData(sce.harm, features =  rownames(sce.harm))
DoHeatmap(sce.harm,features = top10,size=3)
ggsave("Kras_WT_heatmap.pdf",width = 10, height = 8)

Idents(sce.harm) <- "celltype"
VlnPlot(sce.harm,features = top10,split.by = 'orig.ident',idents = 'HSC')
ggsave("Kras_WT_VlnPlot.pdf",width = 10, height = 8)

FeaturePlot(sce.harm,features = top10,split.by = 'orig.ident')
DotPlot(sce.harm,features = top10,split.by ='orig.ident')#默认只有两种颜色
ggsave("Kras_WT_Dotplot.pdf",width = 12, height = 8)

#DotPlot(sce.harm,features = top10,split.by ='orig.ident',cols = c('blue','yellow','pink'))
####提取表达量#######
mymatrix <- as.data.frame(sce.harm@assays$RNA@data)
mymatrix2<-t(mymatrix)%>%as.data.frame()
mymatrix2[,1]<-sce.harm$celltype
colnames(mymatrix2)[1] <- "celltype"

mymatrix2[,ncol(mymatrix2)+1]<-sce.harm$orig.ident
colnames(mymatrix2)[ncol(mymatrix2)] <- "orig.ident"

#绘图
library(ggplot2)
p1<- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=Cd48,fill=orig.ident))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p1
ggsave("Kras_WT_ggplot.pdf",plot=p1,width = 12, height = 8)


DimPlot(object = scedata,label = T,pt.size = 1)+
  labs(x = "UMAP1", y = "UMAP2",title = 'celltype') + 
  theme(legend.position = c(14,0),
        legend.justification = c(0,1),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))



#DEGs
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- FindMarkers(HSC, min.pct = 0.25, 
                               logfc.threshold = 0.25,
                               group.by = "orig.ident",
                               ident.1 ="Kras",
                               ident.2="WT")#结果中的Flodchange为第一种/第二种  
class(HSC)
write.csv(HSC,'HSC_deg.csv')

EnhancedVolcano(HSC,
                lab = rownames(HSC),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'HSC_deg')
ggsave("HSC_deg.pdf",width = 9, height = 8)

vHSC= HSC[,HSC@meta.data$seurat_clusters %in% c(0,2,3,5)]
vHSC <- FindMarkers(vHSC, min.pct = 0.25, 
                   logfc.threshold = 0.25,
                   group.by = "orig.ident",
                   ident.1 ="Kras",
                   ident.2="WT")#结果中的Flodchange为第一种/第二种  
class(HSC)
write.csv(vHSC,'vHSC_deg.csv')

EnhancedVolcano(vHSC,
                lab = rownames(vHSC),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'vHSC_deg')
ggsave("vHSC_deg.pdf",width = 9, height = 8)

tHSC0 = HSC[,HSC@meta.data$seurat_clusters %in% c(0)]
tHSC0 <- FindMarkers(tHSC0, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tHSC0,'tHSC0_deg.csv')

EnhancedVolcano(tHSC0,
                lab = rownames(tHSC0),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tHSC0_deg')
ggsave("tHSC0_deg.pdf",width = 9, height = 8)
##1群
tHSC1 = HSC[,HSC@meta.data$seurat_clusters %in% c(1)]
tHSC1 <- FindMarkers(tHSC1, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tHSC1,'tHSC1_deg.csv')

EnhancedVolcano(tHSC1,
                lab = rownames(tHSC1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tHSC1_deg')
ggsave("tHSC1_deg.pdf",width = 9, height = 8)
##2群
tHSC2 = HSC[,HSC@meta.data$seurat_clusters %in% c(2)]
tHSC2 <- FindMarkers(tHSC2, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tHSC2,'tHSC2_deg.csv')

EnhancedVolcano(tHSC2,
                lab = rownames(tHSC2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tHSC2_deg')
ggsave("tHSC2_deg.pdf",width = 9, height = 8)






#MPP
MPP <- subset(sce.harm, celltype=="MPP")
MPP <- FindMarkers(MPP, min.pct = 0.25, 
                   logfc.threshold = 0.25,
                   group.by = "orig.ident",
                   ident.1 ="Kras",
                   ident.2="WT")#结果中的Flodchange为第一种/第二种  
class(MPP)
write.csv(MPP,'MPP_deg.csv')

EnhancedVolcano(MPP,
                lab = rownames(MPP),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'MPP_deg')
ggsave("MPP_deg.pdf",width = 9, height = 8)

allmarkers <- FindAllMarkers(MPP,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'MPP_allmarkers.csv')







#MPP亚群再分析
##0群
tMPP0 = MPP[,MPP@meta.data$seurat_clusters %in% c(0)]
tMPP0 <- FindMarkers(tMPP0, min.pct = 0.25, 
                   logfc.threshold = 0.25,
                   group.by = "orig.ident",
                   ident.1 ="Kras",
                   ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tMPP0,'MPP_deg.csv')

EnhancedVolcano(tMPP0,
                lab = rownames(tMPP0),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tMPP0_deg')
ggsave("tMPP0_deg.pdf",width = 9, height = 8)
##1群
tMPP1 = MPP[,MPP@meta.data$seurat_clusters %in% c(1)]
tMPP1 <- FindMarkers(tMPP1, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tMPP1,'MPP_deg.csv')

EnhancedVolcano(tMPP1,
                lab = rownames(tMPP1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tMPP1_deg')
ggsave("tMPP1_deg.pdf",width = 9, height = 8)
##2群
tMPP2 = MPP[,MPP@meta.data$seurat_clusters %in% c(2)]
tMPP2 <- FindMarkers(tMPP2, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")#结果中的Flodchange为第一种/第二种  
write.csv(tMPP2,'MPP_deg.csv')

EnhancedVolcano(tMPP2,
                lab = rownames(tMPP2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'tMPP2_deg')
ggsave("tMPP2_deg.pdf",width = 9, height = 8)
