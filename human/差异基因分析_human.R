rm(list=ls())
#DEGs
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
load("D:/aJMML/scRNA2/sc_seurat_integr_HSC.RData")
DimPlot(HSC)
HSC <- FindMarkers(HSC, min.pct = 0.25, 
                   logfc.threshold = 0.25,
                   group.by = "orig.ident",
                   ident.1 ="JMML",
                   ident.2="Normal")#结果中的Flodchange为第一种/第二种  
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


tHSC1 = HSC[,HSC@meta.data$seurat_clusters %in% c(1)]
tHSC1 <- FindMarkers(tHSC1, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="JMML",
                     ident.2="Normal")#结果中的Flodchange为第一种/第二种  
write.csv(tHSC1,'tHSC1_deg_human.csv')

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


#MPP
MPP <- subset(HSC, celltype=="MPP")
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
