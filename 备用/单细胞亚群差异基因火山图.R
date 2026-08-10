library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db) ##加载小鼠
library(org.Hs.eg.db) ##加载人类
library(tidyverse)
sce_monocle<-readRDS("D:/aJMML/scRNA/monocle/sce_monocle.rds")
pbmc <-sce_monocle
table(pbmc$celltype)
object.markers <- FindMarkers(pbmc, ident.1 = 'HSC',ident.2 = 'MPP', 
                              group.by = 'celltype',logfc.threshold = 0,min.pct = 0,pseudocount.use = 0.01)
object.markers$names <- rownames(object.markers)
#sig_dge.all <- subset(object.markers, p_val_adj<0.05&abs(avg_log2FC)>0.15) #所有差异基因
#View(sig_dge.all)
library(dplyr)
object.markers <- object.markers %>%
  mutate(Difference = pct.1 - pct.2)
library(ggplot)
library(ggrepel)
object.markers$group=0
for (i in 1:nrow(object.markers)){
  if (object.markers$avg_log2FC[i] >= 1 & object.markers$Difference[i] >= 0.2 & object.markers$pct.2[i] <= 0.05){
    object.markers$group[i]='up'
  }
  else if(object.markers$avg_log2FC[i] <= -1 & object.markers$Difference[i] <= -0.2 & object.markers$pct.1[i] <= 0.05){
    object.markers$group[i]='down'
  }
  else {
    object.markers$group[i]='no'
  }
}

ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=0.5,aes(color=group)) + 
  scale_color_manual(values=c('blue','grey','red'))+
  geom_label_repel(data=subset(object.markers, group !='no'), aes(label=names), segment.size = 0.25, size=2.5)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()
ggsave("TopMarkerVol2.pdf", height=8, width=8)
