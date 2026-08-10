library(Seurat) 
library(harmony) 
library(tidyverse) 
library(patchwork)
load("D:/aJMML/scRNA/sce_harm_tu.RData")
genes_to_check =c("Cd3e","Cd4","Cd8a","Ptprc","Ly6g","Itgam","Il7r","Cd19","Ly6a",
                  "Kit","Slamf1","Cd48","Fgfr1","Cd34","Fcgr3","Fcgr2b",
                  "Itga2b","Klrb1c","Itgam","Ly6g","Adgre1","Csf1r","Cd3e","Cd27","Ms4a1",
                  "Cd44","Tfrc","Ly76",
                 "CD34","Flt3","Gcnt2","Hlf","Procr")
DefaultAssay(sce.harm) = "RNA"
DotPlot(sce.harm,features = unique(genes_to_check)) + coord_flip()
ggsave("har_gene_show.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show.png",width = 14,height = 12,dpi = 500)


sce <- sce.harm
celltype = data.frame(ClusterID = 0:22,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,2,5,8,9,13),2] = "Granolucyte"
celltype[celltype$ClusterID %in% c(14),2] = "Macrophage"
celltype[celltype$ClusterID %in% c(18),2] = 'Monocyte'
celltype[celltype$ClusterID %in% c(11,22),2] = 'CLP'
celltype[celltype$ClusterID %in% c(7,17),2] = 'LMPP'
celltype[celltype$ClusterID %in% c(12,21),2] = 'Ery' 
celltype[celltype$ClusterID %in% c(10),2] = 'Ery_prog'
celltype[celltype$ClusterID %in% c(15,6),2] = "HSC"
celltype[celltype$ClusterID %in% c(1),2] = "MPP"
celltype[celltype$ClusterID %in% c(16),2] = "CMP"
celltype[celltype$ClusterID %in% c(19),2] = "NK"
celltype[celltype$ClusterID %in% c(3,4),2] = "CMP"
celltype[celltype$ClusterID %in% c(20),2] = "MK"
 


head(celltype)
celltype 
table(celltype$celltype)
sce.harm$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.harm@meta.data[which(sce.harm@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.harm@meta.data$celltype)

Idents(sce.harm) = sce.harm$celltype
Idents(sce.harm) = factor(Idents(sce.harm),levels = c('HSC','MPP',"LMPP","CMP",'CLP','Ery_prog','Ery',
                                                      'Granolucyte', "Macrophage","NK","Monocyte",'MK'))
sce.harm$celltype <- Idents(sce.harm)
Idents(sce.harm) <- factor(Idents(sce.harm),levels = rev(levels(Idents(sce.harm))))

####4.3 看一下注释后的marker####
DotPlot(sce.harm,features = unique(genes_to_check),cols = c("grey","blue")) + theme_bw(base_line_size = 0) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid = element_blank()) + labs(x='',y='')
ggsave("dotplot_anno.png",width = 19,height = 8,dpi = 500)
ggsave("dotplot_anno.pdf",width = 19,height = 8,dpi = 500)

library(RColorBrewer)
color_ct=c(brewer.pal(12, "Set3"),"#b3b3b3",
           brewer.pal(5, "Set1"),
           brewer.pal(3, "Dark2"),
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")

####4.4 注释后可视化####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#devtools::install_github("junjunlab/jjAnno")

library(ComplexHeatmap)
library(jjAnno)
#if(!require(ggunchull))devtools::install_github("sajuukLyu/ggunchull", type = "source")

#devtools::install_local("D:/R/R-4.2.1/library/scRNAtoolVis/scRNAtoolVis-master.zip")
library(scRNAtoolVis)

## umap/tsne
clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'celltype',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(color_ct,0.65)) + NoLegend() +
  scale_fill_manual(values = alpha(color_ct,0.65))
clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'celltype',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) 
ggsave("umap.pdf",width = 10,height = 9,dpi = 500)
ggsave("umap.png",width = 10,height = 9,dpi = 500)


clusterCornerAxes(object = sce.harm,reduction = 'tsne',clusterCol = 'celltype',cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(color_ct,0.65)) +  NoLegend() +
  scale_fill_manual(values = alpha(color_ct,0.65))
ggsave("tsne.pdf",width = 10,height = 9,dpi = 500)
ggsave("tsne.png",width = 10,height = 9,dpi = 500)


clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("umap_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("umap_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sce.harm,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("tsne_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("tsne_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sce.harm,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("tsne_indi_s.pdf",width = 8,height = 4,dpi = 500)
ggsave("tsne_indi_s.png",width = 8,height = 4,dpi = 500)


clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("umap_indi_s.pdf",width = 8,height = 4,dpi = 500)
ggsave("umap_indi_s.png",width = 8,height = 4,dpi = 500)


save(sce.harm,file = 'sce_har_anno.RData')

####4.5 注释后查看细胞数####
table(sce.harm$celltype)
as.data.frame(table(sce.harm$orig.ident,sce.harm$celltype))
