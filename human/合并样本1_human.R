####准备工作####
rm(list = ls())
options(stringsAsFactors = F)
set.seed(220625)
library(Seurat) 
library(harmony) 
library(tidyverse) 
library(patchwork)

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(harmony))


####一、多个单细胞样本的合并####
####1.读取并合并数据==####
####1.1读取数据####
dir = c('D:/aJMML/JMML/JMMLID5',
        'D:/aJMML/JMML/PBM1')
samples_name= c('JMML','Normal')

####1.2 批量创建seurat对象####
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.MT"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  }
}

####1.3 使用merge函数将scRNAlist合成一个Seurat对象####
sc_merge <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

table(sc_merge$orig.ident) 

head(sc_merge@meta.data,5)
save(sc_merge,file = 'sc_merge.Rdata')

####二、数据质控####
####2.1 计算各类基因占比####
# 设置可能用到的主题 
scRNA <- sc_merge
theme.set2 = theme(axis.title.x=element_blank()) 
# 设置绘图元素 
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.MT") 
group = "orig.ident" 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=3) 
ggsave("vlnplot_before_qc2.pdf", plot = violin, width = 9, height = 15)

####2.2 设置质控标准####
minGene=500
maxGene=6000#改成6000
pctMT=10
#maxUMI=15000
#pctMT=10 
#pctHB=1 
####2.3 过滤后可视化####
sc_filt<- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                   percent.MT < pctMT) 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(sc_filt, group.by=group, pt.size = 0, 
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("vlnplot_after_qc2.pdf", plot = violin, width = 10, height = 8)
save(sc_filt, file = "sc_filt.RData")

####三. 数据整合####
####2.1 准备数据####
####2.2 数据标准化（和锚点整合不同，不需拆分样本，直接标准化）####
load("D:/aJMML/scRNA/sc_filt.RData")
sc_n= NormalizeData(sc_filt)
sc_n= FindVariableFeatures(sc_n,selection.method = "vst",nfeatures = 2000, verbose = FALSE)
sc_n = ScaleData(sc_n)
sc_n

####2.3 使用harmony整合数据####
### PCA
sc_n <- RunPCA(sc_n, npcs=50, verbose=FALSE)
ElbowPlot(sc_n, ndims = 50)
ggsave('elbowplot_n2.pdf',width = 6.5,height = 5,dpi = 500)
ggsave('elbowplot_n2.png',width = 6.5,height = 5,dpi = 500)
sce_ha = sc_n
pc.num=1:30
### 整合
sce_har<-RunHarmony(sce_ha,group.by.vars = "orig.ident",project.dim = F,plot_convergence = T)

####2.4 分群####
sce.har <- sce_har
sce.har@active.assay
sce.har = FindNeighbors(sce.har, dims = pc.num,reduction = "harmony")
sce.har = FindClusters(sce.har, graph.name = "RNA_snn", resolution = 0.6, algorithm = 1)
table(sce.har@active.ident)


####2.5 降维及可视化####
sce.harm = RunUMAP(sce.har, dims = pc.num,reduction = "harmony")
save(sce.harm,file = "sce_harm_tu.RData")


DimPlot(sce.harm,label = T,repel = T) + NoLegend()
ggsave('umap2.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap2.png',width = 6,height = 6,dpi = 500)


DimPlot(sce.harm,label = F,group.by = 'orig.ident')
ggsave('umap_indi.pdf',width = 7,height = 6,dpi = 500)
ggsave('umap_indi.png',width = 7,height = 6,dpi = 500)


####四.细胞注释####
####4.1 基因检查####
rm(list = ls())
load("D:/aJMML/scRNA/sce_harm_tu.RData")
#以下使用的基因为张琳琳师姐提供的
genes_to_check =c("Prss34","Cd200r3","Lmo4",
                  "Vpreb3","Cd79a","Cd79b","Cd74","Cd19",
                  "Prg2","Itga2b","Plek","Cebpe","Lmo4","Prss34",
                  "Blvrb","Rhd","Hmbs","Gata1","Car1",
                  "Hba-a1","Hbb-bt","Hbb-bs","Blvrb","Hba-a2","Hmbs",
                  "Cd68","Cd74","Mpeg1","Adgre1","Csf1r","C1qc",
                  "Cd74","Mpeg1","Cd68","Csf1r","Adgre1",
                  "Csf1r","Cd74","Mpeg1","Cd68","S100a4",
                  "Pf4","Ctla2a","Itga2b","Gata2","Plek","Cd9",
                  "Prtn3","Elane","Mpo","Ctsg",
                  "Elane","Mpo","Prtn3",
                  "Ms4a6c","F13a1","Irf8","Csf1r",
                  "Cd34","Flt3","Ctla2a",
                  "Camp","Ngp","Lcn2","Cebpe","Fcnb",
                  "Mpo","Cd177",
                  "Lcn2",
                  "Ngp","Lcn2","Cd177","Cebpe","Ltf",
                  "Ngp","Lcn2","Cd177","Cebpe",
                  "Mmp8","Lcn2","Ngp","Camp","Retnlg",
                  "Lmo4","Runx1",
                  "Slamf1","Kit","Cd34","Sca1","Cd48","Procr","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf","Cd201","Cd47")
DefaultAssay(sce.harm) = "RNA"
DotPlot(sce.harm,features = unique(genes_to_check)) + coord_flip()
ggsave("har_gene_show.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show.png",width = 14,height = 12,dpi = 500)

allmarkers <- FindAllMarkers(sce.harm,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10.csv')

p<-DoHeatmap(sce.harm, features = top10$gene) + NoLegend()
ggsave("markers.heatmap.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap.pdf", plot = p, width = 17, height = 17)

#点图可视化
p <- DotPlot(sce.harm, features = unique(top10$gene) ,
             assay='RNA' )  + coord_flip()
plotd<-p+ theme(axis.text.x = element_text(angle = 45, 
                                           vjust = 0.5, hjust=0.5))
#保存可视化图片
ggsave("markers.dotplot.png", plot = plotd, width = 17, height = 49)
ggsave("markers.dotplot.pdf", plot = plotd, width = 17, height = 49)

####富集####
rm(list=ls())
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)#人的用(org.Hs.eg.db)
allmarkers<-read.csv("D:/aJMML/scRNA/allmarkers.csv")
ids=bitr(allmarkers$gene,'SYMBOL','ENTREZID','org.Mm.eg.db') ## 将SYMBOL转成ENTREZID
allmarkers=merge(allmarkers,ids,by.x='gene',by.y='SYMBOL')
View(allmarkers)

#富集分析
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(allmarkers$ENTREZID, allmarkers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("godotplot1.pdf", plot = p, width = 17, height = 45)



####4.2 注释####
sce <- sce.harm
celltype = data.frame(ClusterID = 0:22,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,2,4,5,8,9),2] = "Neutrophil"
celltype[celltype$ClusterID %in% c(1,13,14,17,18),2] = "Macrophage"
celltype[celltype$ClusterID %in% c(3),2] = 'CMP'
celltype[celltype$ClusterID %in% c(6),2] = 'CLP'
celltype[celltype$ClusterID %in% c(7),2] = 'MPP'
celltype[celltype$ClusterID %in% c(10,21),2] = 'Ery' 
celltype[celltype$ClusterID %in% c(11),2] = "Pre_B" 
celltype[celltype$ClusterID %in% c(12),2] = 'Ery_prog'
celltype[celltype$ClusterID %in% c(15),2] = "HSC"
celltype[celltype$ClusterID %in% c(16),2] = 'Basophil' 
celltype[celltype$ClusterID %in% c(19),2] = "T" 
celltype[celltype$ClusterID %in% c(20),2] = "Unknown"
celltype[celltype$ClusterID %in% c(22),2] = 'B' 


head(celltype)
celltype 
table(celltype$celltype)
sce.harm$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.harm@meta.data[which(sce.harm@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.harm@meta.data$celltype)

Idents(sce.harm) = sce.harm$celltype
Idents(sce.harm) = factor(Idents(sce.harm),levels = c('HSC','MPP','CLP','CMP','Ery_prog','Pre_B','T','B','Neutrophil','Macrophage','Basophil','Ery','Unknown'))
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
