####准备工作####
rm(list = ls())
library(Seurat) 
library(tidyverse) 
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(stringr)

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
ggsave("vlnplot_before_qc3.pdf", plot = violin, width = 9, height = 15)

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
ggsave("vlnplot_after_qc3.pdf", plot = violin, width = 10, height = 8)
save(sc_filt, file = "sc_filt.RData")

####三. 数据整合####
####2.1 准备数据####
####2.2 数据标准化####
load("D:/aJMML/scRNA/sc_filt.RData")
sc_n= NormalizeData(sc_filt)
sc_n= FindVariableFeatures(sc_n,selection.method = "vst",nfeatures = 2000, verbose = FALSE)
sc_n = ScaleData(sc_n)
sc_n
sc_n <- RunPCA(sc_n, npcs=50, verbose=FALSE)
ElbowPlot(sc_n, ndims = 50)
pc.num=1:30
sc_n <- sc_n %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
sc_n <- FindNeighbors(sc_n, dims=pc.num) %>% FindClusters()
DimPlot(sc_n, label = T)
saveRDS(sc_n,"sc_nRDS")
####1.1 读入数据，拆分样本####
rm(list=ls())
library(future) #Seurat并行计算的一个包，不加载这个包不能进行并行计算
options(future.globals.maxSize = 50 * 1024^3) #将全局变量上限调至50G（锚点整合很占内存）

##重新创建没有处理的经过降维等处理的数据 
#scRNA <- readRDS("sc_nRDS")
scRNA <- sc_n
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.MT"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
#做锚点整合需要把样本处理成单个的Seurat对象来两两组合
scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
#也可以按别的指标（metadata中的）来进行拆分，比如可以按不同的分组来拆分样本，再进行整合。
#SCTransform标准化(⚠️使用log标准化还是SCT标准化差别不大)
#如果用log标准化，后面FindIntegrationAnchors和IntegrateData函数的normalization.method参数选'LogNormalize'
scRNAlist <- parallel::mclapply(scRNAlist, FUN=function(x) SCTransform(x)) #10个对象最好写10个核，没有10个核少写几个也可以。top命令可以查看服务器有几个核，mc.core设置为1就每次处理一个对象。
# mclapply是lapply的多核版本

### FindAnchors 
### 每个样本的高变基因不完全一样，SelectIntegrationFeatures可以整合这些高变基因，选出3000个
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000) 
### 将每个样本的高变基因都调整成上一步选出的3000个
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features) 
##寻找锚点，运行速度非常慢，至少需要1-2小时
plan("sequential")
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,
                                        normalization.method = "SCT",  #如果前面是log标准化，这里改成LogNormalize
                                        anchor.features = scRNA.features)
### Integrate 运行速度慢
scRNA.sct.int <- IntegrateData(scRNA.anchors, normalization.method="SCT") #速度慢
#把并行计算改为单核计算

### redunction
scRNA <- RunPCA(scRNA.sct.int, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims=50)
pc.num=1:30
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)

### Visual
p <- DimPlot(scRNA, group.by = "orig.ident")
ggsave("UMAP_Samples_integr.pdf", p, width = 8, height = 6)
p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("UMAP_Samples_Split_integr.pdf", p, width = 18, height = 12)


#scRNA <- readRDS("scRNA_SCT_int.rds")
sc_integr<-scRNA
sc_integr <- FindNeighbors(sc_integr, dims = 1:30) 
sc_integr <- FindClusters(sc_integr, resolution=0.6)
sc_integr <- RunUMAP(sc_integr, reduction = "pca", dims = 1:30)
table(sc_integr@active.ident)
save(sc_integr,file = "sc_seurat_integr.RData")


DimPlot(sc_integr,label = T,repel = T) + NoLegend()
ggsave('umap_sc_integr_cluster.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap_sc_integr_cluster.png',width = 6,height = 6,dpi = 500)


DimPlot(sc_integr,label = F,group.by = 'orig.ident')
ggsave('umap_sc_integr_indi.pdf',width = 7,height = 6,dpi = 500)
ggsave('umap_sc_integr_indi.png',width = 7,height = 6,dpi = 500)


####四.细胞注释####
####4.1 基因检查####
rm(list = ls())

genes_to_check =c("CD2","ITGAM","ITGAL","ANPEP","CD14","FUT4","CD19","CD33","IL6R","IL7R","MME","CD7","CD127",
                  "TFRC","EPOR","ITGA2B","GP1BA","PTPRC","CD34","CD38","THY1","ITGA6","KIT","FLT3","IL3RA")
DefaultAssay(sc_integr) = "RNA"
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_HSC5.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show_HSC5.png",width = 14,height = 12,dpi = 500)

allmarkers <- FindAllMarkers(sc_integr,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_human.csv') 

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_human.csv')

p<-DoHeatmap(sc_integr, features = top10$gene) + NoLegend()
ggsave("markers.heatmap.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap.pdf", plot = p, width = 17, height = 17)

#点图可视化
p <- DotPlot(sc_integr, features = unique(top10$gene) ,
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

####4.2 注释####
rm(list = ls())
load("D:/aJMML/scRNA2/sc_seurat_integr.RData")
sc_integra<- sc_integr
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
sc_integra$celltype = "NA"
for(i in 1:nrow(celltype)){
  sc_integra@meta.data[which(sc_integra@meta.data$integrated_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sc_integra@meta.data$celltype)

Idents(sc_integra) = sc_integra$celltype
Idents(sc_integra) = factor(Idents(sc_integra),levels = c('HSC','MPP',"LMPP",'CLP','CMP',"GMP",'MEP','Ery_prog',"MK",'B','Granulocyte'))
sc_integra$celltype <- Idents(sc_integra)
Idents(sc_integra) <- factor(Idents(sc_integra),levels = rev(levels(Idents(sc_integra))))

DimPlot(sc_integra,label = T,repel = T) + NoLegend()
save(sc_integra,file = "sce_har_anno_human.RData")
ggsave("dotplot_anno.png",width = 19,height = 8,dpi = 500)


library(tidyverse)
library(tidydr)
library(magrittr)
library(Seurat)
library(colorfindr)

# 读取数据
load("D:/aJMML/scRNA2/sce_har_anno_human.RData")
sce <- sc_integra

# 提取 UMAP 坐标并整合到元数据中
meta <- sce@meta.data
umap_coordinates <- Embeddings(sce, "umap")
meta <- cbind(meta, umap_coordinates)

# 设置 celltype 的因子顺序
celltype_order <- c('HSC','MPP',"LMPP",'CLP','CMP',"GMP",'MEP','Ery_prog',"MK",'B','Granulocyte')
meta <- meta %>% filter(celltype %in% celltype_order)
meta$celltype <- factor(meta$celltype, levels = celltype_order)

# 颜色定义并与顺序对应
mycol <- c(
  "#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", "#92b8da",
  "#b5aa82", "#de9d3d", "#347852", "#ca8399", "#296097"
)
mycol <- setNames(mycol, celltype_order)

# 计算每种细胞类型的中位数坐标，用于标签显示
mid_coord_type <- meta %>%
  dplyr::select(celltype, UMAP_1, UMAP_2) %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# 绘制基础 UMAP 图
pumap1 <- ggplot(meta, aes(UMAP_1, UMAP_2)) +
  geom_point(size = 0.2, aes(color = celltype)) +
  geom_text(data = mid_coord_type, size = 5,
            aes(UMAP_1, UMAP_2, label = celltype)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = mycol)

# 添加 UMAP 坐标轴箭头
pumap2 <- pumap1 +
  theme_dr(
    xlength = 0.2, # X 轴箭头长度
    ylength = 0.2, # Y 轴箭头长度
    arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")
  ) +
  theme(panel.grid = element_blank())

# 保存图像（PNG）
png("umap_celltype_with_arrows.png", width = 8, height = 6, units = 'in', res = 300)
print(pumap2)
dev.off()

# 保存图像（PDF）
pdf("umap_celltype_with_arrows.pdf", width = 8, height = 6)
print(pumap2)
dev.off()


rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

load("D:/aJMML/scRNA2/sce_har_anno_human.RData")
#1.读取目标genneset文件#
genesets <- read.csv("D:/aJMML/scRNA1/gsva_mouse_cluster.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

#2.提取分组平均表达矩阵#
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(sce.harm, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

#3.GSVA富集分析#
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res_mouse.csv", row.names = F)

my_colors <- colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100) 
library(pheatmap)
# 绘制热图并应用自定义颜色
pheatmap(
  gsva.res,
  show_colnames = TRUE,
  scale = "row",  # 对行进行标准化
  cluster_row = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  color = my_colors  # 应用自定义颜色
)
























































DefaultAssay(sc_integr) = "RNA"
VlnPlot(sc_integr, 
        features = c("CD69"),
        pt.size = 0,
        )
#split.by = "orig.ident"
ggsave("vln4.pdf",width = 40,height = 20,units = "cm")

HSC <- subset(sc_integr, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.1)
#HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
library(clustree)
clustree(HSC)
Idents(HSC) <- "RNA_snn_res.0.5"
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)


####4.3 看一下注释后的marker####
genes_to_check =c("CD2","ITGAM","ITGAL","ANPEP","CD14","FUT4","CD19","CD33","IL6R","IL7R","MME","CD7","CD127",
                  "TFRC","EPOR","ITGA2B","GP1BA","PTPRC","CD34","CD38","THY1","ITGA6","KIT","FLT3","IL3RA")
DotPlot(sc_integr,features = unique(genes_to_check),cols = c("blue","red")) + theme_bw(base_line_size = 0) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid = element_blank()) + labs(x='',y='')+ coord_flip()
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
clusterCornerAxes(object =sc_integr,reduction = 'umap',clusterCol = 'celltype',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T)
ggsave("umap2.pdf",width = 10,height = 9,dpi = 500)
ggsave("umap2.png",width = 10,height = 9,dpi = 500)

clusterCornerAxes(object =sc_integr,reduction = 'umap',clusterCol = 'seurat_clusters',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T)
ggsave("umap_cluster_number.pdf",width = 10,height = 9,dpi = 500)
ggsave("umap2.png",width = 10,height = 9,dpi = 500)





clusterCornerAxes(object = sc_integr,reduction = 'tsne',clusterCol = 'celltype',cellLabel = T,cellLabelSize = 5,
                  noSplit = T)
ggsave("tsne1.pdf",width = 10,height = 9,dpi = 500)
ggsave("tsne1.png",width = 10,height = 9,dpi = 500)


clusterCornerAxes(object = sc_integr,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("umap_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("umap_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sc_integr,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("tsne_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("tsne_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sc_integr,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("tsne_indi_s.pdf",width = 8,height = 4,dpi = 500)
ggsave("tsne_indi_s.png",width = 8,height = 4,dpi = 500)


clusterCornerAxes(object = sc_integr,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("umap_indi_s.pdf",width = 8,height = 4,dpi = 500)
ggsave("umap_indi_s.png",width = 8,height = 4,dpi = 500)


save(sc_integr,file = 'sc_seurat_integr_anno.RData')

####4.5 注释后查看细胞数####
table(sc_integr$celltype)
as.data.frame(table(sc_integr$orig.ident,sc_integr$celltype))

VlnPlot(HSC, features = c("CD69"),pt.size=0,split.by = "orig.ident")
ggsave("Cd69_huamn_split.pdf",width = 40,height = 20,units = "cm")
VlnPlot(HSC, features = c("CD69","CD36"),pt.size=0)
ggsave("Cd69_huamn.pdf",width = 40,height = 40,units = "cm")
FeaturePlot(HSC,features = c("CD69","CD36"),min.cutoff =1, max.cutoff = 10)

VlnPlot(HSC, features = c("CD69","ERG1","GATA2","JUNB",'MALAT1',"FTL
"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("CD69","ERG1","GATA2","JUNB",'MALAT1',"FTL
"),pt.size=0)
VlnPlot(HSC, features = c("MLLT3","AVP","JUNB",'CRHBP',"CD74"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("MLLT3","AVP","JUNB"),pt.size=0)
ggsave("Cd69_high.pdf",width = 60,height = 20,units = "cm")
ggsave("Cd69_CD36_ID5.pdf",width = 40,height = 20,units = "cm")
