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
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(harmony))


####一、多个单细胞样本的合并####
####1.读取并合并数据==####
####1.1读取数据####
###批量读取数据
###设置数据路径与样本名称
dir = c('D:/aJMML/JMML/Kras', 
        'D:/aJMML/JMML/WT')
samples_name= c('Kras', 'WT')

####1.2 批量创建seurat对象####
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #不设置min.cells过滤基因会导致CellCycleScoring报错
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-") 
  }
}
####1.3 给列表命名并保存数据####
setwd("./Integrate_n")
names(scRNAlist) <- samples_name
#system.time(save(scRNAlist, file = "Integrate/scRNAlist0.Rdata")) 
system.time(saveRDS(scRNAlist, file = "scRNAlist0.rds"))

####1.4 使用merge函数将scRNAlist合成一个Seurat对象####
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA 
table(scRNA$orig.ident) 
save(scRNA,file = 'scRNA_orig.Rdata')
#scRNAlist <- SplitObject(scRNA, split.by = "orig.ident") #分割Seurat对象


####二、数据质控####
####2.1 计算各类基因占比####
# 设置可能用到的主题 
theme.set2 = theme(axis.title.x=element_blank()) 
# 设置绘图元素 
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt") 
group = "orig.ident" 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=2) 
dir.create("QC") 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

####2.2 设置质控标准####
minGene=500
maxGene=6000#改成6000
pctMT=10
#maxUMI=15000
#pctMT=10 
#pctHB=1 
####2.3 过滤后可视化####
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                  percent.mt < pctMT) 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0, 
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8)
save(scRNA,file = "scRNA_after_qc.RData")

####三. 数据整合####
####1. seurat锚点整合####
#参考：https://satijalab.org/seurat/articles/integration_large_datasets.html
#Seurat2的整合主要用的是CCA（canonical correlation analysis，典型关联分析）的方法，Seurat3和Seurat4用的是CCA+MNN（mutual nearest neighbor，互近邻）
#锚点整合操作速度很慢，且常常会过度整合，因此在实际操作中，跨物种整合的时候或不同的数据类型如ATAC、蛋白组的数据和单细胞的数据整合的时候，可以使用锚点整合。单纯单细胞数据的整合，可以使用Harmony。
####1.1 读入数据，拆分样本####
rm(list=ls())
library(future) #Seurat并行计算的一个包，不加载这个包不能进行并行计算
options(future.globals.maxSize = 50 * 1024^3) #将全局变量上限调至50G（锚点整合很占内存）

##重新创建没有处理的经过降维等处理的数据 
scRNA <- readRDS("scRNA.rds")
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt", "percent.rb", "percent.HB"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
#做锚点整合需要把样本处理成单个的Seurat对象来两两组合
scRNAlist <- SplitObject(scRNA, split.by = "orig.ident")
#也可以按别的指标（metadata中的）来进行拆分，比如可以按不同的分组来拆分样本，再进行整合。

####1.2 标准化####
#（由于锚点整合会把单个样本两两组合，所以需要单独做标准化）
#SCTransform标准化(使用log标准化还是SCT标准化差别不大)
#如果用log标准化，后面FindIntegrationAnchors和IntegrateData函数的normalization.method参数选'LogNormalize'
scRNAlist <- parallel::mclapply(scRNAlist, FUN=function(x) SCTransform(x), mc.cores = 10) #10个对象最好写10个核，没有10个核少写几个也可以。top命令可以查看服务器有几个核，mc.core设置为1就每次处理一个对象。
# mclapply是lapply的多核版本

####1.3 选择用于整合的高变基因（三步）####
### FindAnchors 
### 每个样本的高变基因不完全一样，SelectIntegrationFeatures可以整合这些高变基因，选出3000个
scRNA.features <- SelectIntegrationFeatures(scRNAlist, nfeatures = 3000) 
### 将每个样本的高变基因都调整成上一步选出的3000个
scRNAlist <- PrepSCTIntegration(scRNAlist, anchor.features = scRNA.features) 
##寻找锚点，运行速度非常慢，至少需要1-2小时
plan("multisession", workers = 10)
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist,
                                        normalization.method = "SCT",  #如果前面是log标准化，这里改成LogNormalize
                                        anchor.features = scRNA.features)

####1.4 锚点整合####
### Integrate 运行速度慢
scRNA.sct.int <- IntegrateData(scRNA.anchors, normalization.method="SCT") #速度慢
plan("sequential") #把并行计算改为单核计算

####1.5 降维，可视化####
### redunction
scRNA <- RunPCA(scRNA.sct.int, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims=50)
pc.num=1:20
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)

### Visual
p <- DimPlot(scRNA, group.by = "orig.ident")
ggsave("UMAP_Samples_integr.pdf", p, width = 8, height = 6)
p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "orig.ident", ncol = 4)
ggsave("UMAP_Samples_Split_integr.pdf", p, width = 18, height = 12)



####2. harmony整合####
#Harmony整合的官网教程及其原理此前已经介绍过：https://www.jianshu.com/p/7c43dc99c4b1
####2.1 准备数据####
rm(list = ls())
load("D:/aJMML/JMML/scRNA/Integrate/Data/scRNA_after_qc.RData")
cellinfo <- subset(scRNA@meta.data, select = c("orig.ident", "percent.mt"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)

####2.2 数据标准化（和锚点整合不同，不需拆分样本，直接标准化）####
### SCT标准化数据:scRNA <- SCTransform(scRNA)
scRNA= NormalizeData(scRNA)
scRNA = FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000, verbose = FALSE)
scRNA = ScaleData(scRNA)

scRNA

####2.3 使用harmony整合数据####
### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
ElbowPlot(scRNA, ndims = 50)
ggsave('elbowplot_n.pdf',width = 6.5,height = 5,dpi = 500)
ggsave('elbowplot_n.png',width = 6.5,height = 5,dpi = 500)
pc.num=1:30
### 整合
#scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20) 
scRNA<-RunHarmony(scRNA,group.by.vars = "orig.ident",project.dim = F,plot_convergence = T)

# group.by.vars参数是设置按哪个分组来整合
# max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
#⚠️RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）

##save seurat object
#saveRDS(scRNA, "scRNA_SCT_harmony.rds") 
saveRDS(scRNA, "scRNA_N_harmony.rds") 

####2.4 分群####
scRNA <- scRNA
scRNA@active.assay
scRNA <- FindNeighbors(scRNA, dims = pc.num,reduction = "harmony") 
scRNA = FindClusters(scRNA, graph.name = "RNA_snn", resolution = 0.5, algorithm = 1)
table(scRNA@meta.data$seurat_clusters)


####2.5 降维及可视化####
scRNA <-RunUMAP(scRNA,reduction="harmony", dims=pc.num)
save(scRNA,file = "scRNA_tu_n.RData")

DimPlot(scRNA,label = T,repel = T) + NoLegend()
ggsave('umap_n.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap_n.png',width = 6,height = 6,dpi = 500)


DimPlot(scRNA,label = F,group.by = 'orig.ident')
ggsave('umap_indi_n.pdf',width = 7,height = 6,dpi = 500)
ggsave('umap_indi_n.png',width = 7,height = 6,dpi = 500)



####四.细胞注释####
####4.1 基因检查####
rm(list = ls())
load("D:/aJMML/JMML/scRNA/Integrate/scRNA_tu_n.RData")
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
                  "Slamf1","Kit","Cd34","Scal1","Cd48","Procr11","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf")
DefaultAssay(scRNA) = "RNA"
DotPlot(scRNA,features = unique(genes_to_check)) + coord_flip()
ggsave("har_gene_show_n.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show_n.png",width = 14,height = 12,dpi = 500)

allmarkers <- FindAllMarkers(scRNA,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_n.csv')

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10.csv')

p<-DoHeatmap(scRNA, features = top10$gene) + NoLegend()
ggsave("markers.heatmap.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap.pdf", plot = p, width = 17, height = 17)

#点图可视化
p <- DotPlot(scRNA, features = unique(top10$gene) ,
             assay='RNA' )  + coord_flip()
plotd<-p+ theme(axis.text.x = element_text(angle = 45, 
                                           vjust = 0.5, hjust=0.5))
#保存可视化图片
ggsave("markers.dotplot.png", plot = plotd, width = 17, height = 49)
ggsave("markers.dotplot.pdf", plot = plotd, width = 17, height = 49)



####4.2 注释####
sce <- scRNA
celltype = data.frame(ClusterID = 0:32,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,1,4,10,20,31),2] = "Pod"
celltype[celltype$ClusterID %in% c(2),2] = "SSBpod"
celltype[celltype$ClusterID %in% c(3),2] = 'NPCa'
celltype[celltype$ClusterID %in% c(5,22,25,32),2] = 'ICb'
celltype[celltype$ClusterID %in% c(6,28),2] = 'SSBm/d'
celltype[celltype$ClusterID %in% c(7),2] = 'ICa' 
celltype[celltype$ClusterID %in% c(8),2] = "CnT" 
celltype[celltype$ClusterID %in% c(9,24),2] = 'ErPrT'
celltype[celltype$ClusterID %in% c(18),2] = 'SSBpr'
celltype[celltype$ClusterID %in% c(11,15),2] = "NPCd" 
celltype[celltype$ClusterID %in% c(12),2] = "NPCc"
celltype[celltype$ClusterID %in% c(13),2] = 'Unknown' 
celltype[celltype$ClusterID %in% c(14),2] = "RVCSBa" 
celltype[celltype$ClusterID %in% c(16),2] = 'NPCb'
celltype[celltype$ClusterID %in% c(17),2] = "End" 
celltype[celltype$ClusterID %in% c(19),2] = "PTA"
celltype[celltype$ClusterID %in% c(23),2] = "DTLH"
celltype[celltype$ClusterID %in% c(26),2] = 'Leu' 
celltype[celltype$ClusterID %in% c(21),2] = 'IPC' 
celltype[celltype$ClusterID %in% c(27),2] = "RVCSBb" 
celltype[celltype$ClusterID %in% c(29),2] = 'Mes'
celltype[celltype$ClusterID %in% c(30),2] = "UBCD" 

head(celltype)
celltype 
table(celltype$celltype)
sce.har$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.har@meta.data[which(sce.har@meta.data$RNA_snn_res.1.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.har@meta.data$celltype)

Idents(sce.har) = sce.har$celltype
Idents(sce.har) = factor(Idents(sce.har),levels = c('NPCa','NPCb','NPCc','NPCd','PTA','RVCSBa','RVCSBb','SSBm/d','SSBpr',
                                                    'SSBpod','CnT','DTLH','ErPrT','Pod','UBCD','IPC','ICa','ICb','Mes','End','Leu','Unknown'))
sce.har$celltype <- Idents(sce.har)
Idents(sce.har) <- factor(Idents(sce.har),levels = rev(levels(Idents(sce.har))))

####4.3 看一下注释后的marker####
DotPlot(sce.har,features = unique(genes_to_check),cols = c("grey","red")) + theme_bw(base_line_size = 0) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid = element_blank()) + labs(x='',y='')
ggsave("dotplot_anno.png",width = 19,height = 8,dpi = 500)
ggsave("dotplot_anno.pdf",width = 19,height = 8,dpi = 500)


library(RColorBrewer)
color_ct=c(brewer.pal(12, "Set3"),"#b3b3b3",
           brewer.pal(5, "Set1"),
           brewer.pal(3, "Dark2"),
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
if(!require(ggunchull))devtools::install_github("sajuukLyu/ggunchull", type = "source")
if(!require(scRNAtoolVis))devtools::install_github('junjunlab/scRNAtoolVis')

####4.4 注释后可视化####
## umap/tsne
clusterCornerAxes(object = sce.har,reduction = 'umap',clusterCol = 'celltype',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(color_ct,0.65)) + NoLegend() +
  scale_fill_manual(values = alpha(color_ct,0.65))
ggsave("umap.pdf",width = 10,height = 9,dpi = 500)
ggsave("umap.png",width = 10,height = 9,dpi = 500)


clusterCornerAxes(object = sce.har,reduction = 'tsne',clusterCol = 'celltype',cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(color_ct,0.65)) +  NoLegend() +
  scale_fill_manual(values = alpha(color_ct,0.65))
ggsave("tsne.pdf",width = 10,height = 9,dpi = 500)
ggsave("tsne.png",width = 10,height = 9,dpi = 500)


clusterCornerAxes(object = sce.har,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("umap_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("umap_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sce.har,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("tsne_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("tsne_indi.png",width = 5.5,height = 5,dpi = 500)


clusterCornerAxes(object = sce.har,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("tsne_indi_s.pdf",width = 15,height = 4,dpi = 500)
ggsave("tsne_indi_s.png",width = 15,height = 4,dpi = 500)


clusterCornerAxes(object = sce.har,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("umap_indi_s.pdf",width = 15,height = 4,dpi = 500)
ggsave("umap_indi_s.png",width = 15,height = 4,dpi = 500)


save(sce.har,file = 'sce_har_anno.RData')

####4.5 注释后查看细胞数####
table(sce.har$celltype)
as.data.frame(table(sce.har$orig.ident,sce.har$celltype))










#以下为备用









####细胞聚类####
scRNA1 <- FindNeighbors(scRNA1, dims = pc.num) 
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
table(scRNA1@meta.data$seurat_clusters)
metadata <- scRNA1@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

####非线性降维####
#tSNE
scRNA1 = RunTSNE(scRNA1, dims = pc.num)
embed_tsne <- Embeddings(scRNA1, 'tsne')   #提取tsne图坐标
write.csv(embed_tsne,'cluster1/embed_tsne.csv')
#group_by_cluster
plot1 = DimPlot(scRNA1, reduction = "tsne", label=T) 
ggsave("cluster1/tSNE.png", plot = plot1, width = 8, height = 7)
#group_by_sample
plot2 = DimPlot(scRNA1, reduction = "tsne", group.by='orig.ident') 
ggsave("cluster1/tSNE_sample.png", plot = plot2, width = 8, height = 7)
#combinate
plotc <- plot1+plot2
ggsave("cluster1/tSNE_cluster_sample.png", plot = plotc, width = 10, height = 5)

#UMAP
scRNA1 <- RunUMAP(scRNA1, dims = pc.num)
embed_umap <- Embeddings(scRNA1, 'umap')   #提取umap图坐标
write.csv(embed_umap,'cluster1/embed_umap.csv') 
#group_by_cluster
plot3 = DimPlot(scRNA1, reduction = "umap", label=T) 
ggsave("cluster1/UMAP.png", plot = plot3, width = 8, height = 7)
#group_by_sample
plot4 = DimPlot(scRNA1, reduction = "umap", group.by='orig.ident')
ggsave("cluster1/UMAP.png", plot = plot4, width = 8, height = 7)
#combinate
plotc <- plot3+plot4
ggsave("cluster1/UMAP_cluster_sample.png", plot = plotc, width = 10, height = 5)

#合并tSNE与UMAP
plotc <- plot2+plot4+ plot_layout(guides = 'collect')
ggsave("cluster1/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)

#scRNAlist是之前代码运行保存好的seurat对象列表，保存了2个样本的独立数据
#数据整合之前要对每个样本的seurat对象进行数据标准化和选择高变基因

##scRNA2对象的降维聚类参考scRNA1的代码

#scRNAlist是之前代码运行保存好的seurat对象列表，保存了2个样本的独立数据
#数据整合之前要对每个样本的seurat对象进行数据标准化和选择高变基因
for (i in 1:length(scRNAlist)) {
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst")
}
##以VariableFeatures为基础寻找锚点，运行时间较长
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist)
##利用锚点整合数据，运行时间较长
scRNA3 <- IntegrateData(anchorset = scRNA.anchors)
dim(scRNA3)
saveRDS(scRNA3, file = "scRNA3.Rds")
scRNA3<-readRDS(file = "scRNA3.Rds")

####==数据质控==####
scRNA <- scRNA3  #以后的分析使用整合的数据进行
##meta.data添加信息
proj_name <- data.frame(proj_name=rep("demo2",ncol(scRNA)))
rownames(proj_name) <- row.names(scRNA@meta.data)
scRNA <- AddMetaData(scRNA, proj_name)
##切换数据集
DefaultAssay(scRNA) <- "RNA"
##计算线粒体和红细胞基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
#head(scRNA@meta.data)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
##绘制小提琴图
#所有样本一个小提琴图用group.by="proj_name"，每个样本一个小提琴图用group.by="orig.ident"
violin <-VlnPlot(scRNA, group.by = "proj_name",  
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot_before_qc.png", plot = pearplot, width = 12, height = 5)

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=500
maxGene=6000#改成6000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(as.factor(scRNA@meta.data$orig.ident)))
violin <-VlnPlot(scRNA, group.by = "proj_name",
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 12, height = 6)

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst")
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA))
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(scRNA, ndims=30, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster1/pca.png", plot = plotc, width = 8, height = 4)
print(c("请选择哪些pc轴用于后续分析？示例如下：","pc.num=1:15"))
#选取主成分
pc.num=1:30
##细胞聚类
scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution = 0.5)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster1/cell_cluster.csv',row.names = F)

####细胞注释####
Basophil_gene<-c("Prss34","Cd200r3","Lmo4")
B_cell_gene<-c("Vpreb3","Cd79a","Cd79b","Cd74","Cd19")
Eosinophil_progenitor_cell_gene<-c("Prg2","Itga2b","Plek","Cebpe","Lmo4","Prss34")
Erythroblast_Car1_high_gene<-c("Blvrb","Rhd","Hmbs","Gata1","Car1")
Erythroblast_Hba_a1_high_gene<-c("Hba-a1","Hbb-bt","Hbb-bs","Blvrb","Hba-a2","Hmbs")
Macrophage_C1qc_high_gene<-c("Cd68","Cd74","Mpeg1","Adgre1","Csf1r","C1qc")
Macrophage_Cd74_high_gene<-c("Cd74","Mpeg1","Cd68","Csf1r","Adgre1")
Macrophage_S100a4_high_gene<-c("Csf1r","Cd74","Mpeg1","Cd68","S100a4")
Megakaryocyte_progenitor_cell_gene<-c("Pf4","Ctla2a","Itga2b","Gata2","Plek","Cd9")
Monocyte_progenitor_cell_Ctsg_high_gene<-c("Prtn3","Elane","Mpo","Ctsg")
Monocyte_progenitor_cell_Prtn3_high_gene<-c("Elane","Mpo","Ighg1","Prtn3")
Monocyte_progenitor_gene<-c("Ms4a6c","F13a1","Irf8","Csf1r")
Multipotent_progenitor_Ctla2a_high_gene<-c("Cd34","Flt3","Ctla2a")
Neutrophil_Fcnb_high_gene<-c("Camp","Ngp","Lcn2","Cebpe","Fcnb")
Neutrophil_Ighg1_high_gene<-c("Mpo","Cd177","Ighg1")
Neutrophil_Lcn2_high_gene<-c("Lcn2")
Neutrophil_Ltf_high_gene<-c("Ngp","Lcn2","Cd177","Cebpe","Ltf")
Neutrophil_Ngp_high_gene<-c("Ngp","Lcn2","Cd177","Cebpe")
Neutrophil_Retnlg_high_gene<-c("Mmp8","Lcn2","Ngp","Camp","Retnlg")
Basophil_gene<-c("Lmo4","Runx1")

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  
genes_to_check =list(
  Basophil<-c("Prss34","Cd200r3","Lmo4"),
  B_cell<-c("Vpreb3","Cd79a","Cd79b","Cd74","Cd19"),
  Eosinophil_progenitor_cell<-c("Prg2","Itga2b","Plek","Cebpe","Lmo4","Prss34"),
  Erythroblast_Car1<-c("Blvrb","Rhd","Hmbs","Gata1","Car1"),
  Erythroblast_Hbaa1<-c("Hba-a1","Hbb-bt","Hbb-bs","Blvrb","Hba-a2","Hmbs"),
  Macrophage_C1qc<-c("Cd68","Cd74","Mpeg1","Adgre1","Csf1r","C1qc"),
  Macrophage_Cd74<-c("Cd74","Mpeg1","Cd68","Csf1r","Adgre1"),
  Macrophage_S100a4<-c("Csf1r","Cd74","Mpeg1","Cd68","S100a4"),
  Megakaryocyte_progenitor_cell<-c("Pf4","Ctla2a","Itga2b","Gata2","Plek","Cd9"),
  Monocyte_progenitor_cell_Ctsg<-c("Prtn3","Elane","Mpo","Ctsg"),
  Monocyte_progenitor_cell_Prtn3<-c("Elane","Mpo","Ighg1","Prtn3"),
  Monocyte_progenitor<-c("Ms4a6c","F13a1","Irf8","Csf1r"),
  Multipotent_progenitor_Ctla2a<-c("Cd34","Flt3","Ctla2a"),
  Neutrophil_Fcnb<-c("Camp","Ngp","Lcn2","Cebpe","Fcnb"),
  Neutrophil_Lcn2<-c("Lcn2"),
  Neutrophil_Ltf<-c("Ngp","Lcn2","Cd177","Cebpe","Ltf"),
  Neutrophil_Ngp<-c("Ngp","Lcn2","Cd177","Cebpe"),
  Neutrophil_Retnlg<-c("Mmp8","Lcn2","Ngp","Camp","Retnlg"),
  Basophil<-c("Lmo4","Runx1")
)
p_all_markers=DotPlot(scRNA, 
                      features = genes_to_check,
                      scale = T,assay='RNA' )+ th +
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggsave('check_paper_markers.pdf',
       height = 8,width = 10)


DotPlot(scRNA,features = Basophil_gene)
DotPlot(scRNA,features = B_cell_gene)
DotPlot(scRNA,features = Eosinophil_progenitor_cell_gene)
DotPlot(scRNA,features = Erythroblast_Car1_high_gene)
DotPlot(scRNA,features = Erythroblast_Hba_a1_high_gene)
DotPlot(scRNA,features = Macrophage_C1qc_high_gene)
DotPlot(scRNA,features = Macrophage_Cd74_high_gene)
DotPlot(scRNA,features = Macrophage_S100a4_high_gene)
DotPlot(scRNA,features = Megakaryocyte_progenitor_cell_gene)
DotPlot(scRNA,features = Monocyte_progenitor_cell_Ctsg_high_gene)
DotPlot(scRNA,features = Monocyte_progenitor_cell_Prtn3_high_gene)
DotPlot(scRNA,features = Monocyte_progenitor_gene)
DotPlot(scRNA,features = Multipotent_progenitor_Ctla2a_high_gene)
DotPlot(scRNA,features = Neutrophil_Fcnb_high_gene)
DotPlot(scRNA,features = Neutrophil_Ighg1_high_gene)
DotPlot(scRNA,features = Neutrophil_Lcn2_high_gene)
DotPlot(scRNA,features = Neutrophil_Ltf_high_gene)
DotPlot(scRNA,features = Neutrophil_Ngp_high_gene)
DotPlot(scRNA,features = Neutrophil_Retnlg_high_gene)
DotPlot(scRNA,features = Basophil_gene)

#findallmarkers
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(scRNA.markers,"scRNA.markers.csv")
write.table (scRNA.markers, file ="scRNA.markers.txt", sep =",", row.names =TRUE, col.names =TRUE, quote =TRUE)
#（1）方法1：找cluster5的marker
cluster0.markers<-FindMarkers(scRNA,ident.1 = 5,ident.2 = c(0,3),min.pct = 0.25)
#（2）方法3：找某个cluster与其他cluster的差异基因（marker）
cluster5.markers <- FindMarkers(seurat.obj, 
                                ident.1 = 5, # 找cluster5的marker
                                ident.2 = c(0, 3), # 和cluster 0-3对比
                                min.pct = 0.25)
head(cluster5.markers, n = 5)
library(dplyr)
#scRNA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#分别用热图、点图、小提琴图可视化
ploth<-DoHeatmap(scRNA, features = top10$gene) + NoLegend()
plotv<-VlnPlot(scRNA, features = top10$gene[1:20],pt.size=0)
p <- DotPlot(scRNA, features = unique(top10$gene) ,
             assay='RNA' )  + coord_flip()
plotd<-p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
#保存可视化图片
ggsave("CellType/markers.hetmap.png", plot = ploth, width = 17, height = 17)
ggsave("CellType/markers.heatmap.pdf", plot = ploth, width = 17, height = 17)
ggsave("CellType/markers.png", plot = plotv, width = 17, height = 15)
ggsave("CellType/markers.pdf", plot = plotv, width = 17, height = 15)
ggsave("CellType/markers.dotplot.png", plot = plotd, width = 17, height = 49)
ggsave("CellType/markers.dotplot.pdf", plot = plotd, width = 17, height = 49)

VlnPlot(CB35,features = c("NKG7","PF4"),slot = "counts",log = TRUE)#小提琴图可视化cluster marker
VlnPlot(CB35,features = c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A","LY2"))
FeaturePlot(CB35,features = c("MS4A1","GNLY","CD3E","CD14","FCER1A"))#降维图可视化cluster marker


##==鉴定细胞类型==##
library(BiocManager)
BiocManager::install("SingleR")
library(SingleR)

dir.create("CellType")
#直接load下载好的数据库
refdata <-get(load("D:/aJMML/JMML/SingleR_ref/ref_Mouse_imm.RData"))
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
#使用参考数据库鉴定
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_Monaco.csv",row.names = F)
scRNA@meta.data$celltype_Monaco = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_Monaco'] <- celltype$celltype[i]}
p1 = DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='tsne')
p2 = DimPlot(scRNA, group.by="celltype_Monaco", repel=T, label=T, label.size=5, reduction='umap')
p3 = p1+p2+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_Monaco.png", p1, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_Monaco.png", p2, width=7 ,height=6)
ggsave("CellType/celltype_Monaco.png", p3, width=10 ,height=5)
#使用DICE参考数据库鉴定
refdata <- DatabaseImmuneCellExpressionData()
# load('~/database/SingleR_ref/ref_DICE_1561s.RData')
# refdata <- ref_DICE
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"CellType/celltype_DICE.csv",row.names = F)
scRNA@meta.data$celltype_DICE = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_DICE'] <- celltype$celltype[i]}
p4 = DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='tsne')
p5 = DimPlot(scRNA, group.by="celltype_DICE", repel=T, label=T, label.size=5, reduction='umap')
p6 = p3+p4+ plot_layout(guides = 'collect')
ggsave("CellType/tSNE_celltype_DICE.png", p4, width=7 ,height=6)
ggsave("CellType/UMAP_celltype_DICE.png", p5, width=7 ,height=6)
ggsave("CellType/celltype_DICE.png", p6, width=10 ,height=5)
#对比两种数据库鉴定的结果
p8 = p1+p4
ggsave("CellType/Monaco_DICE.png", p8, width=12 ,height=5)

##保存数据
saveRDS(scRNA,'scRNA.rds')