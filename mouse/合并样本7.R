####准备工作####
rm(list = ls())
options(stringsAsFactors = F)
set.seed(220625)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))

####一、多个单细胞样本的合并####
####1.读取并合并数据==####
####1.1读取数据####
dir = c('D:/aJMML/JMML/Kras', 
        'D:/aJMML/JMML/WT')
samples_name= c('Kras', 'WT')

####1.2 批量创建seurat对象####
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i],
                                       min.cells=3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-") 
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
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt") 
group = "orig.ident" 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 4, height =9)

####2.2 设置质控标准####
minGene=500
maxGene=6000#改成6000
pctMT=10
####2.3 过滤后可视化####
sc_filt<- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                   percent.mt < pctMT) 
plots = list() 
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(sc_filt, group.by=group, pt.size = 0, 
                       features = plot.featrures[i]) + theme.set2 + NoLegend()} 
violin <- wrap_plots(plots = plots, nrow=2) 
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 4, height =9)
save(sc_filt, file = "sc_filt.RData")

####三. 数据整合####
####2.1 准备数据####
####2.2 数据标准化（和锚点整合不同，不需拆分样本，直接标准化）####
load("D:/aJMML/scRNA1/sc_filt.RData")
sc_n= NormalizeData(sc_filt)
sc_n= FindVariableFeatures(sc_n,selection.method = "vst",nfeatures = 2000, verbose = FALSE)
sc_n = ScaleData(sc_n)
sc_n

####2.3 使用harmony整合数据####
### PCA
sc_n <- RunPCA(sc_n, npcs=50, verbose=FALSE)
ElbowPlot(sc_n, ndims = 50)
ggsave('elbowplot_n.pdf',width = 6.5,height = 5,dpi = 500)
ggsave('elbowplot_n.png',width = 6.5,height = 5,dpi = 500)
sce_ha = sc_n
pc.num=1:30
### 整合
sce_har<-RunHarmony(sce_ha,group.by.vars = "orig.ident",project.dim = F,plot_convergence = T)

####2.4 分群####
sce.har <- sce_har
sce.har@active.assay
sce.har = FindNeighbors(sce.har, dims = pc.num,reduction = "harmony")
sce.har = FindClusters(sce.har, graph.name = "RNA_snn", resolution = 1.2, algorithm = 1)
table(sce.har@active.ident)

####2.5 降维及可视化####
sce.harm = RunUMAP(sce.har, dims = pc.num,reduction = "harmony")
save(sce.harm,file = "sce_harm_tu.RData")
DimPlot(sce.harm,label = T,repel = T) + NoLegend()
ggsave('umap.pdf',width = 6,height = 6,dpi = 500)
ggsave('umap.png',width = 6,height = 6,dpi = 500)

DimPlot(sce.harm,label = F,group.by = 'orig.ident')
ggsave('umap_indi.pdf',width = 7,height = 6,dpi = 500)
ggsave('umap_indi.png',width = 7,height = 6,dpi = 500)


####四.细胞注释####
####4.1 基因检查####
rm(list = ls())
load("D:/aJMML/scRNA1/sce_harm_tu.RData")

allmarkers <- FindAllMarkers(sce.harm,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_celltype_mouse.csv') 

#allmarkers<-read.csv("allmarkers.csv",header=F)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_celltype_mouse.csv')
top10<-read.csv("top10.csv",header=T)

p<-DoHeatmap(sce.harm, features = top10$gene) + NoLegend()
ggsave("markers.heatmap_celltype_mouse.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap_celltype_mouse_top10.pdf", plot = p, width = 30, height = 30)


####4.2 注释####
sce <- sce.harm
celltype = data.frame(ClusterID = 0:22,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,2,5,8,9,13,16),2] = "Granulocyte"
celltype[celltype$ClusterID %in% c(14),2] = "Macrophage"
celltype[celltype$ClusterID %in% c(18),2] = 'Monocyte'
celltype[celltype$ClusterID %in% c(11,22),2] = 'B'
celltype[celltype$ClusterID %in% c(7,17),2] = 'cDC'
celltype[celltype$ClusterID %in% c(21),2] = 'Ery' 
celltype[celltype$ClusterID %in% c(10,12),2] = 'Erythroblast'
celltype[celltype$ClusterID %in% c(15),2] = "HSC"
celltype[celltype$ClusterID %in% c(6),2] = "MPP"
celltype[celltype$ClusterID %in% c(19),2] = "T"
celltype[celltype$ClusterID %in% c(1,3,4),2] = "GMP"
celltype[celltype$ClusterID %in% c(20),2] = "MK"

head(celltype)
celltype 
table(celltype$celltype)
sce.harm$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.harm@meta.data[which(sce.harm@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.harm@meta.data$celltype)

Idents(sce.harm) = sce.harm$celltype
Idents(sce.harm) = factor(Idents(sce.harm),levels = c('HSC','MPP',"GMP",'Erythroblast','Ery',
                                                      'Granulocyte', "Macrophage","Monocyte",'MK',"B","T","cDC"))
sce.harm$celltype <- Idents(sce.harm)
Idents(sce.harm) <- factor(Idents(sce.harm),levels = rev(levels(Idents(sce.harm))))
DimPlot(sce.harm,label = T,repel = T) + NoLegend()
save(sce.harm,file = "sce_har_anno.RData")
ggsave("dotplot_anno.png",width = 19,height = 8,dpi = 500)


####4.4 注释后可视化####
####4.4.1美化方案一####
# 加载必要的包
library(tidyverse)
library(tidydr)
library(magrittr)
library(Seurat)
library(colorfindr)
# 读取数据
load("D:/aJMML/scRNA1/sce_har_anno.RData")
sce <- sce.harm
# 提取 UMAP 坐标并整合到元数据中
meta <- sce@meta.data
umap_coordinates <- Embeddings(sce, "umap")
meta <- cbind(meta, umap_coordinates)
# 颜色定义
col_df <- data.frame(name = unique(meta$celltype)) %>% arrange(name)
mycol <- c(
  "#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", "#92b8da",
  "#b5aa82", "#de9d3d", "#347852", "#ca8399", "#296097", "#564b84"
)
mycol <- setNames(mycol, col_df$name)
# 计算每种细胞类型的中位数坐标，用于标签显示
mid_coord_type <- meta %>%
  dplyr::select(c('celltype', 'UMAP_1', 'UMAP_2')) %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
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

# 保存图像
png(paste0('umap_celltype_with_arrows.png'), width = 8, height = 6, units = 'in', res = 300)
print(pumap2)
dev.off()

pdf(paste0('umap_celltype_with_arrows.pdf'), width = 8, height = 6)
print(pumap2)
dev.off()


####4.4.1美化方案二####
#载入相关R包：
rm(list = ls())
library(tidydr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cols4all)
#读取scRNA-seq测试数据：
##已完成分析的Seurat对象
load("D:/aJMML/scRNA1/sce_har_anno.RData")
sce<-sce.harm
sce[['celltype']] <- Idents(sce)
dim(sce)
celltype <- Idents(sce)
table(celltype)
#UMAP提取方法相同：
UMAP <- as.data.frame(sce@reductions$umap@cell.embeddings)
UMAP <- cbind(UMAP,celltype)
head(UMAP)
#建立自定义主题：
mytheme <- theme_void() + #空白主题，便于我们后期添加UMAP箭头
  theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整
#建立映射，添加散点：
p <- ggplot(data = UMAP, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = celltype),
             size = 0.4,
             alpha = 0.8)
p
#添加实线椭圆置信区间：
p1 <- p +
  stat_ellipse(aes(color = celltype),
               level = 0.95, linetype = 1, show.legend = F) +
  mytheme
p1

#添加虚线椭圆置信区间：
p2 <- p +
  stat_ellipse(aes(color = celltype),
               level = 0.95, linetype = 2, show.legend = F) +
  mytheme
p2
#添加填充型置信区间：
p3 <- p +
  stat_ellipse(aes(color = celltype, fill = celltype),
               level = 0.95, linetype = 1, show.legend = F,
               geom = 'polygon', alpha = 0.1) +
  mytheme
p3
#添加UMAP坐标轴箭头：
p4 <- p3 +
  theme_dr(xlength = 0.2, #x轴长度
           ylength = 0.2, #y轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), #箭头大小/长度
                               ends = 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())
p4
#在图中增加亚群标签:
##计算每个亚群散点的中位数，作为图中标签的坐标
label <- UMAP %>%
  group_by(celltype)%>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
head(label)
p5 <- p4 +
  geom_text(data = label,
            aes(x = UMAP_1, y = UMAP_2, label = celltype),
            fontface = "bold", #粗体强调
            color = 'black', size = 4)
p5
p6 <- p5 +
  guides(color = guide_legend(override.aes = list(size = 8))) #放大图例中的散点
p6
#自定义颜色：
mycol <-c(
  "#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da",
  "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84")
p7 <- p6 +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol)
p7

####五、验证注释准确性####
####5.1使用FIndmarkers中的基因####
allmarkers <- FindAllMarkers(sce.harm,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_celltype_mouse.csv') 
allmarkers<-read.csv("allmarkers.csv",header=F)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_celltype_mouse.csv')
top5 <- allmarkers %>% group_by(cluster) %>% top_n(5,wt = avg_log2FC)
write.csv(top5,'top5_celltype_mouse.csv')

# 设置细胞类型的显示顺序，将HSC放在最前面
cell_types <- c('HSC','MPP',"GMP",'Erythroblast','Ery',
                'Granulocyte', "Macrophage","Monocyte",'MK',"B","T","cDC")
# 确保 `sce.harm` 中的细胞类型顺序与自定义顺序一致
sce.harm$celltype <- factor(sce.harm$celltype, levels = cell_types)
# 定义颜色
mycol <- c("#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", 
               "#92b8da", "#b5aa82", "#de9d3d", "#347852", 
               "#ca8399", "#296097", "#564b84")
names(mycol) <- cell_types
#以下基因均从allmarkers_celltype_mouse中挑选展示
genes_to_check <- rev(c(
  "H2-Eb1", "H2-Aa", "H2-DMb1", "H2-DMa", "Cd74", #cDC
  "Trbc1", "Trbc2", "Cd3d", "Cd3g", "Cd247", #T
  "Vpreb3", "Cd79a", "Cd79b", "Cd74", "Cd19",#B
  "Cxcl12", "Gas6", "Kitl", "Pmp22", "Gpm6b",#MK
  "S100a4", "Clec4a1", "Clec4a3", "Fn1", "Mafb", #Monocyte
  "Cd302", "Cd68", "Csf1r", "Fcgr1", "Adgre4",#Macrophage
  "Cd177", "Ltf", "Retnlg", "Ly6g", "Fpr2",#Neutrophil
  "Blvrb", "Hba-a1", "Hbb-bt", "Hbb-bs", "Hba-a2", "Gata1",#Erythroblast&Ery
  "Ctsg", "Elane","Gfi1", "Mpo", "Ms4a3", #GMP
  "Meis1", "Mef2c", "Flt3", "Hlf","Cd34",#MPP
  "Cdk6","Tal1", "Meis1", "Gata2", "Slamf1"#HSC
))

p <- DoHeatmap(
  object = sce.harm,
  features = genes_to_check,
  group.colors = mycol,
  group.by = "celltype"
) +
  scale_fill_gradientn(colors = c("#92b8da", "white", "#9f2b39"))  # 调整基因表达量的颜色梯度
print(p)

ggsave("markers.heatmap_celltype_mouse_final.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap_celltype_mouse_final.pdf", plot = p, width = 5, height = 5)

####5.2使用文献中的genesets验证####
rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

load("D:/aJMML/scRNA1/sce_har_anno.RData")
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

####5.3Featureplot验证####
library(Seurat)
genes_to_plot <- c("Ly6a", "Kit", "Slamf1", "Cd48", "Gata1", "Cd34")

p <- FeaturePlot(
  sce.harm,  # Seurat对象
  features = genes_to_plot,  # 需要展示的基因
  ncol = 3,  # 每行显示两列
  cols = c("lightgrey", "#9f2b39")  # 设置颜色，从浅灰色到深红色（#9f2b39）
)

ggsave("FeaturePlot_Ly6a_Kit_Slamf1_Cd48_Gata1_Cd34.pdf.png", plot = p, width = 15, height = 10)
ggsave("FeaturePlot_Ly6a_Kit_Slamf1_Cd48_Gata1_Cd34.pdf", plot = p, width = 15, height = 10)


####5.4迁移注释####
#代码在scRNA_mouse_R0.6_6文件夹里
rm(list = ls()) 
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#构建参考数据集
sc<-readRDS("D:/aJMML/scRNA/scRNA_mouse_R0.6_6/GSE122465/notlabel.RDS")#GSE124822
sc= UpdateSeuratObject(object = sc)#seurta版本不一致会报错
str(sc)
slotNames(sc)
data_matrix <- read.table("D:/aJMML/scRNA/scRNA_mouse_R0.6_6/GSE122465/metaInfo.txt", header = TRUE, row.names = 1, sep = "\t")
head(data_matrix)
# 提取 CellType 列
celltype_data <- data_matrix[rownames(sc@meta.data), "CellType..based.on.scmap."]
# 将 CellType 添加到 Seurat 对象的 meta.data 中
sc@meta.data$CellType <- celltype_data
head(sc@meta.data)
sc= UpdateSeuratObject(object = sc)
saveRDS(sc, file = "healthy_mouse_GSE122465.rds")
#构建询问数据集#
load("D:/aJMML/scRNA1/sce_harm_tu.RData")
anchors <- FindTransferAnchors(reference = sc, query = sce.harm, dims = 1:30)
refdata <- sc$CellType
predictions <- TransferData(anchorset = anchors, refdata = refdata, dims = 1:30)
sc_q <- AddMetaData(sce.harm, metadata = predictions)
sc_q$celltype <- sc_q$predicted.id
sc_q$celltype <- as.factor(sc_q$celltype)
head(sc_q$celltype)
table(sc_q$celltype)
DimPlot(sc_q, reduction = "umap",group.by = "celltype", label = TRUE)
ggsave('umap_sc_q_GSE122465.png',width = 6,height = 6,dpi = 500)
ggsave('umap_sc_q_GSE122465.pdf',width = 20,height = 20,dpi = 500)



####六、细胞比例####
rm(list = ls())
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
ggsave('percentage_HSC_mouse.pdf',width = 7,height = 6,dpi = 500)

library(gplots)
tab.1=table(sce.harm$orig.ident,sce.harm$celltype) 
balloonplot(tab.1)
ggsave('balloonplot.png',width = 7,height = 6,dpi = 500)


####七、HSC亚群再分析####
####7.1HSC降维聚类####
rm(list = ls())
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.8)
#HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
HSC_markers <- FindAllMarkers(HSC,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(HSC_markers,'HSC_markers_mouse.csv') 

Idents(HSC)<-HSC$orig.ident
HSC_deg <- FindMarkers(HSC,ident.1 = 'Kras',ident.2 = 'WT', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(HSC_deg,'HSC_mouse_deg.csv') 

CD69high<-subset(HSC,idents = "0")
Idents(CD69high)<-CD69high$orig.ident
CD69high_HSC_deg <- FindMarkers(CD69high,ident.1 = 'Kras',ident.2 = 'WT', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
write.csv(CD69high_HSC_deg,'CD69high_HSC_deg.csv')

VlnPlot(CD69high, features = c("Klf2","Egr3","Fos"),pt.size=0,split.by = "orig.ident",group.by = "orig.ident")
ggsave("CD69highHSC转录因子.pdf",width = 40,height = 20,units = "cm")
VlnPlot(CD69high, features = c("Ccnd1","Gem"),pt.size=0,split.by = "orig.ident",group.by = "orig.ident")
ggsave("CD69highHSC增殖相关基因.pdf",width = 40,height = 20,units = "cm")
VlnPlot(CD69high, features = c("Mmp9","S100a6"),pt.size=0,split.by = "orig.ident",group.by = "orig.ident")
ggsave("CD69highHSC黏附相关基因.pdf",width = 40,height = 20,units = "cm")


HSC$seurat_clusters <- HSC@active.ident
HSC_col <- c("#9e2a2f","#4f85b8","#d87e2d")
save(HSC,file = "mouse_HSC.RData")

p <- DimPlot(HSC, label = TRUE, pt.size = 1) + 
  scale_color_manual(values = HSC_col) +  # 使用自定义颜色
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.grid = element_blank(),  # 去掉网格线
    axis.text = element_blank(),   # 去掉坐标轴文字
    axis.title = element_blank(),  # 去掉坐标轴标题
    legend.position = "right"      # 设置图例位置
  ) 

# 保存图片
ggsave('umap_HSC.png', plot = p, width = 6, height = 6, dpi = 500)
ggsave('umap_HSC.pdf', plot = p, width = 20, height = 20, dpi = 500)

####7.2HSC比例####
table(HSC$orig.ident)#查看各组细胞数
prop.table(table(Idents(HSC)))
table(Idents(HSC), HSC$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(HSC), HSC$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

library(ggplot2)
HSC_col <- c("#d87e2d", "#4f85b8", "#9e2a2f")
# 计算颜色数量
colourCount = length(unique(Cellratio$Var1))
#两组位置互换
Cellratio$Var2 <- factor(Cellratio$Var2, levels = rev(levels(Cellratio$Var2)))
# 绘制柱状图并应用配色
ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.2, size = 0.5, colour = '#222222') + 
  scale_fill_manual(values = HSC_col) +  # 使用自定义颜色
  labs(x = 'Sample', y = 'Ratio', fill = "subpopulation") +  # 修改图例标题为 "subpopulation"
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

# 保存图片
ggsave('percentage_HSC.png', width = 7, height = 6, dpi = 500)
ggsave('percentage_HSC.pdf', width = 7, height = 6, dpi = 500)

DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Cd69"),pt.size=0)


####7.3小提琴图展示CD69####
# 加载所需的库
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsignif)

# 配置颜色
HSC_col <- c("#9e2a2f","#4f85b8","#d87e2d")
HSC_col_deep <- c("#6c1a1f", "#2e5b8a","#9f4c1a")

# 提取 HSC 对象的表达数据并添加到元数据中
HSC$Cd69 <- HSC[["RNA"]]@data["Cd69", ]  # 从 RNA assay 中提取 Cd69 表达值并添加到 HSC 元数据

# 绘制小提琴图，并加入箱线图
p <- VlnPlot(HSC, features = "Cd69", pt.size = 0, group.by = "seurat_clusters") + 
  geom_boxplot(width = .2, col = "black", fill = HSC_col_deep) +  # 添加箱线图
  scale_fill_manual(values = HSC_col) +  # 使用HSC_col着色
  NoLegend() +  # 不显示图例
  theme_classic() +  # 使用经典主题
  labs(title = "CD69 Expression across Clusters")  # 添加标题
print(p)

# 计算群体间显著性差异
# 获取HSC数据的群体和CD69表达数据
group_data <- HSC@meta.data %>% 
  dplyr::select(seurat_clusters, Cd69)  # 提取聚类和Cd69表达数据

# 计算显著性差异，假设我们使用Wilcoxon检验，您可以根据数据选择适当的检验方法
# 比如，这里是计算0群，1群，2群的显著性差异
group_data$seurat_clusters <- factor(group_data$seurat_clusters, levels = c(0, 1, 2))

# 用Wilcoxon检验比较CD69在不同群体中的差异
res <- pairwise.wilcox.test(group_data$Cd69, group_data$seurat_clusters, p.adjust.method = "BH")

# 输出检验结果
print(res)

ggsave("VlnPlot_Cd69_with_significance.png", plot = p, width = 6, height = 6, dpi = 300)
ggsave("VlnPlot_Cd69_with_significance.pdf", plot = p, width = 6, height = 6, dpi = 300)

####7.4 3个亚群功能####
rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
load("D:/aJMML/scRNA1/mouse_HSC.RData")
#1.读取目标genneset文件#
genesets <- read.csv("D:/aJMML/scRNA1/function_mouse_HSC.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

#2.提取分组平均表达矩阵#
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(HSC, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

#3.GSVA富集分析#
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva_res_mouse_HSC.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res_mouse_HSC.csv", row.names = F)

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

# 修改横坐标的列名
colnames(gsva.res) <- gsub("g", "", colnames(gsva.res))

# 重新绘制热图
pheatmap(
  gsva.res,
  show_colnames = TRUE,
  scale = "row",  # 对行进行标准化
  cluster_row = FALSE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  color = my_colors,  # 应用自定义颜色
  filename = "gsva_heatmap_mouse_HSC.pdf"  # 保存文件路径
)





####八、信号通路分析####
rm(list = ls())
####安装及加载R包
library(Seurat)
library(tidyverse)
BiocManager::install("tidyverse")
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
library(org.Mm.eg.db)#人类则用library(org.Hs.eg.db)
####读取数据
load("D:/aJMML/scRNA1/mouse_HSC.RData")
tHSC0<-subset(HSC,ident=0)
####筛选差异基因
## 在寻找差异基因之前，把默认的assay切换为RNA。
DefaultAssay(tHSC0) <- 'RNA'
## 定义好你想要在哪一个分群基础上找差异表达基因
head(tHSC0@meta.data)
Idents(tHSC0) <- 'orig.ident' 
## 在不同cluster/或者celltype中找差异表达基因
#注意这里是Kras/WT，ident.1 = "Kras"

markers <- FindMarkers(tHSC0, ident.1 = "Kras", ident.2 = "WT", only.pos = FALSE, logfc.threshold = 0.25)
head(markers)
markers = markers %>% rownames_to_column('gene') %>% filter(p_val_adj < 0.05)
head(markers)
write.csv(markers,'tHSC0_mouse.csv')
####ID 转换
OrgDb = "org.Mm.eg.db" # 根据物种来指定
gene_convert <- bitr(markers$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)
markers = markers%>%inner_join(gene_convert,by=c("gene"="SYMBOL"))

####8.1 GO富集分析####
ont = "BP" # posiible value: BP, CC, MF, all
go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = FALSE)
head(go.results)

go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = TRUE)
head(go.results)
write.csv(go.results@result, file = "tHSC0_GO_results_BP.csv", row.names = TRUE)
dotplot(go.results,showCategory = 20,label_format=10000)
ggsave("GO_dotplot_BP_tHSC0_mouse.png",plot=p1, width = 10, height = 8, dpi = 300)
ggsave("GO_dotplot_BP_tHSC0_mouse.pdf", plot=p1,width = 10, height = 8)

emapplot(pairwise_termsim(go.results))
ggsave("GO_emapplot_BP_HSC_mouse.png", width = 10, height = 8, dpi = 300)
ggsave("GO_emaplot_BP_HSC_mouse.pdf", width = 10, height = 8)

####8.2 KEGG 富集####
organism = "mmu"  # 小鼠的 KEGG 物种缩写
kegg.results <- enrichKEGG(markers$ENTREZID, organism = organism)
write.csv(kegg.results@result, file = "tHSC0_KEGG_results_mmu.csv", row.names = TRUE)

# 2. 设置可读性
kegg.results <- setReadable(kegg.results, OrgDb = OrgDb, keyType = 'ENTREZID')
head(kegg.results)
# 3. 绘制并保存 KEGG 图形
# 获取 KEGG 富集分析结果中的 Description 列
desc <- kegg.results@result$Description
# 使用 gsub 去掉 " - Mus musculus (house mouse)" 部分
desc_cleaned <- gsub(" - Mus musculus \\(house mouse\\)", "", desc)
# 将修改后的 Description 列重新赋值给 kegg.results
kegg.results@result$Description <- desc_cleaned
dotplot(kegg.results, showCategory = 20, label_format = 10000)
ggsave("KEGG_dotplot_BP_tHSC0_mouse.png", width = 10, height = 8, dpi = 300)
ggsave("KEGG_dotplot_BP_tHSC0_mouse.pdf", width = 10, height = 8)

emapplot(pairwise_termsim(kegg.results))
ggsave("KEGG_emaplot_BP_tHSC0_mouse.png", width = 10, height = 8, dpi = 300)
ggsave("KEGG_emaplot_BP_tHSC0_mouse.pdf", width = 10, height = 8)
# 4. 绘制并保存 Pathview 图形
diff_genes_avg_logFC = markers$avg_log2FC
names(diff_genes_avg_logFC) = markers$ENTREZID
pathview(gene.data = diff_genes_avg_logFC, species = organism, pathway.id = kegg.results@result$ID[3])
ggsave("Pathview_pathway_mmu.png", width = 10, height = 8, dpi = 300)
ggsave("Pathview_pathway_mmu.pdf", width = 10, height = 8)

####8.3 GSEA富集分析####
min.pct = 0.01 ## 至少多少比例的细胞表达这个基因，过滤一些只在极少数细胞中有表达的基因
logfc.threshold = 0.01 ## 过滤掉在两组中几乎没有差异的基因
markers.for.gsea <- FindMarkers(tHSC0, ident.1 = "Kras", ident.2 = "WT", min.pct = min.pct, logfc.threshold=logfc.threshold)
write.csv(markers.for.gsea, "markers.for.gsea_mouse.csv",row.names = TRUE)
# GSEA 要求输入的是一个排好序的列表
Markers_genelist <- markers.for.gsea$avg_log2FC
names(Markers_genelist)= rownames(markers.for.gsea)
head(Markers_genelist)
Markers_genelist <- sort(Markers_genelist, decreasing = T)
#m_df_categories <- msigdbr(species = 'Mus musculus')
#unique_categories <- unique(m_df_categories$gs_cat)
#print(unique_categories)
# 导入MSigDB
categories <- msigdbr_collections()
print(categories)
m_df <- msigdbr(species = 'Mus musculus', category = "C2")#人m_df = msigdbr(species = 'Homo sapiens' , category = "C2")
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea_results_C2 <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea_results_C2)
write_csv(gsea_results_C2 %>% data.frame, "tHSC0_gsea_results_C2.csv")
head(gsea_results_C2)

setid_C2<-c("REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
            "BROWNE_INTERFERON_RESPONSIVE_GENES",
            "SANA_RESPONSE_TO_IFNG_UP",
            "BOSCO_INTERFERON_INDUCED_ANTIVIRAL_MODULE",
            "RADAEVA_RESPONSE_TO_IFNA1_UP",
            "TAKEDA_TARGETS_OF_NUP98_HOXA9_FUSION_3D_DN",
            "WANG_NEOPLASTIC_TRANSFORMATION_BY_CCND1_MYC"
)
gseap2 <- gseaplot2(gsea_results_C2,
                    setid_C2,#富集的ID编号                    
                    title = "GSEA",#标题                    
                    color = mycol,#GSEA线条颜色                    
                    base_size = 20,#基础字体大小                    
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度                    
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图                    
                    ES_geom = "line",#enrichment score用线还是用点"dot"                    
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap2, filename = "GSEA_C2_mouse.pdf",width =12,height =12)


m_df <- msigdbr(species = 'Mus musculus', category = "H")#人m_df = msigdbr(species = 'Homo sapiens' , category = "C2")
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea_results_H <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea_results_H)

# 保存富集结果
write_csv(gsea_results_H %>% data.frame, "gsea_results_H.csv")
head(gsea_results_H)
setid_H<-c("HALLMARK_KRAS_SIGNALING_UP")
gseap_H <- gseaplot2(gsea_results_H,
                    setid_H,#富集的ID编号                    
                    title = "GSEA",#标题                    
                    color = mycol,#GSEA线条颜色                    
                    base_size = 20,#基础字体大小                    
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度                    
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图                    
                    ES_geom = "line",#enrichment score用线还是用点"dot"                    
                    pvalue_table = T) #显示pvalue等信息
ggsave(gseap_H, filename = "GSEA_H_mouse.pdf",width =12,height =12)


rm(list=ls())
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)

load("D:/aJMML/scRNA1/mouse_HSC.RData")
#1.读取目标genneset文件#
msgdC2 = msigdbr(species = "Mus musculus", category = "H")
genesets= msgdC2 %>% split(x = .$gene_symbol, f = .$gs_name)
genesets <- msgdC2 %>% 
  split(x = .$gene_symbol, f = gsub("HALLMARK_", "", .$gs_name))

# 查看结果
head(genesets)

#2.提取分组平均表达矩阵#
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(HSC, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

#3.GSVA富集分析#
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "HSC_gsva.res_mouse.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "HSC_gsva_res_mouse.csv", row.names = F)

my_colors <- colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100) 
library(pheatmap)
# 绘制热图并应用自定义颜色
p<-pheatmap(
  gsva.res,
  show_colnames = TRUE,
  scale = "row",  # 对行进行标准化
  cluster_row = TRUE,  # 不聚类行
  cluster_cols = FALSE,  # 不聚类列
  color = my_colors  # 应用自定义颜色
)
ggsave("HSC_mouse_gsea.png", plot = p, width = 17, height = 17)
ggsave("HSC_mouse_gsea.pdf", plot = p, width = 30, height = 30)


