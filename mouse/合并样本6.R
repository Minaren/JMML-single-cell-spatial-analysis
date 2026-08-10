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
####4.3 看一下注释后的marker####
DotPlot(sce.harm,features = unique(genes_to_check),cols = c("grey","#9f2b39")) + theme_bw(base_line_size = 0) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid = element_blank()) + labs(x='',y='')
ggsave("dotplot_anno.png",width = 19,height = 8,dpi = 500)
ggsave("dotplot_anno.pdf",width = 19,height = 8,dpi = 500)

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
install.packages('tidydr')
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
#allmarkers<-read.csv("allmarkers.csv",header=F)
top10 <- allmarkers %>% group_by(cluster) %>% top_n(10,wt = avg_log2FC)
write.csv(top10,'top10_celltype_mouse.csv')
top5 <- allmarkers %>% group_by(cluster) %>% top_n(5,wt = avg_log2FC)
write.csv(top5,'top5_celltype_mouse.csv')

# 设置细胞类型的显示顺序，将HSC放在最前面
cell_types <- c('HSC', 'MPP', 'GMP', 'Erythroblast', 'Ery', 
                'Granolucyte', 'Macrophage', 'Monocyte', 'MK', 
                'B', 'T', 'cDC')
# 确保 `sce.harm` 中的细胞类型顺序与自定义顺序一致
sce.harm$celltype <- factor(sce.harm$celltype, levels = cell_types)
# 定义颜色
mycol <- rev(c("#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", 
               "#92b8da", "#b5aa82", "#de9d3d", "#347852", 
               "#ca8399", "#296097", "#564b84"))
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
  "Ctsg", "Elane", "Myb", "Gfi1", "Mpo", "Ms4a3", #GMP
  "Cd34", "Klf4", "Mef2c", "Pou2f2", "Zeb2",#MPP
  "Cd34", "Mecom", "Hoxa9", "Hlf", "Slamf1"#HSC
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
ggsave("markers.heatmap_celltype_mouse_final.pdf", plot = p, width = 30, height = 30)

####5.2使用文献中的genesets验证####


####细胞比例####
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




####HSC亚群再分析####
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
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)
FeaturePlot()

table(HSC$orig.ident)#查看各组细胞数
prop.table(table(Idents(HSC)))
table(Idents(HSC), HSC$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(HSC), HSC$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

library(ggplot2)
colourCount = length(unique(Cellratio$Var1))
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage_HSC.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_HSC.pdf',width = 7,height = 6,dpi = 500)

celltype = data.frame(ClusterID = 0:2,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0),2] = "HSC_Cd69negative"
celltype[celltype$ClusterID %in% c(2),2] = "HSC_Cd69high"
celltype[celltype$ClusterID %in% c(2),2] = "HSC_Cd69low"
head(celltype)
celltype 
table(celltype$celltype)
HSC$celltype = "NA"
for(i in 1:nrow(celltype)){
  HSC@meta.data[which(HSC@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(HSC@meta.data$celltype)
Idents(HSC) = HSC$celltype
Idents(HSC) = factor(Idents(HSC),levels = c('Bright_HSC','Dark_HSC'))
HSC$celltype <- Idents(HSC)
Idents(HSC) <- factor(Idents(HSC),levels = rev(levels(Idents(HSC))))

DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Cd69"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69"),pt.size=0)

VlnPlot(HSC, features = "Cd69",pt.size=0)+ 
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()


VlnPlot(HSC, features = c("Cd69"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Cd69"),pt.size=0)
VlnPlot(HSC, features = c("Cd52"),pt.size=0)
FeaturePlot(HSC,features = c("Cd34"),min.cutoff =0, max.cutoff =10)
FeaturePlot(HSC,features = c("CD74"),min.cutoff =0, max.cutoff =10)
FeaturePlot(sce.harm,features = c("Cd69"),min.cutoff =0, max.cutoff = 10)
FeaturePlot(HSC, features = "Cd69", 
            min.cutoff = "q10", max.cutoff = "q90", 
            cols = c("lightblue", "darkblue", "red"))
VlnPlot(HSC, features = c("Cd69","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                          "Cd34","Cd63",'Cd9',"Cd117","Cd48","Cd74"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Vim","Sox4","Lgals1","Ddah2","Apoe","Jun","Cmtm7","Dusp1","Cpa3",
"H3f3b","Junb",'Gata2',"Fam117a","Plac8","Cd69"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Mp1","Tek","Gfi1b","Egr1"," Tal1","Gata2","Erg","Pbx1","Meis1",
                          " Hox9","Gata1",'Spi1',"Runx3","Kif1","Kras"),pt.size=0)#,split.by = "orig.ident"
VlnPlot(HSC, features = c("Egr1","Gata2","Junb","Pf4","Car1","Car2"),pt.size=0)
VlnPlot(HSC, features = c("Egr1"),pt.size=0)+
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
ggsave("Cd69.pdf",width = 20,height = 20,units = "cm")
VlnPlot(HSC, features = c("Gata2"),pt.size=0)+
  #width控制箱体宽度，col控制边框颜色，fill控制填充颜色  
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
ggsave("Cd69_mouse_共享特征_gata2.pdf",width = 20,height = 20,units = "cm")
VlnPlot(HSC, features = c("Sox4"),pt.size=0,split.by = "orig.ident")
ggsave("Cd69_mouse_共享特征_Sox4.pdf",width = 40,height = 20,units = "cm")

VlnPlot(HSC, features = c("Egr1","Gata2","Junb"),pt.size=0)

ggsave("Cd69_mouse_共享特征_split.pdf",width = 40,height = 20,units = "cm")
VlnPlot(HSC, features = c("Sox4","Kras","Vim"),pt.size=0,split.by = "orig.ident")
ggsave("HSC_mouse_共享特征_split.pdf",width = 30,height = 20,units = "cm")

sce.harm@meta.data$seurat_clusters
mycolor <- c('lightgrey', "red")
FeaturePlot(sce.harm,features = c("Ly6a","Kit","Slamf1","Cd48","Fgfr1","Cd34","Itgam1")
            ,min.cutoff =0, max.cutoff = 10,cols = mycolor,ncol = 2)

RidgePlot(sce.harm,features ="Kras")
FeatureScatter(HSC, feature1 = "Cd69", feature2 = "Gem")
FeatureScatter(MPP, feature1 = "Cd69", feature2 = "Kras")

FeaturePlot(HSC,features = c("Slamf1"),min.cutoff =0, max.cutoff = 10)
ggsave("HSC_标记.pdf",width = 20,height = 30,units = "cm")


VlnPlot(sce.harm, features = c("Kras","Cd69"),pt.size=0,split.by = "orig.ident")

VlnPlot(sce.harm, features = c("Ccnd1"),pt.size=0,split.by = "orig.ident")
VlnPlot(HSC, features = c("Ctsg","Klf2","Tnf","Fos","Prtn3","Gem","Egr1","Sox4"),pt.size=0,split.by = "orig.ident",ncol = 2)
VlnPlot(HSC, features = c("Snhg3","Egr1"),pt.size=0,split.by = "orig.ident")
ggsave("HSC_标记_差异.pdf",width = 20,height = 40,units = "cm")


HSC <- subset(sce.harm, celltype=="HSC")
MPP<-subset(sce.harm, celltype=="MPP")


VlnPlot(HSC, features = c("Kras","Cd74"),pt.size=0,split.by = "orig.ident")
