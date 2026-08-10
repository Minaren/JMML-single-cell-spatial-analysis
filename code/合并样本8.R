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
sce.harm<-JoinLayers(sce.harm)
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
#载入相关R包：
rm(list = ls())
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
ggsave("markers.heatmap_celltype_mouse_final.pdf", plot = p, width = 17, height = 17)

####5.2使用文献中的genesets验证####
rm(list=ls())
load("D:/aJMML/scRNA1/sce_har_anno.RData")
library(GSVA)
library(pheatmap)

# 1. 读取目标 geneset 文件
genesets <- read.csv("D:/aJMML/scRNA1/gsva_mouse_cluster.csv", header = FALSE)
genesets <- subset(genesets, select = c("V1", "V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)   # 生成 list

# 2. 提取分组平均表达矩阵
expr <- AverageExpression(sce.harm, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ]   # 选取非零基因
expr <- as.matrix(expr)

# 3. GSVA 富集分析（GSVA ≥ 1.50.0 的写法）
# 创建 ssGSEA 专用的参数对象
param <- ssgseaParam(expr, genesets)   # 可以在此添加 minSize, maxSize 等参数
gsva.res <- gsva(param)

# 保存和写出结果
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = FALSE)
write.csv(gsva.df, "gsva_res_mouse.csv", row.names = FALSE)

# 4. 绘制热图
my_colors <- colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100)
pheatmap(
  gsva.res,
  filename = "heatmap_gsva.pdf",  # 或 .png、.jpeg，根据后缀自动识别格式
  show_colnames = TRUE,
  scale = "row",
  cluster_row = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  width = 10,                     # 图片宽度（英寸）
  height = 6                      # 图片高度（英寸）
)


####5.3Featureplot验证####
library(Seurat)
library(patchwork)

genes_to_plot <- c("Ly6a", "Kit", "Slamf1", "Cd48", "Cd34", "Flt3")

# 生成每个基因的 FeaturePlot 并移除坐标轴和图例
plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(sce.harm, features = gene, cols = c("lightgrey", "#9f2b39")) + NoAxes() + NoLegend()
})

# 使用 patchwork 拼接图
p <- wrap_plots(plots, ncol = 2)

# 保存图像
ggsave("FeaturePlot_custom.pdf", plot = p, width = 10, height = 15)


####六、细胞比例（简洁分组柱状图）####
rm(list = ls())
library(Seurat)
library(tidyverse)

load("D:/aJMML/scRNA1/sce_har_anno.RData")

# 细胞类型顺序与颜色
cell_types <- c('HSC','MPP',"GMP",'Erythroblast','Ery',
                'Granulocyte', "Macrophage","Monocyte",'MK',"B","T","cDC")
sce.harm$celltype <- factor(sce.harm$celltype, levels = cell_types)

# 1. 计算每个样本 × 细胞类型的绝对计数
count_df <- sce.harm@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n(), .groups = "drop")

# 2. 设置样本顺序：WT 在前，Kras 在后
count_df$orig.ident <- factor(count_df$orig.ident, levels = c("WT", "Kras"))

# 3. 单张并列柱状图
p <- ggplot(count_df, aes(x = celltype, y = count, fill = orig.ident)) +
  geom_col(position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("WT" = "#4f85b8", "Kras" = "#9f2b39")) +  # WT蓝色，Kras红色
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Cell count", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# 保存
ggsave("cellcount_bar.pdf", p, width = 10, height = 5, dpi = 500)
ggsave("cellcount_bar.png", p, width = 10, height = 5, dpi = 500)


####重点基因表达####
library(Seurat)
library(ggplot2)
library(patchwork)

# 如果之前没有加载 sce.harm，请先加载
load("D:/aJMML/scRNA1/sce_har_anno.RData")

genes_to_plot <- c("Ly6a","Cd34", "Slamf1", "Cd48", "Kit","Flt3")
df <- FetchData(sce.harm, vars = c("UMAP_1", "UMAP_2", genes_to_plot))

# HSC 细胞坐标
hsc_cells <- WhichCells(sce.harm, idents = "HSC")
highlight_df <- df[hsc_cells, ]


plot_gene_clean <- function(gene) {
  colnames(df)[colnames(df) == gene] <- "expr"
  cutoff <- quantile(df$expr, 0.95)
  if (cutoff == 0) cutoff <- 0.01
  
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey90", size = 0.2, alpha = 0.6) +
    geom_point(data = subset(df, expr > cutoff),
               aes(color = expr, size = expr), alpha = 0.7) +
    # 蓝色渐变（从你提供的颜色中选取）
    scale_color_gradientn(colors = c("#FFF5E1", "#de9d3d", "#b2675e")) +
    scale_size_continuous(range = c(0.4, 1.8)) +
    # 红色椭圆标记 HSC
    stat_ellipse(data = highlight_df, aes(UMAP_1, UMAP_2),
                 level = 0.95, geom = "polygon",
                 fill = "#9f2b39", alpha = 0.15,       # 红色填充稍明显
                 color = "#9f2b39", linewidth = 0.8) +
    theme_void() + NoLegend() +
    labs(title = gene) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  colnames(df)[colnames(df) == "expr"] <<- gene
  return(p)
}

# 出图并保存
plots <- lapply(genes_to_plot, plot_gene_clean)
p <- wrap_plots(plots, ncol = 2)
ggsave("FeaturePlot_HSC_ellipse.pdf", p, width = 10, height = 15)



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
cell_types <- c('HSC','MPP',"GMP",'Erythroblast','Ery',
                'Granulocyte', "Macrophage","Monocyte",'MK',"B","T","cDC")

sce.harm$celltype <- factor(sce.harm$celltype, levels = cell_types)

mycol <- c("#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", 
           "#92b8da", "#b5aa82", "#de9d3d", "#347852", 
           "#ca8399", "#296097", "#564b84")
names(mycol) <- cell_types

# 保证 Cellratio$Var1 也是 factor 并有正确顺序
Cellratio$Var1 <- factor(Cellratio$Var1, levels = cell_types)

# 画图 + 自定义颜色
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "identity", width = 0.2, size = 0.5, colour = '#222222') + 
  scale_fill_manual(values = mycol) +  # 👈 加入这一行
  theme_classic() +
  labs(x='Sample', y = 'Ratio') +
  theme(panel.border = element_rect(fill=NA, color="black", size=0.5, linetype="solid"))

# 保存图像
ggsave('percentage.png', width = 7, height = 6, dpi = 500)
ggsave('percentage_HSC_mouse.pdf', width = 7, height = 6, dpi = 500)


library(gplots)
tab.1=table(sce.harm$orig.ident,sce.harm$celltype) 
balloonplot(tab.1)
ggsave('balloonplot.png',width = 7,height = 6,dpi = 500)


# ================================================================
#  HSC 亚群再分析
#  包含：降维聚类、UMAP、比例、功能热图、表面标志物筛选、拟时序
# ================================================================
library(Seurat)
library(tidyverse)
library(patchwork)
library(tidydr)
library(pheatmap)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(slingshot)
library(viridis)
library(mgcv)

# ---- 0. 全局配色与主题 ----
hsc_cluster_cols <- c("0" = "#9f2b39", "1" = "#296097", "2" = "#d6873b")
heat_ann_colors <- list(Cluster = c("Cluster 0" = "#9f2b39",
                                    "Cluster 1" = "#296097",
                                    "Cluster 2" = "#d6873b"))
# 热图渐变：浅蓝-白-红
heat_gradient <- c("#92b8da", "white", "#9f2b39")
# 点图渐变：白-红
dot_gradient <- c("white", "#9f2b39")
# 全局透明度
global_alpha <- 0.8

# 主题：无网格，细线
theme_clean <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(size = 0.4, color = "grey30"),
          axis.ticks = element_line(size = 0.4, color = "grey30"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5))
}


# ================================================================
# 单细胞转录组分析：HSC 亚群功能与拟时序
# 图A：分样本 UMAP
# 图B：功能模块热图
# 图C：表面标志物表达谱
# 图D：Cd69 小提琴图
# 图E：细胞比例柱状图
# 图F：拟时序组合图（密度/CD69/模块/TF）
# ================================================================

# ---- 加载必要的包 ----
library(Seurat)
library(tidyverse)
library(pheatmap)
library(slingshot)
library(mgcv)
library(patchwork)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggpubr)
library(viridis)

# ---- 0. 自定义主题（若环境无 theme_dr / theme_clean 可临时定义） ----
# theme_dr 和 theme_clean 常源于 ggsci 或自定义脚本，此处若缺失可注释掉或自行实现
# 例：theme_dr <- function(...) theme_classic()   # 简化版

# ================================================================
# 1. 数据加载与降维聚类
# ================================================================
rm(list = ls())
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype == "HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.8)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)

Idents(HSC) <- "seurat_clusters"
HSC$seurat_clusters <- factor(HSC$seurat_clusters, levels = c("0", "1", "2"))

# 定义亚群颜色（全局使用）
hsc_cluster_cols <- c("0" = "#9f2b39", "1" = "#296097", "2" = "#d6873b")

# ================================================================
# 2. 图A：分样本 UMAP
# ================================================================
p_WT <- DimPlot(subset(HSC, orig.ident == "WT"), group.by = "seurat_clusters",
                label = TRUE, pt.size = 1.2, cols = hsc_cluster_cols, alpha = 0.7) +
  ggtitle("WT") + 
  theme_dr(xlength = 0.15, ylength = 0.15,
           arrow = grid::arrow(length = unit(0.08, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())

p_Kras <- DimPlot(subset(HSC, orig.ident == "Kras"), group.by = "seurat_clusters",
                  label = TRUE, pt.size = 1.2, cols = hsc_cluster_cols, alpha = 0.7) +
  ggtitle("Kras") + 
  theme_dr(xlength = 0.15, ylength = 0.15,
           arrow = grid::arrow(length = unit(0.08, "inches"), type = "closed")) +
  theme(panel.grid = element_blank())

p_umap <- p_WT | p_Kras
ggsave("Fig_HSC_UMAP.pdf", p_umap, width = 10, height = 5, device = "pdf")

# ================================================================
# 3. 图B：功能模块热图（Hallmark 公认基因集）
# ================================================================
# 基因集定义（保持不变）
inflam_genes <- c("CCL20","CCL4","CCL5","CCL7","CCL8","CD69","CD74","CSF1","CSF2","CSF3",
                  "CXCL1","CXCL10","CXCL2","CXCL3","CXCL8","ICAM1","IL12B","IL15","IL18",
                  "IL1A","IL1B","IL6","IRF7","LTA","PTGS2","SELE","SELP","TNF","VCAM1")
ifna_genes <- c("BST2","EIF2AK2","GBP1","GBP2","IFI27","IFI35","IFI6","IFIT1","IFIT2",
                "IFIT3","IFITM1","IFITM2","IFITM3","IRF1","IRF7","IRF9","ISG15","MX1",
                "MX2","OAS1","OAS2","OAS3","RSAD2","STAT1","STAT2","USP18")
ifng_genes <- c("CASP1","CCL2","CCL5","CCL7","CD274","CIITA","CXCL10","CXCL9","GBP1",
                "GBP2","GBP4","GBP5","ICAM1","IFNG","IFNGR1","IFNGR2","IRF1","IRF7",
                "IRF9","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","RSAD2","SOCS1",
                "STAT1","STAT2","TAP1","TAP2","TNF","USP18","VCAM1")
ifn_all <- unique(c(ifna_genes, ifng_genes))
kras_genes <- c("ADAMTS1","AKT1","APAF1","AREG","ATF3","BCL2","BMP2","BTG2","CCL2",
                "CCL20","CCL4","CCL5","CCND1","CD44","CDKN1A","CFLAR","CLU","CSF1",
                "CSF2","CSF3","CTGF","CXCL1","CXCL10","CXCL2","CXCL3","CXCL8","CYR61",
                "DUSP1","DUSP4","DUSP5","DUSP6","EDN1","EFNB1","EGR1","EGR2","EGR3",
                "EIF4EBP1","EPHA2","EREG","ETV5","FGF2","FOS","FOSB","FOSL1","FOSL2",
                "FST","GADD45A","GADD45B","GAPDH","GEM","HBEGF","HMGA2","HMOX1",
                "ICAM1","IER2","IGF1R","IL13RA1","IL1A","IL1B","IL6","IL7R","IRF1",
                "JUN","JUNB","KDM5B","KLF2","KLF4","KLF6","LIF","MAP2K3","MAPK8",
                "MCL1","MEST","MUC1","MYC","NFKB1","NFKB2","NR4A1","NR4A2","NR4A3",
                "PDGFB","PLAUR","PLK2","PLK3","PMEPA1","PPP1R15A","PTGER4","PTGS2",
                "PTPRE","RAC1","RAP1A","REL","RELB","RHOB","SAT1","SERPINB2","SERPINE1",
                "SKP2","SLC2A1","SLC2A3","SLC4A7","SNHG12","SPRY1","SPRY2","SPRY4",
                "SYK","TBX2","TFPI2","TGFBR2","THBS1","TIMP1","TNF","TRIB3","UPP1",
                "VEGFA","VIM","WT1")
myeloid_genes <- c("AIF1","ASXL1","BCL11A","BMP4","CCR1","CCR2","CD34","CEBPA",
                   "CEBPB","CEBPE","CLEC4D","CSF1","CSF1R","CSF2","CSF2RA","CSF2RB",
                   "CSF3","CSF3R","CTSG","ELANE","FCER2","FLT3","GATA1","GATA2",
                   "GFI1","HOXB4","ID2","IFNG","IKZF1","IL3","IL5","IRF8","ITGAM",
                   "JAG1","KIT","KITLG","LY6G","LYL1","MAFB","MEIS1","MPO","MYB",
                   "NFE2","NOTCH1","PRTN3","PU1","RUNX1","S100A8","S100A9","SPI1",
                   "STAT3","STAT5A","STAT5B","TAL1","TLR2","TNF","ZFP36")
stemness_genes <- c("Procr", "Ly6a", "Mllt3", "Mecom", "Hlf", "Gfi1",
                    "Tek", "Fgd5", "Hoxb5", "Mpl", "Cd34", "Kit",
                    "Hoxa9", "Meis1", "Fli1", "Gata2", "Tcf15", "Cdkn1c")

# 组装模块列表
module_list <- list(
  Stemness              = stemness_genes,
  Myeloid_Priming       = myeloid_genes,
  Inflammatory_Response = inflam_genes,
  IFN_Response          = ifn_all,
  KRAS_Signaling        = kras_genes
)

# 清除旧评分列，计算新评分
old_cols <- grep("_Score$", colnames(HSC@meta.data), value = TRUE)
if (length(old_cols)) HSC@meta.data <- HSC@meta.data[, !colnames(HSC@meta.data) %in% old_cols]

for (nm in names(module_list)) {
  HSC <- AddModuleScore(HSC, features = module_list[nm],
                        name = paste0(nm, "_Score"), ctrl = min(100, nrow(HSC)))
  colnames(HSC@meta.data)[ncol(HSC@meta.data)] <- paste0(nm, "_Score")
}

# 提取平均评分矩阵并绘制热图
score_names <- paste0(names(module_list), "_Score")
avg <- FetchData(HSC, vars = c(score_names, "seurat_clusters")) %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("seurat_clusters")
rownames(avg) <- paste0("Cluster ", rownames(avg))
colnames(avg) <- str_remove(colnames(avg), "_Score$") %>% str_replace_all("_", " ")
mat <- t(avg)

ann_col <- data.frame(Cluster = colnames(mat), row.names = colnames(mat))
ann_colors <- list(Cluster = c("Cluster 0" = "#9f2b39", "Cluster 1" = "#296097", "Cluster 2" = "#d6873b"))

pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
         color = colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100),
         border_color = NA, treeheight_row = 12, treeheight_col = 8,
         annotation_col = ann_col, annotation_colors = ann_colors,
         show_colnames = TRUE, show_rownames = TRUE, fontsize = 11, fontsize_row = 10,
         cellwidth = 45, cellheight = 22, main = "Core Functional Modules (Hallmark)",
         filename = "Fig_HSC_FunctionalHeatmap.pdf")

# ================================================================
# 4. 图C：表面标志物筛选与组合图
# ================================================================
# 获取 cell surface 基因符号 (GO:0009986)
surface_symbols <- unique(AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0009986",
                                                columns = "SYMBOL", keytype = "GOALL")$SYMBOL)
surface_symbols <- str_to_title(surface_symbols)

target_mods <- c("Inflammatory_Response", "KRAS_Signaling", "Myeloid_Priming")
module_genes <- lapply(target_mods, function(m) str_to_title(module_list[[m]]))
names(module_genes) <- target_mods

surf_inflam  <- intersect(module_genes$Inflammatory_Response, surface_symbols)
surf_kras    <- intersect(module_genes$KRAS_Signaling, surface_symbols)
surf_myeloid <- intersect(module_genes$Myeloid_Priming, surface_symbols)

all_surface <- unique(c(surf_inflam, surf_kras, surf_myeloid))
all_surface <- all_surface[all_surface %in% rownames(HSC)]

gene_info <- data.frame(gene = all_surface) %>%
  mutate(Inflammatory = gene %in% surf_inflam,
         KRAS         = gene %in% surf_kras,
         Myeloid      = gene %in% surf_myeloid)

expr_data <- FetchData(HSC, vars = c(all_surface, "seurat_clusters"))
expr_long <- expr_data %>%
  pivot_longer(-seurat_clusters, names_to = "gene", values_to = "expr") %>%
  group_by(seurat_clusters, gene) %>%
  summarise(pct = mean(expr > 0) * 100, avg_expr = mean(expr), .groups = "drop") %>%
  left_join(gene_info, by = "gene")

gene_order <- expr_long %>% filter(seurat_clusters == "0") %>% arrange(desc(pct)) %>% pull(gene)
expr_long$gene <- factor(expr_long$gene, levels = gene_order)

mat_left <- gene_info %>%
  pivot_longer(-gene, names_to = "Module", values_to = "present") %>%
  mutate(Module = factor(Module, levels = c("Inflammatory", "KRAS", "Myeloid")),
         gene = factor(gene, levels = gene_order))

p_left <- ggplot(mat_left, aes(x = Module, y = gene, fill = present)) +
  geom_tile(color = "white", size = 0.6) +
  scale_fill_manual(values = c("TRUE" = "grey60", "FALSE" = "white"), guide = "none") +
  labs(x = NULL, y = NULL, title = "Module") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11),
        axis.text.y = element_text(face = "italic", size = 10),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.ticks = element_blank())

p_right <- ggplot(expr_long, aes(x = seurat_clusters, y = gene)) +
  geom_hline(yintercept = seq_len(length(all_surface)), color = "grey85", size = 0.3) +
  geom_point(aes(size = pct, fill = avg_expr), shape = 21, color = "black", stroke = 0.2) +
  scale_size_continuous(range = c(2, 8), name = "% expressed") +
  scale_fill_gradientn(colors = c("white", "#9f2b39"), name = "Avg expr") +
  labs(x = NULL, y = NULL, title = "Expression") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "right")

p_surface <- p_left + p_right + plot_layout(widths = c(0.52, 1)) +
  plot_annotation(title = "Cell surface markers from three functional modules",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))
ggsave("Fig_SurfaceMarkers_Module_Expression.pdf", p_surface,
       width = 8, height = max(6, length(all_surface) * 0.42), device = "pdf")

# ================================================================
# 图D：Cd69 表达小提琴图
# ================================================================
p_cd69_vln <- VlnPlot(HSC, features = "Cd69", group.by = "seurat_clusters",
                      pt.size = 0.2, cols = hsc_cluster_cols) +
  stat_compare_means(comparisons = list(c("0", "1"), c("0", "2"), c("1", "2")),
                     method = "wilcox.test", label = "p.format", size = 3.5) +
  labs(x = "Cluster", y = "Cd69 Expression",
       title = "Cd69 expression across HSC subclusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "none")
ggsave("Fig_Cd69_Violin.pdf", p_cd69_vln, width = 5, height = 5, device = "pdf")

# ================================================================
# 图E：细胞比例（WT vs Kras 各亚群占比）
# ================================================================
cell_tab <- table(HSC$orig.ident, HSC$seurat_clusters)
cell_ratio <- prop.table(cell_tab, margin = 1)
prop_df <- as.data.frame(cell_ratio)
names(prop_df) <- c("Sample", "Cluster", "Proportion")
prop_df$Sample <- factor(prop_df$Sample, levels = c("WT", "Kras"))
prop_df$Cluster <- factor(prop_df$Cluster, levels = c("0", "1", "2"))

p_prop <- ggplot(prop_df, aes(Sample, Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.3) +
  scale_fill_manual(values = hsc_cluster_cols) +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
  labs(x = NULL, y = "Proportion", title = "HSC Subcluster Proportion") +
  theme_clean() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Fig_CellProportion.pdf", p_prop, width = 5, height = 5, device = "pdf")

# ================================================================
# 5. 拟时序分析
# ================================================================
rd <- Embeddings(HSC, "umap")
cl <- HSC$seurat_clusters
sds <- slingshot(rd, cl, start.clus = "1")
HSC$pseudotime <- slingPseudotime(sds)[, 1]

curves <- slingCurves(sds)
curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
  df <- as.data.frame(curves[[i]]$s)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$lineage <- i
  df
}))
umap_df <- as.data.frame(rd) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  mutate(cluster = cl, pseudotime = HSC$pseudotime)

# 轨迹图（亚群着色）
p_traj_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cluster), size = 0.4, alpha = 0.6) +
  scale_color_manual(values = hsc_cluster_cols) +
  geom_path(data = curve_df, aes(group = lineage), color = "grey50", size = 0.9, alpha = 0.5,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed", angle = 20)) +
  labs(title = "HSC Subclusters on Trajectory") +
  theme_dr(xlength = 0.15, ylength = 0.15) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Fig_Trajectory_by_Cluster.pdf", p_traj_cluster, width = 8, height = 7, device = "pdf")

# 轨迹图（拟时间着色）
p_traj_pseudotime <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = pseudotime), size = 0.4) +
  scale_color_viridis_c(option = "inferno", name = "Pseudotime") +
  geom_path(data = curve_df, aes(group = lineage), color = "grey50", size = 0.8, alpha = 0.5,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed", angle = 20)) +
  labs(title = "Pseudotime along HSC Trajectory") +
  theme_dr(xlength = 0.15, ylength = 0.15) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Fig_Trajectory_by_Pseudotime.pdf", p_traj_pseudotime, width = 8, height = 7, device = "pdf")

# 拟时序密度分布（可单独用）
p_pseudotime_dist <- ggplot(umap_df, aes(x = pseudotime, fill = cluster)) +
  geom_density(alpha = 0.6, size = 0.3, color = "black") +
  scale_fill_manual(values = hsc_cluster_cols, name = "Cluster") +
  labs(x = "Pseudotime", y = "Density",
       title = "HSC subcluster distribution along pseudotime") +
  theme_clean() +
  theme(legend.position = "right")
ggsave("Fig_Pseudotime_Distribution.pdf", p_pseudotime_dist, width = 5, height = 4.5, device = "pdf")

# 模块评分与 CD69 沿拟时序的动态
dynamic_data <- FetchData(HSC, vars = c(
  "KRAS_Signaling_Score", "Inflammatory_Response_Score",
  "Myeloid_Priming_Score", "Cd69", "pseudotime")) %>% drop_na()

dynamic_long <- dynamic_data %>%
  pivot_longer(-pseudotime, names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = recode(Variable,
                           KRAS_Signaling_Score = "KRAS Signaling",
                           Inflammatory_Response_Score = "Inflammatory Response",
                           Myeloid_Priming_Score = "Myeloid Priming",
                           Cd69 = "CD69")) %>%
  mutate(Variable = factor(Variable, levels = c("KRAS Signaling", "Inflammatory Response", "Myeloid Priming", "CD69")))

var_colors <- c("KRAS Signaling" = "#9f2b39", "Inflammatory Response" = "#296097",
                "Myeloid Priming" = "#d6873b", "CD69" = "grey30")

p_dynamics <- ggplot(dynamic_long, aes(pseudotime, Value, color = Variable, fill = Variable)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.1, alpha = 0.15) +
  scale_color_manual(values = var_colors, guide = "none") +
  scale_fill_manual(values = var_colors, guide = "none") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(x = "Pseudotime", y = NULL, title = "Module scores & CD69 along pseudotime") +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Fig_Dynamics_Modules_Cd69.pdf", p_dynamics, width = 7, height = 9, device = "pdf")

# ================================================================
# 图F：拟时序分析：密度 → CD69 → 功能模块(合并) → 关键 TF(合并)
# ================================================================

# ---- TF 活性计算（若已存在可跳过） ----
tf_targets <- list(
  Nfkb1 = c("Tnf", "Il1b", "Il6", "Ccl4", "Cxcl10", "Icam1", "Vcam1"),
  Fos   = c("Junb", "Egr1", "Fosb", "Ccl4", "Il6", "Cxcl10"),
  Stat3 = c("Socs3", "Il6", "Ccl4", "Bcl2", "Myc", "Junb"),
  Jun   = c("Fos", "Egr1", "Junb", "Ccl4", "Il1b", "Tnf")
)
old_tf_cols <- grep("_TF$", colnames(HSC@meta.data), value = TRUE)
if (length(old_tf_cols)) HSC@meta.data <- HSC@meta.data[, !colnames(HSC@meta.data) %in% old_tf_cols]
for (tf in names(tf_targets)) {
  HSC <- AddModuleScore(HSC, features = tf_targets[tf],
                        name = paste0(tf, "_TF"), ctrl = min(100, nrow(HSC)))
  colnames(HSC@meta.data)[ncol(HSC@meta.data)] <- paste0(tf, "_TF")
}

# ---- 密度图 ----
p_density <- ggplot(umap_df, aes(x = pseudotime, fill = cluster)) +
  geom_density(alpha = 0.6, size = 0.3, color = "black") +
  scale_fill_manual(values = hsc_cluster_cols, name = "Cluster") +
  labs(x = NULL, y = "Density") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ---- CD69 单独一行 ----
cd69_traj <- FetchData(HSC, vars = c("Cd69", "pseudotime")) %>% drop_na()
p_cd69 <- ggplot(cd69_traj, aes(pseudotime, Cd69)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
              color = "#9f2b39", fill = "#9f2b39", size = 1, alpha = 0.2) +
  labs(x = NULL, y = "CD69") +
  theme_clean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ---- 功能模块合并（三条曲线同框） ----
module_dynamic <- FetchData(HSC, vars = c(
  "KRAS_Signaling_Score", "Inflammatory_Response_Score",
  "Myeloid_Priming_Score", "pseudotime"
)) %>% drop_na()

module_long <- module_dynamic %>%
  pivot_longer(-pseudotime, names_to = "Module", values_to = "Score") %>%
  mutate(Module = recode(Module,
                         KRAS_Signaling_Score = "KRAS Signaling",
                         Inflammatory_Response_Score = "Inflammatory Response",
                         Myeloid_Priming_Score = "Myeloid Priming"),
         Module = factor(Module, levels = c("KRAS Signaling", "Inflammatory Response", "Myeloid Priming")))

mod_colors <- c("KRAS Signaling" = "#9f2b39", "Inflammatory Response" = "#296097", "Myeloid Priming" = "#d6873b")

p_modules <- ggplot(module_long, aes(pseudotime, Score, color = Module, fill = Module)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.1, alpha = 0.15) +
  scale_color_manual(values = mod_colors) +
  scale_fill_manual(values = mod_colors) +
  labs(x = NULL, y = "Module Score", title = "Functional modules") +
  theme_clean() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold"))

# ---- 关键 TF 合并（四条曲线同框） ----
tf_names <- paste0(c("Nfkb1", "Fos", "Stat3", "Jun"), "_TF")
tf_traj_data <- FetchData(HSC, vars = c(tf_names, "pseudotime")) %>% drop_na()

tf_traj_long <- tf_traj_data %>%
  pivot_longer(-pseudotime, names_to = "TF", values_to = "Activity") %>%
  mutate(TF = str_remove(TF, "_TF"),
         TF = factor(TF, levels = c("Nfkb1", "Fos", "Stat3", "Jun")))

tf_colors <- c("Nfkb1" = "#9f2b39", "Fos" = "#c65341", "Stat3" = "#d6873b", "Jun" = "#e67e22")

p_tf <- ggplot(tf_traj_long, aes(pseudotime, Activity, color = TF, fill = TF)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.1, alpha = 0.15) +
  scale_color_manual(values = tf_colors) +
  scale_fill_manual(values = tf_colors) +
  labs(x = "Pseudotime", y = "TF Activity", title = "Transcription factors") +
  theme_clean() +
  theme(legend.position = "right",
        plot.title = element_text(size = 10, face = "bold"))

# ---- 统一 x 轴范围并拼接 ----
x_range <- range(HSC$pseudotime, na.rm = TRUE)
p_density <- p_density + xlim(x_range)
p_cd69    <- p_cd69 + xlim(x_range)
p_modules <- p_modules + xlim(x_range)
p_tf      <- p_tf + xlim(x_range)

p_trajectory_all <- p_density / p_cd69 / p_modules / p_tf +
  plot_layout(heights = c(0.2, 0.25, 0.35, 0.35)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("Fig_Pseudotime_Combined_Compact.pdf", p_trajectory_all, width = 8, height = 8, device = "pdf")

# ================================================================
# 附图：模块评分与拟时序的相关性棒棒图
# ================================================================
library(corrplot)

# 提取所有模块评分 + pseudotime
score_vars <- c(paste0(names(module_list), "_Score"), "pseudotime")
cor_data <- FetchData(HSC, vars = score_vars) %>% drop_na()

# 计算 Spearman 相关系数矩阵
cor_mat <- cor(cor_data, method = "spearman")

# 提取各模块与 pseudotime 的相关系数
pt_cor <- cor_mat["pseudotime", grep("_Score$", colnames(cor_mat)), drop = FALSE]
pt_cor_df <- data.frame(
  Module = str_remove(colnames(pt_cor), "_Score$") %>% str_replace_all("_", " "),
  Correlation = as.numeric(pt_cor)
) %>%
  mutate(Module = factor(Module, levels = Module[order(Correlation)]))  # 按相关系数排序

# 绘制水平条形图（棒棒图）
p_cor_bar <- ggplot(pt_cor_df, aes(x = Correlation, y = Module, fill = Correlation)) +
  geom_col(width = 0.6, color = "black", size = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Correlation),
                hjust = ifelse(Correlation > 0, -0.2, 1.2)),
            size = 3.5) +
  scale_fill_gradient2(low = "#92b8da", mid = "white", high = "#9f2b39",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = "Spearman Correlation with Pseudotime",
       y = NULL,
       title = "Module–Pseudotime Association") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.x = element_line(color = "grey90", size = 0.3))

ggsave("Fig_S1_Module_Correlation_Barplot.pdf", p_cor_bar, width = 6, height = 4, device = "pdf")




library(tidyverse)
library(ggrepel)

# ---- 数据准备 ----
jmml <- tribble(
  ~celltype, "proportion",
  "Adipo-CAR", 2.53164557,
  "Arteriolar ECs", 1.265822785,
  "CMP", 2.53164557,
  "Ery", 11.39240506,
  "Granolucyte", 12.65822785,
  "Macrophage", 44.30379747,
  "Monocyte", 5.063291139,
  "MPP", 1.265822785,
  "NK", 2.53164557,
  "Sinusoidal ECs", 7.594936709,
  "Stromal fibro.", 3.797468354,
  "T cells", 5.063291139
) %>% mutate(group = "JMML")

wt <- tribble(
  ~celltype, "proportion",
  "Adipo-CAR", 2.127659574,
  "Arteriolar ECs", 3.546099291,
  "B cell", 9.929078014,
  "CLP", 2.127659574,
  "CMP", 1.063829787,
  "Ery", 2.127659574,
  "Granolucyte", 10.28368794,
  "LMPP", 0.354609929,
  "Macrophage", 40.78014184,
  "Monocyte", 2.127659574,
  "MPP", 1.418439716,
  "NK", 0.709219858,
  "Sinusoidal ECs", 6.737588652,
  "Stromal fibro.", 6.028368794,
  "T cells", 10.63829787
) %>% mutate(group = "WT")

# 15 色向量
mycol15 <- c(
  "#9f2b39", "#7d6450", "#5e7f66", "#409079", "#4aa2aa",
  "#52a5c1", "#b46754", "#c65341", "#d6873b", "#c1b47d",
  "#92b8da", "#7ca1c2", "#564b84", "#5e6ca0", "#347852"
)

# ---- 添加百分比列和中点位置 ----
add_label <- function(df){
  df <- df %>%
    mutate(perc = proportion / sum(proportion) * 100) %>%
    arrange(desc(celltype)) %>%
    mutate(ypos = cumsum(perc) - 0.5*perc)
  return(df)
}

jmml <- add_label(jmml)
wt <- add_label(wt)

# ---- JMML 饼图（百分比保留1位小数，引线标注） ----
p_jmml <- ggplot(jmml, aes(x = 1, y = perc, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, color = "#222222") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = mycol15) +
  geom_text_repel(
    aes(y = ypos, label = paste0(round(perc, 0), "%")),  # 保留1位小数
    nudge_x = 1.2,
    segment.color = "grey50",
    segment.size = 0.5,
    size = 3,
    direction = "y",
    hjust = 0
  ) +
  labs(title = "JMML HSC Surrounding Cell Composition", fill = "Cell type") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  )

ggsave("JMML_HSC_pie_label_line.pdf", plot = p_jmml, width = 7, height = 7)

# ---- WT 饼图（百分比保留1位小数，引线标注） ----
p_wt <- ggplot(wt, aes(x = 1, y = perc, fill = celltype)) +
  geom_bar(stat = "identity", width = 1, color = "#222222") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = mycol15) +
  geom_text_repel(
    aes(y = ypos, label = paste0(round(perc, 0), "%")),  # 保留1位小数
    nudge_x = 1.2,
    segment.color = "grey50",
    segment.size = 0.5,
    size = 3,
    direction = "y",
    hjust = 0
  ) +
  labs(title = "WT HSC Surrounding Cell Composition", fill = "Cell type") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  )

ggsave("WT_HSC_pie_label_line.pdf", plot = p_wt, width = 7, height = 7)

# ===============================
# Spatial HSC Niche & Cell Interaction
# Full Analysis Pipeline with Figure Annotations (Fig4A–F)
# ===============================

## ===============================
## 0. Environment & dependencies
## ===============================
suppressMessages({
  library(Seurat)
  library(spacexr)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(Rcpp)
  library(tidyr)
  library(stringr)
  library(reshape2)
  library(CellChat)
  library(ktplots)
})

set.seed(123)

## ===============================
## 1. Parameters
## ===============================
scmeta      <- "../sc_meta.txt"
anno        <- "../ref_cell_anno"
stmat       <- "./mtx/"
od          <- "./outdir"
number      <- 2
min_cells   <- 5
min_features<- 100
point_size  <- 0.01
HSC_name    <- "HSC"
radius      <- 100  # HSC neighborhood radius
if(!dir.exists(od)) dir.create(od, recursive = TRUE)

## ===============================
## 2. Single-cell reference (RCTD reference)
## ===============================
# Fig4A: Prepare reference for spatial annotation of WT and Kras samples
sc_counts <- read.table(scmeta, header = TRUE, row.names = 1, check.names = FALSE)
sc_nUMI   <- colSums(sc_counts)
cellType <- read.table(anno, header = FALSE, sep = "\t", check.names = FALSE)
colnames(cellType) <- c("barcode", "cell_type")

cellType$cell_type <- as.factor(cellType$cell_type)
cellType <- cellType %>%
  filter(cell_type %in% names(table(cell_type)[table(cell_type) > 25]))
cell_types <- setNames(as.character(cellType$cell_type), cellType$barcode)
cell_types <- as.factor(cell_types)

sc_counts <- sc_counts[, colnames(sc_counts) %in% names(cell_types)]
sc_nUMI   <- sc_nUMI[names(cell_types)]

reference <- Reference(sc_counts, cell_types, sc_nUMI)

## ===============================
## 3. Sparse to dense matrix conversion
## ===============================
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  IntegerMatrix mat(nrows, ncols);
  for (int i = 0; i < z.size(); i++){
    mat(rp[i], cp[i]) = z[i];
  }
  return mat;
}
')

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x) - 1, mat@p[-1])
  tmp <- asMatrix(row_pos, col_pos, mat@x,
                  mat@Dim[1], mat@Dim[2])
  rownames(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  tmp
}

## ===============================
## 4. Load spatial data
## ===============================
# Fig4A: Spatial coordinates
coords <- read.table(gzfile(paste0(stmat,"/barcodes_pos.tsv.gz")),
                     sep="\t", header=FALSE)
colnames(coords) <- c("barcode","x","y")
coords$y <- -coords$y
rownames(coords) <- coords$barcode
coords$barcode <- NULL

expr <- Read10X(stmat, gene.column = number)

sp_obj <- CreateSeuratObject(expr,
                             assay = "Spatial",
                             min.cells = min_cells,
                             min.features = min_features)

sp_counts <- as_matrix(sp_obj@assays$Spatial@counts)
sp_nUMI   <- colSums(sp_counts)

puck <- SpatialRNA(coords, sp_counts, sp_nUMI)

## ===============================
## 5. RCTD integration
## ===============================
# Fig4A: Annotate spatial spots using single-cell reference
myRCTD <- create.RCTD(puck, reference,
                      max_cores = 8,
                      CELL_MIN_INSTANCE = 20)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
saveRDS(myRCTD, file = paste0(od,"/RCTD.rds"))

## ===============================
## 6. Organize spatial annotation
## ===============================
# Fig4A: Save spatial cell-type annotations for WT and Kras
res_df <- myRCTD@results$results_df
valid_bc <- rownames(res_df[res_df$spot_class!="reject" & puck@nUMI>=1,])
anno_df <- puck@coords[valid_bc,]
anno_df$cell_type <- res_df[valid_bc,"first_type"]

write.table(
  anno_df %>% rownames_to_column("barcode"),
  paste0(od,"/Spatial_CellType.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

## ===============================
## 7. Extract HSC neighbors
## ===============================
# Fig4B: Identify cells located near HSCs in WT and Kras
hsc_coords <- anno_df %>% filter(cell_type == HSC_name)
dist_mat <- as.matrix(dist(anno_df[,c("x","y")]))
neighbor_idx <- which(apply(dist_mat[, rownames(hsc_coords)], 1,
                            function(x) any(x <= radius)))
neighbor_cells <- anno_df[neighbor_idx, ]
neighbor_cells <- neighbor_cells %>% filter(cell_type != HSC_name)

prop_df <- neighbor_cells %>%
  count(cell_type) %>%
  mutate(prop = n / sum(n))

write.table(prop_df,
            paste0(od,"/HSC_neighbor_cell_proportion.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

## ===============================
## 8. Donut plots
## ===============================
# Fig4C: WT HSC neighborhood composition
p_donut_WT <- ggplot(prop_df, aes(x=2, y=prop, fill=cell_type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y") +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title="HSC Niche Composition (WT)")

ggsave(paste0(od,"/HSC_neighbor_donut_WT.pdf"), p_donut_WT, width=5, height=5)

# Fig4D: Kras HSC neighborhood composition (repeat above with Kras dataset)
# Assuming neighbor_cells_Kras is prepared similarly
# p_donut_Kras <- ...
# ggsave(...)

## ===============================
## 9. CellPhoneDB DotPlot function
## ===============================
# Fig4E–F: HSC-centered ligand–receptor interactions
cellphoneDB_Dotplot <- function(pvals.data, means.data, key,
                                target.cells_1, target.cells_2=NA,
                                p.cutoff=0.05){
  
  colnames(pvals.data) <- str_replace_all(colnames(pvals.data),"\\.","_")
  colnames(means.data) <- str_replace_all(colnames(means.data),"\\.","_")
  
  kp <- Reduce(`|`, lapply(target.cells_1, grepl, x=colnames(pvals.data)))
  pos <- which(kp)
  
  pvals <- pvals.data[,c(1,2,5,6,8,9,pos)]
  means <- means.data[,c(1,2,5,6,8,9,pos)]
  
  pvals <- pvals[rowSums(pvals[,7:ncol(pvals)] < p.cutoff)>0,]
  means <- means[means$id_cp_interaction %in% pvals$id_cp_interaction,]
  
  df <- merge(
    melt(pvals, id.vars="interacting_pair"),
    melt(means, id.vars="interacting_pair"),
    by=c("interacting_pair","variable")
  )
  
  ggplot(df, aes(variable, interacting_pair)) +
    geom_point(aes(size=-log10(value.x+1e-4),
                   color=log2(value.y+1))) +
    scale_colour_gradientn(colors=c("#3A5978","#F6B31D","#DA2328")) +
    theme_bw() +
    labs(x="", y="", title=paste("CellPhoneDB:", key))
}

## ===============================
## 10. CellPhoneDB networks and dotplots
## ===============================
# Fig4E: Interaction network
setwd("/data/project/E297/WT/")
df.net <- read.table("count_network.txt", header=TRUE, sep="\t")
df.net <- spread(df.net, TARGET, count)
rownames(df.net) <- df.net$SOURCE
df.net <- as.matrix(df.net[,-1])

pvals_stat <- read.delim("pvalues.txt", check.names=FALSE)
means_stat <- read.delim("means.txt", check.names=FALSE)

pdf("CellPhoneDB_net_circle.pdf")
netVisual_circle(df.net, weight.scale=TRUE, label.edge=FALSE)
dev.off()

# Fig4F: HSC-centered DotPlot
pdf("CellPhoneDB_DotPlot_HSC.pdf", width=8, height=10)
cellphoneDB_Dotplot(pvals_stat, means_stat,
                    key=HSC_name,
                    target.cells_1=c(HSC_name))
dev.off()

# 需要检查的包（请根据实际使用情况微调）
pkgs <- c(
  "Seurat", "harmony", "tidyverse", "patchwork", "dplyr", "ggplot2",
  "ggpubr", "cowplot", "stringr", "GSVA", "pheatmap", "slingshot",
  "monocle", "org.Mm.eg.db", "AnnotationDbi", "viridis", "mgcv",
  "gplots", "corrplot"
)

# 提取版本信息
pkg_info <- do.call(rbind, lapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    data.frame(Package = p, Version = as.character(packageVersion(p)),
               stringsAsFactors = FALSE)
  } else {
    data.frame(Package = p, Version = "NOT INSTALLED", stringsAsFactors = FALSE)
  }
}))

# 输出为表格（控制台）
print(pkg_info, row.names = FALSE)

# 若需保存为 CSV
write.csv(pkg_info, "package_versions.csv", row.names = FALSE)
