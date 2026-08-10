#相关R包载入：
rm(list = ls())
library(Seurat)
library(dplyr)
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)

load("D:/aJMML/scRNA1/sce_har_anno.RData")
sce<-sce.harm
sce[['celltype']] <- Idents(sce)
dim(sce)
#Marker基因提取：
#markers <- FindAllMarkers(sce,logfc.threshold = 0.6,min.pct = 0.25,only.pos = T)
#markers<-read.csv("allmarkers_celltype_mouse.csv",header=T)
#head(markers)

#提取各亚群Top5Marker：
#sig_markers <- markers %>%group_by(cluster)%>%top_n(n = 5, wt = avg_log2FC)
#head(sig_markers)
sig_markers<-read.csv("top5_celltype_mouse.csv",header=T)

#各亚群平均表达量提取：
genes <- unique(sig_markers$gene)
aver_dt <- AverageExpression(sce,
                             features = genes,
                             group.by = 'celltype',
                             slot = 'data')
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt <- aver_dt[, ncol(aver_dt):1]
aver_dt[1:6,1:6]

#基础热图绘制：
pheatmap(as.matrix(aver_dt),
         scale = 'row',
         cluster_rows = FALSE,
         cluster_cols = FALSE)


#行列注释信息色块添加
#行注释：Top5marker对应celltype
gene_anno <- data.frame(gene_anno = sig_markers$cluster,
                        row.names = sig_markers$gene)
gene_anno <- gene_anno[rownames(gene_anno) %in% rownames(aver_dt), , drop = FALSE]
# 重新排序 gene_anno 以匹配 aver_dt 的行顺序
gene_anno <- gene_anno[rownames(aver_dt), , drop = FALSE]
rownames(aver_dt)
rownames(gene_anno)
#列注释：celltype
cell_anno <- data.frame(cell_anno = colnames(aver_dt),
                        row.names = colnames(aver_dt))
head(gene_anno);head(cell_anno)


#将行列注释信息添加到热图：
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno, #列注释色块添加
         annotation_row = gene_anno) #行注释色块添加

#热图美化：
#热图配色自定义：
mycol <- colorRampPalette(c("#5E3C99", "white", "#E66101"))(50)

#行列注释配色自定义：
#celltype_col <- c4a('10',9)
celltype_col = c(
  "#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da",
  "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84")
names(celltype_col) <- cell_anno$cell_anno

anno_col <- list(cell_anno = celltype_col,
                 gene_anno = celltype_col)
anno_col

#绘图：
pheatmap(as.matrix(aver_dt),
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = cell_anno,
         annotation_row = gene_anno,
         annotation_colors = anno_col, #注释配色
         color = mycol, #热图配色
         border_color = 'white') #描边颜色

#自定义配色：
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
#归一化：
aver_dtt <- t(scale(t(aver_dt)))

#添加行列注释：

#cols <- c4a('classic10light', 9)
cols<-c(
  "#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da",
  "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84")
names(cols) <- cell_anno$cell_anno

#列注释：
cell <- data.frame(colnames(aver_dtt))
colnames(cell) <- 'cell'

col_anno <- HeatmapAnnotation(df = cell,
                              show_annotation_name = F,
                              gp = gpar(col = 'white', lwd = 2),
                              col = list(cell = cols))

#行注释：
row_cols <- setNames(rep(cols, each = 5), rownames(aver_dtt))
row_cols

row_anno <- rowAnnotation(foo = anno_text(rownames(aver_dtt),
                                          location = 0,
                                          just = "left",
                                          gp = gpar(fill = row_cols,
                                                    col = "black",
                                                    fontface = 'italic'),
                                          width = max_text_width(rownames(aver_dtt))*1.25))

#热图绘制：
Heatmap(aver_dtt,
        name = 'expression',
        col = mycol2,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_names_side = 'top',
        column_names_rot = 60,
        column_names_gp = gpar(fontsize =20, fontface = 'italic'),  # 调整列名字体大小
        row_names_gp = gpar(fontsize = 20, fontface = 'italic'),  # 调整行名字体大小
        rect_gp = gpar(col = "white", lwd = 1.5),
        top_annotation = col_anno) + row_anno
