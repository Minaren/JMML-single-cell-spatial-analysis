suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(stringr))

load(pbmc_filt)
allmarkers <- FindAllMarkers(pbmc_filt,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
write.csv(allmarkers,'allmarkers_mouse_CD45+CD3+CD4+.csv') 
table(Idents(pbmc_filt))
CD4Tcells <- subset(pbmc_filt, seurat_clusters %in% c(1, 2))


load("D:/aJMML/scRNA1/mouse_CD4Tcells.RData")
table(CD4Tcells$orig.ident)
CD4Tcells$celltype <- "Other"  # 先初始化
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 1] <- "Treg"
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 2] <- "Tconv"






####Fug5B####
library(ggplot2)
library(grid)

# 基础图 + facet
pumap_facet <- ggplot(meta, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 1) +
  facet_wrap(~orig.ident, nrow = 1) +
  scale_color_manual(values = mycol) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# 获取坐标范围，用于箭头位置
xrange <- range(meta$UMAP_1)
yrange <- range(meta$UMAP_2)
arrow_len_x <- diff(xrange) * 0.1  # 箭头长度可调
arrow_len_y <- diff(yrange) * 0.1

# 添加坐标轴箭头（从左下角画出 X、Y 轴箭头）
pumap_facet_arrows <- pumap_facet +
  geom_segment(aes(x = xrange[1], y = yrange[1],
                   xend = xrange[1] + arrow_len_x, yend = yrange[1]),
               arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
               inherit.aes = FALSE) +
  geom_segment(aes(x = xrange[1], y = yrange[1],
                   xend = xrange[1], yend = yrange[1] + arrow_len_y),
               arrow = arrow(length = unit(0.15, "inches"), type = "closed"),
               inherit.aes = FALSE)

# 保存图像
png("umap_facet_celltype_arrow.png", width = 10, height = 5, units = "in", res = 300)
print(pumap_facet_arrows)
dev.off()

pdf("umap_facet_celltype_arrow.pdf", width = 10, height = 5)
print(pumap_facet_arrows)
dev.off()












####Fig5C####

library(ggplot2)
table(CD4Tcells$orig.ident)#查看各组细胞数
prop.table(table(Idents(CD4Tcells)))
table(Idents(CD4Tcells), CD4Tcells$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(CD4Tcells), CD4Tcells$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
CD4Tcells_col <- c("#9f2b39","#92b8da")
# 计算颜色数量
colourCount = length(unique(Cellratio$Var1))
#两组位置互换
Cellratio$Var2 <- factor(Cellratio$Var2, levels = rev(levels(Cellratio$Var2)))
# 绘制柱状图并应用配色
ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.2, size = 0.5, colour = '#222222') + 
  scale_fill_manual(values = CD4Tcells_col) +  # 使用自定义颜色
  labs(x = 'Sample', y = 'Ratio', fill = "subpopulation") +  # 修改图例标题为 "subpopulation"
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

# 保存图片
ggsave('percentage_CD4Tcells.png', width = 7, height = 6, dpi = 500)
ggsave('percentage_CD4Tcells.pdf', width = 7, height = 6, dpi = 500)



####Fig5D####
library(Seurat)
library(patchwork)

# Treg核心标志性基因（在您的列表中）
treg_core_markers <- c(
  "Foxp3",       # 主调控转录因子
  "Ctla4",       # 关键抑制受体
  "Tnfrsf4",     # OX40
  "Il2ra",       #Cd25
  "Il10",       # 介导免疫耐受的核心细胞因子
  "Cd69"
)

# 生成每个基因的 FeaturePlot 并移除坐标轴和图例
plots <- lapply(treg_core_markers, function(gene) {
  FeaturePlot(CD4Tcells, features = gene, cols = c("lightgrey", "#9f2b39")) + NoAxes() + NoLegend()
})

# 使用 patchwork 拼接图
p <- wrap_plots(plots, ncol = 3)

# 保存图像
ggsave("FeaturePlot_Treg.pdf", plot = p, width = 15, height = 10)

FeaturePlot(CD4Tcells, features = treg_core_markers, cols = c("lightgrey", "#9f2b39"))
ggsave("FeaturePlot_Treg.pdf", width = 10, height = 15)




FeaturePlot(Treg, features = "Cd69", cols = c("lightgrey", "#9f2b39"), split.by = "orig.ident")






####Fig5E####
Treg <- subset(CD4Tcells, seurat_clusters == 1)
FeaturePlot(Treg, features = "Cd69", cols = c("lightgrey", "#9f2b39"), split.by = "orig.ident")#Cd2和Cd5效果也很好
ggsave("FeaturePlot_Treg_Cd691.pdf", width = 10, height = 5)
FeaturePlot(
  Treg,
  features = c("Cd69"),
  cols = c("#f0f0f0", "#9f2b39"),  # 柔和的渐变色
  ncol = 2,
  pt.size = 0.7,                   # 点的大小
  split.by = "orig.ident"
) & 
  theme_bw() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


FeaturePlot(
  CD4Tcells,
  features = c("Foxp3", "Ctla4", "Fgl2", "Tnfrsf4"),
  cols = c("lightgrey", "#9f2b39"),
  ncol = 2,
  pt.size = 0.8,      # 点大小
  combine = TRUE
) & 
  theme_bw() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

ggsave("FeaturePlot_Treg1.pdf", width = 10, height = 10)
# 然后在作图中使用 celltype
DimPlot(CD4Tcells, group.by = "celltype", label = TRUE)
DimPlot(CD4Tcells, group.by = "seurat_clusters", label = TRUE)
save(CD4Tcells,file = "mouse_CD4Tcells.RData")

# 整合元数据和UMAP坐标
meta <- CD4Tcells@meta.data
umap_coordinates <- Embeddings(CD4Tcells, "umap")
meta <- cbind(meta, umap_coordinates)

# 定义 celltype 的颜色
col_df <- data.frame(name = unique(meta$celltype)) %>% arrange(name)
mycol <- c("#92b8da", "#9f2b39")  # 根据 celltype 数量调整
mycol <- setNames(mycol, col_df$name)

# 设置 orig.ident 顺序（如对照在左，实验在右）
meta$orig.ident <- factor(meta$orig.ident, levels = c("WT", "CK")) # 示例

# 基础图 + facet
pumap_facet <- ggplot(meta, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 1) +
  facet_wrap(~orig.ident, nrow = 1) +
  scale_color_manual(values = mycol) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# 添加坐标轴箭头
pumap_facet_arrows <- pumap_facet +
  theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")
  )

png("umap_facet_celltype_arrow.png", width = 10, height = 5, units = "in", res = 300)
print(pumap_facet_arrows)
dev.off()

pdf("umap_facet_celltype_arrow.pdf", width = 10, height = 5)
print(pumap_facet_arrows)
dev.off()









####Fig5F分析####
####三、信号通路分析####
rm(list = ls())
####安装及加载R包
# 正确筛选 Treg
load("D:/aJMML/scRNA1/mouse_CD4Tcells.RData")
Treg <- subset(CD4Tcells, subset = celltype == "Treg")
Idents(Treg) <- "orig.ident"

# DEGs
markers <- FindMarkers(Treg, ident.1 = "CK", ident.2 = "WT", only.pos = FALSE, logfc.threshold = 0.25)
markers <- markers %>% rownames_to_column("gene") %>% filter(p_val_adj < 0.05)

# ID转换
gene_convert <- bitr(markers$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
markers <- inner_join(markers, gene_convert, by = c("gene" = "SYMBOL"))

# GO富集
go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID", ont="BP", OrgDb=org.Mm.eg.db)
p1 <- dotplot(go.results, showCategory = 20, label_format = 10000)
ggsave("GO_dotplot_BP_Treg_mouse.pdf", plot = p1, width = 10, height = 8)
write.csv(go.results@result, file = "Treg_GO_results_BP.csv", row.names = TRUE)




## 1. 加载包
library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(limma)
library(pheatmap)

## 2. 加载数据并设定 CD69 状态
load("D:/aJMML/scRNA1/CD4Tcell.RData")
Treg <- subset(CD4Tcells, celltype=="Treg")
Idents(Treg) <- "orig.ident"
## 3. 提取表达矩阵 & 元信息
Treg_counts <- as.matrix(Treg@assays$RNA@counts)
meta <- Treg@meta.data
cells_to_use <- WhichCells(Treg, idents = c("CK", "WT"))
Treg_counts_filtered <- Treg_counts[, cells_to_use]

## 4. 获取 GO:BP 基因集
mouse_GO_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
mouse_GO_bp_Set <- split(mouse_GO_bp$gene_symbol, mouse_GO_bp$gs_name)

## 5. GSVA 分析
Treg_gsva <- gsva(expr = Treg_counts_filtered,
                  gset.idx.list = mouse_GO_bp_Set,
                  kcdf = "Poisson",
                  parallel.sz = 4)
write.table(Treg_gsva, "Treg_gsva.xls", sep = "\t")

## 6. 差异分析
group <- factor(meta[cells_to_use, "orig.ident"])
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(Treg_gsva, design)
cont.matrix <- makeContrasts(CK_vs_WT = CK - WT, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, adjust = "fdr", number = Inf)
diff$ID <- rownames(diff)
write.csv(diff, "CD69_GSVA_diff.csv")

## 7. 自定义感兴趣的通路
# 细胞增殖相关

diff <- read.csv("CD69_GSVA_diff.csv", header = TRUE, stringsAsFactors = FALSE)

# 查看前几行确认内容
head(diff)

# 1. 从 diff 中筛选感兴趣的免疫与炎症通路
immune_inflammation_pathways <- c(
  "leukocyte aggregation",       # 白细胞聚集：炎症反应中白细胞在损伤或感染区域的集合
  "leukocyte migration",         # 白细胞迁移：指白细胞从血管向炎症部位的运动
  "leukocyte chemotaxis",        # 白细胞趋化：白细胞对化学信号（如趋化因子）的定向移动
  "leukocyte cell-cell adhesion",# 白细胞细胞间粘附：帮助白细胞贴附并穿过血管壁
  "cell chemotaxis",             # 细胞趋化：广义的细胞定向运动过程
  "chemotaxis",                  # 趋化反应：细胞对外界刺激的定向移动
  "taxis"                      # 趋向性运动：总称细胞或有机体对刺激方向移动的现象
)



# 2. 构建 plot_df：只保留这些通路的数据
plot_df <- diff[diff$ID %in% immune_inflammation_related, ]

# 添加模块分类列
plot_df$module <- "Immune_Inflammation"

# 3. 按 module 和 t 值排序，准备 factor 水平
plot_df <- plot_df[order(plot_df$module, plot_df$t), ]
plot_df$ID <- factor(plot_df$ID, levels = plot_df$ID)

# 4. 绘图
library(ggplot2)

p <- ggplot(plot_df, aes(ID, t, fill = module)) +
  geom_bar(stat = 'identity', alpha = 0.9, color = "grey30", width = 0.7) +
  coord_flip() +
  facet_wrap(~ module, scales = "free_y", ncol = 1) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold", color = "#4B4B4B"),
    panel.grid.major = element_blank(),    # 去除主网格线
    panel.grid.minor = element_blank(),    # 去除次网格线
    panel.spacing = unit(1, "lines"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "GSVA t value", title = "Immune/Inflammation Related Pathways") +
  scale_fill_manual(values = c("Immune_Inflammation" = "#9f2b39"))


####Fig5F####
ggsave("gsva_barplot_immune_inflammation.pdf", plot = p, width = 12, height = 8)


