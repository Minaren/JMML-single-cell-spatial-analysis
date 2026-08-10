load("D:/aJMML/scRNA/sce_har_anno.RData")
saveRDS(sce.harm, file = "pmbc.rds")

#安装Monocle 3
#保证自己的环境符合这些版本要求，不然报错很麻烦
#R >- 4.1.0，Bioconductor - 3.14,这样装上的monocle3 >= 1.2.7
#没有达到理想版本
####1. 安装####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
package.version("monocle3")

####2. 数据准备，创建CDS对象并进行降维。####
#注意：该数据集使用的是pbmc3k的数据集，由于pbmc都是分化成熟的免疫细胞，理论上并不存在直接的分化关系，因此不适合用来做拟时轨迹分析。这里仅作为学习演示。

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
rm(list=ls())
#dir.create("Monocle3")
setwd("Monocle3")

##创建CDS对象并预处理数据

pbmc <- readRDS("pbmc.rds")
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
####3. 预处理####
#3.1 标准化和PCA降维
#RNA-seq是使用PCA，如果是处理ATAC-seq的数据用Latent Semantic Indexing)

#⚠️preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#3.2 可视化
#umap降维
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
plot_cells(cds)
colnames(colData(cds))
#以之前的Seurat分群来添加颜色，和原有的Seurat分群对比
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('int.umap')
p = p1|p2
ggsave("Reduction_Compare.pdf", plot = p, width = 10, height = 5)

#可视化指定基因
ciliated_genes <- c( "Slamf1","Kit","Cd34","Scal1","Cd48","Procr11","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ggsave("ciliated.pdf",width = 10, height = 10)

#也可以使用tSNE降维
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")
plot_cells(cds, reduction_method="UMAP", color_cells_by="cell_type") 

####4. 构建细胞轨迹####
#4.1 轨迹学习Learn the trajectory graph（使用learn_graph()函数）
## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Trajectory.pdf", plot = p, width = 8, height = 6)
plot_cells(cds, color_cells_by = "cell_type", label_groups_by_cluster=FALSE,
           +                label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)

##手动选择root
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)
saveRDS(cds, file = "cds.rds")


####5. 差异表达分析####
#5.1 寻找拟时轨迹差异基因
#graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
#空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)
#5.2 挑选top10画图展示
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
#基因表达趋势图

p <- plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters", 
                              min_expr=0.5, ncol = 2)
ggsave("Genes_Jitterplot.pdf", plot = p, width = 8, height = 6)
#FeaturePlot图

p <- plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
                label_cell_groups=FALSE,  label_leaves=FALSE)
p$facet$params$ncol <- 5
ggsave("Genes_Featureplot.pdf", plot = p, width = 20, height = 8)

#寻找共表达基因模块

Track_genes <- read.csv("Trajectory_genes.csv")
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-1, cores = 6)
write.csv(gene_module, "Genes_Module.csv", row.names = F)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
ggsave("Genes_Module.pdf", plot = p, width = 8, height = 8)

#提取拟时分析结果返回seurat对象

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
pseudotime <- pseudotime[rownames(pbmc@meta.data)]
pbmc$pseudotime <- pseudotime
p = FeaturePlot(pbmc, reduction = "umap", features = "pseudotime")
# pseudotime中有无限值，无法绘图。
ggsave("Pseudotime_Seurat.pdf", plot = p, width = 8, height = 6)
saveRDS(pbmc, file = "sco_pseudotime.rds")