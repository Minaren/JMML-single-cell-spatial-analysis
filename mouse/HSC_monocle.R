##提取HSC，HSC
#处理注释过的数据
load("D:/aJMML/scRNA1/sce_har_anno.RData")
table(Idents(sce.harm))
table(sce.harm@meta.data$seurat_clusters) 
table(sce.harm@meta.data$orig.ident) 
str(sce.harm)#观察数据信息

##取子集
levels(Idents(sce.harm))
sce_HSC = sce.harm[, Idents(sce.harm) %in% c( "HSC"  )] #sce.s是sce.subset的简写
levels(Idents(sce_HSC))
saveRDS(sce_HSC, file = "sce_HSC.rds")

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
rm(list=ls())
#dir.create("Monocle3")
setwd("Monocle3")

##创建CDS对象并预处理数据
sce_HSC <- readRDS("sce_HSC.rds")
pbmc<-sce_HSC
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
cds <- align_cds(cds, alignment_group = "plate",
                 residual_model_formula_str = '~Size_Factor')

#3.2 可视化
#降维、聚类
cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
#plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")#seurat_clusters
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters")
ggsave("Reduction_HSC.pdf",width = 10, height = 10)
#可视化指定基因
ciliated_genes <- c( "Slamf1","Kit","Cd34","Scal1","Cd48","Procr11","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ggsave("ciliated.pdf",width = 10, height = 10)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

####4. 构建细胞轨迹####
#4.1 轨迹学习Learn the trajectory graph（使用learn_graph()函数）
## 识别轨迹
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Trajectory.pdf", plot = p, width = 8, height = 6)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)#展示分支节点
##手动选择root
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = TRUE,  label_branch_points = FALSE)
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


cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds <- order_cells(cds,root_pr_nodes=myselect(cds,select.classify = 'cao_cell_type',
                                              my_select = "Body wall muscle"))
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="cao_cell_type")
cds_3d_plot_obj
