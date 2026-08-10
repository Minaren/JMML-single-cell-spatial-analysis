
levels(Idents(MPP))
MPP_KRAS<-subset(x = MPP, subset = orig.ident == "Kras")#sce.s是sce.subset的简写
levels(Idents(MPP_KRAS))

pbmc<-MPP_KRAS
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
ggsave("Reduction_MPP_KRAS.pdf",width = 10, height = 10)
#可视化指定基因
ciliated_genes <- c( "Slamf1","Kit","Cd34","Scal1","Cd48","Procr11","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
ggsave("ciliated_MPP_KRAS.pdf",width = 10, height = 10)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

####4. 构建细胞轨迹####
#4.1 轨迹学习Learn the trajectory graph（使用learn_graph()函数）
## 识别轨迹
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=TRUE,
           label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Trajectory_MPP_Kras.pdf", plot = p, width = 8, height = 6)
plot_cells(cds, color_cells_by = "seurat_clusters", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)#展示分支节点
##手动选择root
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = TRUE,  label_branch_points = FALSE)
ggsave("Trajectory_Pseudotime_MPP_KRAS.pdf", plot = p, width = 8, height = 6)
saveRDS(cds, file = "cds.rds")
