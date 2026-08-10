load("D:/aJMML/scRNA/sce_harm_tu.RData")
sce <- sce.harm
celltype = data.frame(ClusterID = 0:22,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0,2,4,5,8,9),2] = "Neutrophil"
celltype[celltype$ClusterID %in% c(1,13,14,17,18),2] = "Macrophage"
celltype[celltype$ClusterID %in% c(3),2] = 'CMP'
celltype[celltype$ClusterID %in% c(6),2] = 'CLP'
celltype[celltype$ClusterID %in% c(7,15),2] = 'HSPC'
celltype[celltype$ClusterID %in% c(10,21),2] = 'Ery' 
celltype[celltype$ClusterID %in% c(11),2] = "Pre_B" 
celltype[celltype$ClusterID %in% c(12),2] = 'Ery_prog'
celltype[celltype$ClusterID %in% c(16),2] = 'Basophil' 
celltype[celltype$ClusterID %in% c(19),2] = "T" 
celltype[celltype$ClusterID %in% c(20),2] = "Unknown"
celltype[celltype$ClusterID %in% c(22),2] = 'B' 

head(celltype)
celltype 
table(celltype$celltype)
sce.harm$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.harm@meta.data[which(sce.harm@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.harm@meta.data$celltype)

Idents(sce.harm) = sce.harm$celltype
Idents(sce.harm) = factor(Idents(sce.harm),levels = c('HSPC','CLP','CMP','Ery_prog','Pre_B','T','B','Neutrophil','Macrophage','Basophil','Ery','Unknown'))
sce.harm$celltype <- Idents(sce.harm)
Idents(sce.harm) <- factor(Idents(sce.harm),levels = rev(levels(Idents(sce.harm))))


library(RColorBrewer)
color_ct=c(brewer.pal(12, "Set3"),"#b3b3b3",
           brewer.pal(5, "Set1"),
           brewer.pal(3, "Dark2"),
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")


library(ComplexHeatmap)
library(jjAnno)
#if(!require(ggunchull))devtools::install_github("sajuukLyu/ggunchull", type = "source")

#devtools::install_local("D:/R/R-4.2.1/library/scRNAtoolVis/scRNAtoolVis-master.zip")
library(scRNAtoolVis)
clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'celltype',pSize = 0.1,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(color_ct,0.65)) + NoLegend() +
  scale_fill_manual(values = alpha(color_ct,0.65))
ggsave("umap.pdf",width = 10,height = 9,dpi = 500)
ggsave("umap.png",width = 10,height = 9,dpi = 500)

clusterCornerAxes(object = sce.harm,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("umap_indi.pdf",width = 5.5,height = 5,dpi = 500)
ggsave("umap_indi.png",width = 5.5,height = 5,dpi = 500)

saveRDS(sce.harm, file = "sce_try.rds")
sce_try<- readRDS("sce_try.rds")
pbmc<-sce_try
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
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "celltype")#seurat_clusters
ggsave("Reduction.pdf",width = 10, height = 10)
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
p = plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
               label_leaves=FALSE, label_branch_points=FALSE)
ggsave("Trajectory.pdf", plot = p, width = 8, height = 6)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,
           label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=1.5)#展示分支节点
ggsave("Trajectory_points.pdf", plot = p, width = 8, height = 6)
##手动选择root
# 解决order_cells(cds)报错"object 'V1' not found"
# rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
# colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
ggsave("Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)
saveRDS(cds, file = "cds.rds")
