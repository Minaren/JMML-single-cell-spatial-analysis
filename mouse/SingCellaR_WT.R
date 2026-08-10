#if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager")
#BiocManager::install("AUCell")

#devtools::install_github('supatt-lab/SingCellaR')
#install.packages("AUCell")

rm(list=ls())
library(SingCellaR)
library(AUCell)
data_matrices_dir<-"D:/aJMML/JMML/WT"
Sample2<-new("SingCellaR")
Sample2@dir_path_10x_matrix<-data_matrices_dir
Sample2@sample_uniq_id<-"Sample2"

load_matrices_from_cellranger(Sample2,cellranger.version = 3)

Sample2

process_cells_annotation(Sample2,mito_genes_start_with ="mt-")

####QC####
plot_cells_annotation(Sample2,type="histogram")
plot_cells_annotation(Sample2,type="boxplot")
plot_UMIs_vs_Detected_genes(Sample2)
DoubletDetection_with_scrublet(Sample2)
table(Sample2@meta.data$Scrublet_type)

filter_cells_and_genes(Sample2,
                       min_UMIs=1000,
                       max_UMIs=60000,
                       min_detected_genes=500,
                       max_detected_genes=6000,
                       max_percent_mito=15,
                       genes_with_expressing_cells = 10,isRemovedDoublets = F)
normalize_UMIs(Sample2,use.scaled.factor = T)
remove_unwanted_confounders(Sample2,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(Sample2,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)
plot_variable_genes(Sample2)

runPCA(Sample2,use.components=50,use.regressout.data = T)
plot_PCA_Elbowplot(Sample2)
runUMAP(Sample2,dim_reduction_method = "pca",n.dims.use = 20,n.neighbors = 30,
        uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(Sample2,feature = "UMI_count",point.size = 0.1,mark.feature = FALSE)
identifyClusters(Sample2,n.dims.use = 30,n.neighbors = 30,knn.metric = "euclidean")
plot_umap_label_by_clusters(Sample2,show_method = "louvain",point.size = 0.80)

Sample2

runFA2_ForceDirectedGraph(Sample2,n.dims.use = 20,
                          n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)
plot_forceDirectedGraph_label_by_clusters(Sample2,show_method = "louvain",vertex.size = 0.85,
                                          background.color = "black")

mouse <- read.csv("mouse_genesets.csv",header=F)
write.table(mouse,file = "mouse_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

plot_umap_label_by_multiple_gene_sets(Sample2,gmt.file = "mouse_genesets.gmt",
                                      show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","orange","cyan","purple","green"),
                                      isNormalizedByHouseKeeping = T,point.size = 1,background.color = "black")

plot_forceDirectedGraph_label_by_multiple_gene_sets(Sample2,gmt.file = "mouse_genesets.gmt",
                                                    show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                                    custom_color = c("red","orange","cyan","purple","green"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1,edge.size = 0.1,
                                                    background.color = "black")
findMarkerGenes(Sample2,cluster.type = "louvain")
plot_heatmap_for_marker_genes(Sample2,cluster.type = "louvain",n.TopGenes = 10,rowFont.size = 5,use_raster = FALSE)
save(Sample2,file="Sample2.SingCellaR.rdata")
