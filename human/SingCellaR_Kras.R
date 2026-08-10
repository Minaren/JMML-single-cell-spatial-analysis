#if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager")
#BiocManager::install("AUCell")

#devtools::install_github('supatt-lab/SingCellaR')
#install.packages("AUCell")

rm(list=ls())
library(SingCellaR)
library(AUCell)
data_matrices_dir<-"D:/aJMML/JMML/Kras"
Sample1<-new("SingCellaR")
Sample1@dir_path_10x_matrix<-data_matrices_dir
Sample1@sample_uniq_id<-"Sample1"

load_matrices_from_cellranger(Sample1,cellranger.version = 3)

Sample1

process_cells_annotation(Sample1,mito_genes_start_with ="mt-")

####QC####
plot_cells_annotation(Sample1,type="histogram")
plot_cells_annotation(Sample1,type="boxplot")
plot_UMIs_vs_Detected_genes(Sample1)
DoubletDetection_with_scrublet(Sample1)
table(Sample1@meta.data$Scrublet_type)

filter_cells_and_genes(Sample1,
                       min_UMIs=1000,
                       max_UMIs=60000,
                       min_detected_genes=500,
                       max_detected_genes=6000,
                       max_percent_mito=15,
                       genes_with_expressing_cells = 10,isRemovedDoublets = F)
normalize_UMIs(Sample1,use.scaled.factor = T)
remove_unwanted_confounders(Sample1,residualModelFormulaStr="~UMI_count+percent_mito")
get_variable_genes_by_fitting_GLM_model(Sample1,mean_expr_cutoff = 0.05,disp_zscore_cutoff = 0.05)
plot_variable_genes(Sample1)

runPCA(Sample1,use.components=50,use.regressout.data = T)
plot_PCA_Elbowplot(Sample1)
runUMAP(Sample1,dim_reduction_method = "pca",n.dims.use = 20,n.neighbors = 30,
        uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(Sample1,feature = "UMI_count",point.size = 0.1,mark.feature = FALSE)
identifyClusters(Sample1,n.dims.use = 30,n.neighbors = 30,knn.metric = "euclidean")
plot_umap_label_by_clusters(Sample1,show_method = "louvain",point.size = 0.80)

Sample1
class(Sample1)
runFA2_ForceDirectedGraph(Sample1,n.dims.use = 20,
                          n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)
plot_forceDirectedGraph_label_by_clusters(Sample1,show_method = "louvain",vertex.size = 0.85,
                                          background.color = "black")

mouse <- read.csv("mouse_genesets.csv",header=F)
write.table(mouse,file = "mouse_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
