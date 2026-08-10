library(SingCellaR)
human_HSPCs <- new("SingCellaR")
human_HSPCs@dir_path_SingCellR_object_files<-"D:/aJMML/scRNA1/SingCellaR_objects/"
human_HSPCs@SingCellR_object_files=c("Sample1.SingCellaR.rdata",
                                     "Sample2.SingCellaR.rdata")

preprocess_integration(human_HSPCs)

human_HSPCs

filter_cells_and_genes(human_HSPCs,
                       min_UMIs=1000,
                       max_UMIs=80000,
                       min_detected_genes=500,
                       max_detected_genes=8000,
                       max_percent_mito=15,
                       genes_with_expressing_cells = 10,isRemovedDoublets = FALSE)

cell_anno.info<-get_cells_annotation(human_HSPCs)
head(cell_anno.info)
table(cell_anno.info$sampleID)

human_HSPCs@meta.data$status[human_HSPCs@meta.data$sampleID=="1_Sample1"]<-"Kras"
human_HSPCs@meta.data$status[human_HSPCs@meta.data$sampleID=="1_Sample2"]<-"WT"

cell_anno.info<-get_cells_annotation(human_HSPCs)
head(cell_anno.info)

normalize_UMIs(human_HSPCs,use.scaled.factor = T)

get_variable_genes_by_fitting_GLM_model(human_HSPCs,mean_expr_cutoff = 0.1,disp_zscore_cutoff = 0.1)
plot_variable_genes(human_HSPCs)

runPCA(human_HSPCs,use.components=50,use.regressout.data = F)
plot_PCA_Elbowplot(human_HSPCs)

SingCellaR::runUMAP(human_HSPCs,dim_reduction_method = "pca",n.dims.use = 20,n.neighbors = 30,
                    uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "sampleID",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "status",point.size = 0.5)


library(harmony)
runSupervised_Harmony(human_HSPCs,covariates = c("sampleID"),n.dims.use = 20,
                      hcl.height.cutoff = 0.25,harmony.max.iter = 20,n.seed = 1)
SingCellaR::runUMAP(human_HSPCs,useIntegrativeEmbeddings = T, integrative_method = "supervised_harmony",n.dims.use = 20,
                    n.neighbors = 30,uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "sampleID",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "status",point.size = 0.5)

mouse <- read.csv("mouse_genesets.csv",header=F)
write.table(mouse,file = "mouse_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
plot_umap_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                      show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","orange","cyan","purple","green"),
                                      isNormalizedByHouseKeeping = T,point.size = 1,background.color = "black")
runFA2_ForceDirectedGraph(human_HSPCs,n.dims.use = 20, useIntegrativeEmbeddings = T,integrative_method = "supervised_harmony",n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)
plot_forceDirectedGraph_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                                    show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                                    custom_color = c("red","orange","cyan","purple","green"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1,edge.size = 0.05,
                                                    background.color = "black")
runHarmony(human_HSPCs,covariates = c("sampleID"),n.dims.use = 20,harmony.max.iter = 20,n.seed = 1)
SingCellaR::runUMAP(human_HSPCs,useIntegrativeEmbeddings = T, integrative_method = "harmony",n.dims.use = 20,
                    n.neighbors = 30,uwot.metric = "euclidean")
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "sampleID",point.size = 0.5)
plot_umap_label_by_a_feature_of_interest(human_HSPCs,feature = "status",point.size = 0.5)
plot_umap_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                      show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","orange","cyan","purple","green"),
                                      isNormalizedByHouseKeeping = T,point.size = 1,background.color = "black")
runFA2_ForceDirectedGraph(human_HSPCs,n.dims.use = 20, useIntegrativeEmbeddings = T,integrative_method = "harmony",n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)
plot_forceDirectedGraph_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                                    show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                                    custom_color = c("red","orange","cyan","purple","green"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1,edge.size = 0.1,
                                                    background.color = "black")

save(human_HSPCs,file="SingCellaR_intergrtion.rdata")
