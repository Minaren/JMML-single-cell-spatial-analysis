
library(SingCellaR)
human_HSPCs<-local(get(load(file="SingCellaR_intergrtion.rdata")))
human_HSPCs
identifyClusters(human_HSPCs,useIntegrativeEmbeddings = T,integrative_method = "harmony", n.dims.use = 20,n.neighbors = 50, knn.metric = "euclidean")
plot_umap_label_by_clusters(human_HSPCs,show_method = "louvain")
findMarkerGenes(human_HSPCs,cluster.type = "louvain")
plot_heatmap_for_marker_genes(human_HSPCs,cluster.type = "louvain",n.TopGenes = 5,rowFont.size = 5)

mouse <- read.csv("mouse_genesets.csv",header=F)
write.table(mouse,file = "mouse_genesets.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
plot_umap_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                      show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","orange","cyan","purple","green"),
                                      isNormalizedByHouseKeeping = T,point.size = 1,background.color = "black")


runFA2_ForceDirectedGraph(human_HSPCs,n.dims.use = 20, useIntegrativeEmbeddings = T,integrative_method = "harmony",n.neighbors = 5,n.seed = 1,fa2_n_iter = 1000)
plot_forceDirectedGraph_label_by_clusters(human_HSPCs,show_method = "louvain")
plot_forceDirectedGraph_label_by_multiple_gene_sets(human_HSPCs,gmt.file = "mouse_genesets.gmt",
                                                    show_gene_sets = c("Erythroid","Lymphoid","Myeloid","Megakaryocyte","HSPC"),
                                                    custom_color = c("red","orange","cyan","purple","green"),
                                                    isNormalizedByHouseKeeping = T,vertex.size = 1,edge.size = 0.05,
                                                    background.color = "black")


genesets<-get_gmtGeneSets("mouse_genesets.gmt")
plot_3D_knn_graph_label_by_a_signature_gene_set(human_HSPCs,gene_list = genesets[["HSPC"]],vertext.size = 0.40)

pre_rankedGenes_for_GSEA<-identifyGSEAPrerankedGenes_for_all_clusters(human_HSPCs,
                                                                      cluster.type = "louvain")
save(pre_rankedGenes_for_GSEA,file="human_HSPCs_preRankedGenes_for_GSEA.rdata")


load("human_HSPCs_preRankedGenes_for_GSEA.rdata")
mouse_celltype <- read.csv("mouse_celltype.csv",header=F)
write.table(mouse_celltype,file = "mouse_celltype.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
fgsea_Results<-Run_fGSEA_for_multiple_comparisons(GSEAPrerankedGenes_list = pre_rankedGenes_for_GSEA,nPermSimple=2000,
                                                  gmt.file = "mouse_celltype.gmt")

plot_heatmap_for_fGSEA_all_clusters(fgsea_Results,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row = 5,
                                    adjusted_pval = 0.10,
                                    show_only_NES_positive_score = T,format.digits = 3,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)
pdf(filename = "heatmap_for_GSEA.pdf",width = 480,height = 480,units = "px",bg = "white",res = 72) 



plot_heatmap_for_fGSEA_all_clusters(fgsea_Results,isApplyCutoff = TRUE,
                                    use_pvalues_for_clustering=T,
                                    show_NES_score = T,fontsize_row = 4,
                                    adjusted_pval = 0.05,
                                    show_only_NES_positive_score = T,format.digits = 2,
                                    clustering_method = "ward.D",
                                    clustering_distance_rows = "euclidean",
                                    clustering_distance_cols = "euclidean",show_text_for_ns = F)

