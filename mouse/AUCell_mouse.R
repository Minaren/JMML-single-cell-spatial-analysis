library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
#devtools::install_github("chuiqin/irGSEA")
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.5)
#HSC <- FindClusters(HSC, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)

mouse <- read.csv("LT_HSC2_mouse.csv",header=F)
write.table(mouse,file = "LT_HSC2_mouse.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
H <- read.gmt("LT_HSC2_mouse.gmt")
cells_rankings <- AUCell_buildRankings(HSC@assays$RNA@data,splitByBlocks=TRUE)# 关键一步
unique(H$term)
geneSets <- lapply(unique(H$term), function(x){H$gene[H$term == x]})
names(geneSets) <- unique(H$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
length(rownames(cells_AUC@assays@data$AUC))

#grep("INFLAMMATORY",rownames(cells_AUC@assays@data$AUC),value = T)#正则表达式
# [1] "HALLMARK_INFLAMMATORY_RESPONSE"
#geneSet <- "HALLMARK_INFLAMMATORY_RESPONSE"
geneSet<-"LT_HSC"
getAUC(cells_AUC)
#aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
aucs <- as.numeric(getAUC(cells_AUC))
HSC$AUC <- aucs
df<- data.frame(HSC@meta.data, HSC@reductions$umap@cell.embeddings)
head(df)
class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
ggplot(df, aes(UMAP_1, UMAP_2))+
  geom_point(aes(colour= AUC)) + 
  viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label =seurat_clusters),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + 
  theme_bw()

class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
ggplot(df, aes(UMAP_1, UMAP_2))+
  geom_point(aes(colour= AUC)) + 
  viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + 
  theme_bw()

ggplot(df, aes(UMAP_1, UMAP_2))+
  geom_point(aes(colour= AUC)) +
  scale_colour_gradientn(colours = c("#FFFFFF","red"))+
  ggrepel::geom_label_repel(aes(label = seurat_clusters),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + 
  theme_bw()

ggplot(df,aes(UMAP_1,UMAP_2,col= AUC))+
  geom_point(size=0.1)+
  scale_colour_gradientn(colours = c("gray40","#1E90FF","orange","red"))+
  ggrepel::geom_text_repel(data = df,
                           aes(median.1,median.2, label =labels),
                           size=4,
                           color= "gray20",
                           min.segment.length = 0)
 
ggsave("AUCell_LTHSC_mouse.pdf",height = 6,width = 6)


geneSet<-"LT_HSC_Rodriguez_et_al"
getAUC(cells_AUC)
#aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
aucs <- as.numeric(getAUC(cells_AUC))
HSC$AUC <- aucs
df<- data.frame(HSC@meta.data, HSC@reductions$umap@cell.embeddings)
head(df)
class_avg <- df %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
ggplot(df, aes(UMAP_1, UMAP_2))+
  geom_point(aes(colour= AUC)) + 
  viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label =seurat_clusters),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = "none") + 
  theme_bw()