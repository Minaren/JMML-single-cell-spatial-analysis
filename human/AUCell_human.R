library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
devtools::install_github("chuiqin/irGSEA")
load("D:/aJMML/scRNA1/sc_seurat_integr.RData")
mouse <- read.csv("human_HSC_signature.csv",header=F)
write.table(mouse,file = "human_HSC_signature.gmt",sep = "\t",row.names = F,col.names = F,quote = F)
H <- read.gmt("human_HSC_signature.gmt")
cells_rankings <- AUCell_buildRankings(sc_integr@assays$RNA@data,splitByBlocks=TRUE)# 关键一步
unique(H$term)
geneSets <- lapply(unique(H$term), function(x){H$gene[H$term == x]})
names(geneSets) <- unique(H$term)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
length(rownames(cells_AUC@assays@data$AUC))

#grep("INFLAMMATORY",rownames(cells_AUC@assays@data$AUC),value = T)#正则表达式
# [1] "HALLMARK_INFLAMMATORY_RESPONSE"
#geneSet <- "HALLMARK_INFLAMMATORY_RESPONSE"

geneSet<-"HSC_signature"
getAUC(cells_AUC)
#aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
aucs <- as.numeric(getAUC(cells_AUC))
sc_integr$AUC <- aucs
df<- data.frame(sc_integr@meta.data, sc_integr@reductions$umap@cell.embeddings)
head(df)
class_avg <- df %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
ggplot(df, aes(UMAP_1, UMAP_2))+
  geom_point(aes(colour= AUC)) + 
  viridis::scale_color_viridis(option="A") +
  ggrepel::geom_label_repel(aes(label = cell_type),
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
  scale_colour_gradientn(colours = c("#FFFFFF","orange","red"))+
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
 
ggsave("AUCell_UMAP.pdf",height = 6,width = 6)
