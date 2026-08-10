rm(list=ls())
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)#人的用(org.Hs.eg.db)
allmarkers<-read.csv("D:/aJMML/scRNA/allmarkers.csv")
ids=bitr(allmarkers$gene,'SYMBOL','ENTREZID','org.Mm.eg.db') ## 将SYMBOL转成ENTREZID
allmarkers=merge(allmarkers,ids,by.x='gene',by.y='SYMBOL')
View(allmarkers)

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(allmarkers$ENTREZID, allmarkers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
## GO
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("godotplot.pdf", plot = p, width = 17, height = 42)
