library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

####读取差异表达基因，将基因ID从GENE_SYMBOL转换为ENTREZ_ID：####
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.table("mir21.csv",col.names="ENTREZID")

#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

gene <- info
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)

GO1<- as.data.frame(GO)
write.csv(GO1,file = "GO.csv",row.names = T)

KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)          

#使用David软件中富集的KEGG和GO作图，在微信收藏夹里
library(ggplot2)

#读取作图数据，这是一组上调基因的 GO 富集结果，展示了 top 富集的条目
kegg_enrich <- read.csv('KEGG.csv',header = TRUE)

#柱形图，纵坐标是 GO Term，横坐标是各 GO Term 的富集得分（Enrichment_score），颜色按 p 值着色
ggplot(kegg_enrich, aes(Term, Fold.Enrichment)) +
  geom_col(aes(fill = P.Value), width = 0.5) +
  scale_fill_gradient(high = 'blue', low = 'red') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Enrichment Score')





library(ggplot2)

#读取作图数据，这是一组上调基因的 GO 富集结果，展示了 top 富集的条目
go2<- read.csv('go2.csv', header = TRUE)
go3<-go2[go2$PValue<0.05,]
go_enrich<-go3[go3$FDR<0.05,]
#按富集得分由低到高排个序
go_enrich <- go_enrich[order(go_enrich$Fold.Enrichment, decreasing = FALSE), ]
go_enrich$Term <- factor(go_enrich$Term, levels = go_enrich$Term)

#气泡图，纵坐标是 GO Term，横坐标是各 GO Term 的富集得分（Enrichment_score）
#按各 GO Term 中富集的基因数量（Count）赋值气泡图中点的大小，颜色按 p 值着色
ggplot(go_enrich, aes(Term, Fold.Enrichment)) +
  geom_col(aes(fill = PValue), width = 0.5) +
  scale_fill_gradient(high = 'blue', low = 'red') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  coord_flip() +
  labs(x = '', y = 'Enrichment Score')

#或者根据 Category 绘制分面图
p + facet_grid(Category, scale = 'free_y', space = 'free_y')

library(ggplot2)

#读取作图数据，这是一组上调基因的 GO 富集结果，展示了 top 富集的条目

#按富集得分由低到高排个序
go_enrich <- go_enrich[order(go_enrich$Enrichment_score, decreasing = FALSE), ]
go_enrich$Term <- factor(go_enrich$Term, levels = go_enrich$Term)

#气泡图，纵坐标是 GO Term，横坐标是各 GO Term 的富集得分（Enrichment_score）
#按各 GO Term 中富集的基因数量（Count）赋值气泡图中点的大小，颜色按 p 值着色
p <- ggplot(go_enrich, aes(Term, Enrichment_score)) +
  geom_point(aes(size = Count, color = pValue)) +
  scale_size(range = c(2, 6)) +
  scale_color_gradient(high = 'blue', low = 'red') +
  theme(panel.background = element_rect(color = 'black', fill = 'transparent'), 
        panel.grid = element_blank(), legend.key = element_blank()) +
  coord_flip() +
  labs(x = '', y = 'Enrichment Score')
