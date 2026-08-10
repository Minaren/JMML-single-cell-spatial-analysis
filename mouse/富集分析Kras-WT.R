rm(list = ls())
####安装及加载R包####
library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(msigdbr)
library(org.Mm.eg.db)#人类则用library(org.Hs.eg.db)
####读取数据####
load("D:/aJMML/YX/annotation_after/sce_har_anno.RData")
####筛选差异基因####
## 在寻找差异基因之前，把默认的assay切换为RNA。
HSC<- subset(sce.harm, idents = "HSC")
DefaultAssay(HSC) <- 'RNA'
## 定义好你想要在哪一个分群基础上找差异表达基因
head(HSC@meta.data)
Idents(HSC) <- 'orig.ident' 
## 在不同cluster/或者celltype中找差异表达基因
only.pos = TRUE # 只找上调的差异表达基因
logfc.threshold = 0.25 # 差异基因的avg_log2FC必须要大于0.25
#all.clusters.markers <- FindAllMarkers(sce.harm, only.pos = T, logfc.threshold = logfc.threshold)
#head(all.clusters.markers)
markers <- FindMarkers(HSC, ident.1 = "Kras", ident.2 = "WT", only.pos = TRUE, logfc.threshold = 0.25)
head(markers)
markers = markers %>% rownames_to_column('gene') %>% filter(p_val_adj < 0.05)
head(markers)
####ID 转换####
OrgDb = "org.Mm.eg.db" # 根据物种来指定
gene_convert <- bitr(markers$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = OrgDb)
markers = markers%>%inner_join(gene_convert,by=c("gene"="SYMBOL"))

####GO富集分析####
ont = "BP" # posiible value: BP, CC, MF, all
go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = FALSE)
head(go.results)

go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="BP",OrgD = OrgDb, readable = TRUE)
head(go.results)
write.csv(go.results@result, file = "GO_enrichment_results_BP.csv", row.names = FALSE)
dotplot(go.results,label_format=10000)
ggsave("GO_dotplot_BP.png", width = 10, height = 8, dpi = 300)
ggsave("GO_dotplot_BP.pdf", width = 10, height = 8)

barplot(go.results, showCategory = 10,label_format=10000)
ggsave("GO_barplot_BP.png", width = 10, height = 8, dpi = 300)
ggsave("GO_barplot_BP.pdf", width = 10, height = 8)

emapplot(pairwise_termsim(go.results))
ggsave("GO_emapplot_BP.png", width = 10, height = 8, dpi = 300)
ggsave("GO_emaplot_BP.pdf", width = 10, height = 8)

# KEGG 富集
organism = "mmu"  # 小鼠的 KEGG 物种缩写
kegg.results <- enrichKEGG(markers$ENTREZID, organism = organism)
write.csv(kegg.results@result, file = "KEGG_enrichment_results_mmu.csv", row.names = FALSE)

# 2. 设置可读性
kegg.results <- setReadable(kegg.results, OrgDb = OrgDb, keyType = 'ENTREZID')

# 3. 绘制并保存 KEGG 图形
dotplot(kegg.results, label_format = 10000)
ggsave("KEGG_dotplot_mmu.png", width = 10, height = 8, dpi = 300)
ggsave("KEGG_dotplot_mmu.pdf", width = 10, height = 8)

barplot(kegg.results, showCategory = 10, label_format = 10000)
ggsave("KEGG_barplot_mmu.pdf", width = 10, height = 8)
ggsave("KEGG_barplot_mmu.pdf", width = 10, height = 8)

emapplot(pairwise_termsim(kegg.results))
ggsave("KEGG_emapplot_mmu.png", width = 10, height = 8, dpi = 300)
ggsave("KEGG_emapplot_mmu.pdf", width = 10, height = 8)
# 4. 绘制并保存 Pathview 图形
diff_genes_avg_logFC = markers$avg_log2FC
names(diff_genes_avg_logFC) = markers$ENTREZID
pathview(gene.data = diff_genes_avg_logFC, species = organism, pathway.id = kegg.results@result$ID[3])
ggsave("Pathview_pathway_mmu.png", width = 10, height = 8, dpi = 300)
ggsave("Pathview_pathway_mmu.pdf", width = 10, height = 8)

####GSEA富集分析####
min.pct = 0.01 ## 至少多少比例的细胞表达这个基因，过滤一些只在极少数细胞中有表达的基因
logfc.threshold = 0.01 ## 过滤掉在两组中几乎没有差异的基因
markers.for.gsea <- FindMarkers(HSC, ident.1 = "Kras", ident.2 = "WT", min.pct = min.pct, logfc.threshold=logfc.threshold)
# GSEA 要求输入的是一个排好序的列表
Markers_genelist <- markers.for.gsea$avg_log2FC
names(Markers_genelist)= rownames(markers.for.gsea)
head(Markers_genelist)
Markers_genelist <- sort(Markers_genelist, decreasing = T)
#m_df_categories <- msigdbr(species = 'Mus musculus')
#unique_categories <- unique(m_df_categories$gs_cat)
#print(unique_categories)
# 导入MSigDB
m_df <- msigdbr(species = 'Mus musculus', category = "C2")#人m_df = msigdbr(species = 'Homo sapiens' , category = "C2")
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea.results <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea.results)
gseaplot(gsea.results, gsea.results@result$ID[1])
# 保存富集结果
write_csv(gsea.results %>% data.frame, "gsea_results.csv")
head(Markers_genelist)
