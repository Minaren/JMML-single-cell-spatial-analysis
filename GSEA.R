mouse <- read.csv("mouse_geneset_GSEA.csv",header=F)
write.table(mouse,file = "mouse_geneset_GSEA.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

library(Seurat)
library(clusterProfiler)
load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.7)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
library(clustree)
clustree(HSC)
Idents(HSC) <- "RNA_snn_res.0.7"
HSC$seurat_clusters <- HSC@active.ident
DimPlot(HSC, label = T,pt.size = 1)

tHSC1<-subset(HSC, idents = 1)
tHSC1 <- FindMarkers(tHSC1, min.pct = 0.25, 
                   logfc.threshold = 0.25,
                   group.by = "orig.ident",
                   ident.1 ="Kras",
                   ident.2="WT")
write.csv(tHSC1,'tHSC1_deg.csv')

tHSC2<-subset(HSC, idents = 2)
tHSC2 <- FindMarkers(tHSC2, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")
write.csv(tHSC2,'tHSC2_deg.csv')

tHSC0<-subset(HSC, idents = 0)
tHSC0 <- FindMarkers(tHSC0, min.pct = 0.25, 
                     logfc.threshold = 0.25,
                     group.by = "orig.ident",
                     ident.1 ="Kras",
                     ident.2="WT")
write.csv(tHSC0,'tHSC0_deg.csv')
saveRDS(HSC, "HSC_r0.7.rds")


alldiff <- tHSC1[order(tHSC1$avg_log2FC,decreasing = T),]
genelist <- alldiff$avg_log2FC
names(genelist) <- rownames(alldiff)
head(genelist)

hallmark <- read.gmt("mouse_geneset_GSEA.gmt")
gsea.re1<- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = hallmark,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH')  #指定 p 值校正方法
write.table(gsea.re1, 'gsea_tHSC1_custom.txt', sep = '\t', row.names = FALSE, quote = FALSE)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

library(ggsci)
library(enrichplot)
col_gsea1<-pal_simpsons()(16)

num1=1
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num1],
          title = "tHSC1",#标题
          color = col_gsea1[1:num1],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = TRUE,#p值表格
          ES_geom = "line"#line or dot
)

#展示多个基因集，设置多个名称即可，比如这里展示4个，如果有特定基因集，定义完成后传入geneSetID参数即可
num2=4
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)


#tHSC2
alldiff <- tHSC2[order(tHSC2$avg_log2FC,decreasing = T),]
genelist <- alldiff$avg_log2FC
names(genelist) <- rownames(alldiff)
head(genelist)
hallmark <- read.gmt("mouse_geneset_GSEA.gmt")
gsea.re1<- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = hallmark,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH')  #指定 p 值校正方法
write.table(gsea.re1, 'gsea_tHSC2_custom.txt', sep = '\t', row.names = FALSE, quote = FALSE)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

col_gsea1<-pal_simpsons()(16)


#展示多个基因集，设置多个名称即可，比如这里展示4个，如果有特定基因集，定义完成后传入geneSetID参数即可
num2=1
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "tHSC2",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = TRUE,#p值表格
          ES_geom = "line"#line or dot
)

alldiff <- tHSC0[order(tHSC0$avg_log2FC,decreasing = T),]
genelist <- alldiff$avg_log2FC
names(genelist) <- rownames(alldiff)
head(genelist)
hallmark <- read.gmt("mouse_geneset_GSEA.gmt")
gsea.re1<- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = hallmark,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH')  #指定 p 值校正方法
write.table(gsea.re1, 'gsea_tHSC0_custom.txt', sep = '\t', row.names = FALSE, quote = FALSE)
g1<-as.data.frame(gsea.re1)
g1<-subset(g1,p.adjust<0.05)
g1<-g1[order(g1$NES,decreasing = T),]

col_gsea1<-pal_simpsons()(16)


#展示多个基因集，设置多个名称即可，比如这里展示4个，如果有特定基因集，定义完成后传入geneSetID参数即可
num2=2
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "tHSC0",#标题
          color = col_gsea1[1:num2],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = TRUE,#p值表格
          ES_geom = "line"#line or dot
)


