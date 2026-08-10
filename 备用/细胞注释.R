if(T){rm(list = ls())
if (!require("Seurat"))install.packages("Seurat")
if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager")
if (!require("multtest", quietly = TRUE))install.packages("multtest")
if (!require("dplyr", quietly = TRUE))install.packages("dplyr")
download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','pbmc3k_filtered_gene_bc_matrices.tar.gz')
library(R.utils)
gunzip('pbmc3k_filtered_gene_bc_matrices.tar.gz')
untar('pbmc3k_filtered_gene_bc_matrices.tar')
pbmc.data <- Read10X('filtered_gene_bc_matrices/hg19/') 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)      
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10) 
pbmc <- ScaleData(pbmc, features =  rownames(pbmc)) 
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) 
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) 
pbmc <- JackStraw(pbmc, num.replicate = 100) 
pbmc <- ScoreJackStraw(pbmc, dims = 1:20) 
JackStrawPlot(pbmc, dims = 1:15)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) 
head(Idents(pbmc), 5) 
pbmc <- RunUMAP(pbmc, dims = 1:10) 
pbmc <- RunTSNE(pbmc, dims = 1:10) 
DimPlot(pbmc, reduction = "umap", label = TRUE)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
library(dplyr)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
}


####方法一： 查数据库
library(dplyr)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()#展示前10个标记基因的热图
## Warning in DoHeatmap(pbmc, features = top10$gene): The following features were
## omitted as they were not found in the scale.data slot for the RNA assay: CD8A,
## VPREB3, CD40LG, PIK3IP1, PRKCQ-AS1, NOSIP, LEF1, CD3E, CD3D, CCR7, LDHB, RPS3A
VlnPlot(pbmc, features = top10$gene[1:20])#可以看到小提琴图每一个表达值
VlnPlot(pbmc, features = top10$gene[1:20],pt.size=0)
#可以看出CD14是非常好的marker,而PRS3A则非常糟糕
DimPlot(pbmc,label = T)
#通过标记基因及文献，可以人工确定各分类群的细胞类型，则可以如下手动添加细胞群名称
bfreaname.pbmc <- pbmc
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
#帮助单细胞测序进行注释的数据库：
#http://mp.weixin.qq.com/s?__biz=MzI5MTcwNjA4NQ==&mid=2247502903&idx=2&sn=fd21e6e111f57a4a2b6c987e391068fd&chksm=ec0e09bddb7980abf038f62d03d3beea6249753c8fba69b69f399de9854fc208ca863ca5bc23&mpshare=1&scene=24&srcid=1110SJhxDL8hmNB5BThrgOS9&sharer_sharetime=1604979334616&sharer_shareid=853c5fb0f1636baa0a65973e8b5db684#rd
#cellmarker: http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

####方法二：通过singleR进行注释,这一方法很鸡肋，主要是由于很难找到合适的参考数据，但是对于免疫细胞的注释还是有可观效果的
if(!require(SingleR))BiocManager::install(SingleR)
if(!require(matrixStats))BiocManager::install('matrixStats')
if(!require(celldex))BiocManager::install('celldex')#有一些版本冲突，需要重新安装一些包。
assay_for_SingleR <- GetAssayData(bfreaname.pbmc, slot="data")#取出样本中的表达序列
predictions <- SingleR(test=assay_for_SingleR, 
                       ref=cg, labels=cg$label.main)
#以kidney中提取的阵列为输入数据，以小鼠的阵列作为参考，predict细胞类型
table(predictions$labels)#看看都注释到了哪些细胞
cellType=data.frame(seurat=bfreaname.pbmc@meta.data$seurat_clusters,
                    predict=predictions$labels)#得到seurat中编号与预测标签之间的关系
sort(table(cellType[,1]))
table(cellType[,1:2])#访问celltyple的2~3列

#当你有合适的参考数据集时可以用方法三方法四进行注释
####方法三：自定义singleR的注释
library(SingleR)
library(Seurat)
library(ggplot2)
## Warning: package 'ggplot2' was built under R version 4.0.5
if(!require(textshape))install.packages('textshape')
if(!require(scater))BiocManager::install('scater')
if(!require(SingleCellExperiment))BiocManager::install('SingleCellExperiment')
library(dplyr)
myref <- pbmc##这里为了检测，我们将参考数据集与目标数据集用同一个数据进行测试
myref$celltype <- Idents(myref)
table(Idents(myref))

Refassay <- log1p(AverageExpression(myref, verbose = FALSE)$RNA)#求
#Ref <- textshape::column_to_rownames(Ref, loc = 1)#另一种得到参考矩阵的办法
head(Refassay)#看看表达矩阵长啥样
ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=Refassay))
#参考数据集需要构建成一个SingleCellExperiment对象
ref_sce=scater::logNormCounts(ref_sce)
logcounts(ref_sce)[1:4,1:4]
colData(ref_sce)$Type=colnames(Refassay)
ref_sce#构建完成

#提取自己的单细胞矩阵
testdata <- GetAssayData(bfreaname.pbmc, slot="data")
pred <- SingleR(test=testdata, ref=ref_sce, 
                labels=ref_sce$Type,
                #clusters = scRNA@active.ident
)
table(pred$labels)
head(pred) 
as.data.frame(table(pred$labels))
cellType=data.frame(seurat=bfreaname.pbmc@meta.data$seurat_clusters,
                    predict=pred$labels)#得到seurat中编号与预测标签之间的关系
sort(table(cellType[,1]))
table(cellType[,1:2])#访问celltyple的2~3列
lalala <- as.data.frame(table(cellType[,1:2]))
finalmap <- lalala %>% group_by(seurat) %>% top_n(n = 1, wt = Freq)#找出每种seurat_cluster注释比例最高的对应类型
finalmap <-finalmap[order(finalmap$seurat),]$predict#找到seurat中0：8的对应预测细胞类型
print(finalmap)

testname <- bfreaname.pbmc
new.cluster.ids <- as.character(finalmap)
names(new.cluster.ids) <- levels(testname)
testname <- RenameIdents(testname, new.cluster.ids)

p1 <- DimPlot(pbmc,label = T)
p2 <- DimPlot(testname,label = T)#比较一下测试数据与参考数据集之间有没有偏差
p1|p2#完美，无差别注释，当然了，我们这个参考数据用的比较极端


#####方法4：利用seurat内置的原先用于细胞整合的功能，将参考数据与待注释数据进行映射处理
library(Seurat)
pancreas.query <- bfreaname.pbmc#待注释数据
pancreas.anchors <- FindTransferAnchors(reference = pbmc, query = pancreas.query,
                                        dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
#把注释加回原来的数据集
pancreas.query$prediction.match <- pancreas.query$predicted.id 
table(pancreas.query$prediction.match)
Idents(pancreas.query)<- 'prediction.match'

p1 <- DimPlot(pbmc,label = T) 
p2 <- DimPlot(pancreas.query,label = T)#比较一下测试数据与参考数据集之间有没有偏差
p1|p2
