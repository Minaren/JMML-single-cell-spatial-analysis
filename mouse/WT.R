if(!require(multtest))install.packages("multtest")
if(!require(R.utils))install.packages("R.utils")
if(!require(Seurat))install.packages("Seurat")
if(!require(R.utils))install.packages("dplyr")
if(!require(mindr))install.packages("tidyverse")
if(!require(mindr))install.packages("mindr")
library(patchwork)
####1.准备工作####
#1.1读取数据
rm(list = ls())
counts1 <- Read10X(data.dir = "WT") 
#load the PBMC datasets
#1.2创建Seurat对象
WT<-CreateSeuratObject(counts = counts1, project = "WT", min.cells = 3,min.features=200
)
WT#查看数据
WTdataframe<-as.data.frame(WT[["RNA"]]@counts)#保存数据
write.table(WTdataframe,'WTcounts.txt',sep = '\t')
####2.QC####
#2.1细胞过滤即QC条件设置

WT[["percent.mt"]]<-PercentageFeatureSet(WT,pattern = "mt")#计算与线粒体基因对应的转录物百分比,鼠源的pattern需要换成mt

head(WT@meta.data,50)#查看数据前5行
WT<-subset(WT,subset = nFeature_RNA>200&nFeature_RNA<6000&percent.mt<10)

#2.2可视化QC结果
#(1)小提琴图

VlnPlot(WT,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)#nfeature_RNA是指每个barcode中检测到的基因的数目，nCount-RNA是指每个barcode中识别的RNA数量
#(2)散点图查看特征之间的关系
plot1<-FeatureScatter(WT,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(WT,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1,plot2))

####3.normalizing the data标准化数据####
ncol(as.data.frame(WT[["RNA"]]@counts))
WT<-NormalizeData(WT,normalization.method = "LogNormalize",scale.factor = 2000)
#默认情况下，使用全局缩放规范化方法LogNormalize，该方法通过总表达式对每个单元格的特征表达式度量进行标准化，并将其乘以一个缩放因子(默认为10,000)，然后对结果进行log转换。标准化值存储在pbmc[["RNA"]]@data中。

####4.Identification of highly variable features寻找高差异基因####
WT<-FindVariableFeatures(WT,selection.method = "vst",nfeatures = 2000)
#for PCA DoHeatmap
top10<-head(VariableFeatures(WT),10)
#identity the 10 most highly variable genes
top10

plot1<-VariableFeaturePlot(WT)
plot2<-LabelPoints(plot=plot1,points = top10,repel=TRUE)
plot1+plot2


####5.缩放数据####
#5.1缩放的标准
all.genes<-rownames(WT)
WT<-ScaleData(WT,features = all.genes)
#pbmc<-ScaleData(pbmc)##如果觉得scale所有的基因比较缓慢的话，可以只选择高变基因

#5.2删除不需要的数据(校正协变量，矫正线粒体基因比例的影响)
#pbmc<-SclaeData(pbmc,vars.to.regress="percent.mt")
#非常耗时

####6.PCA主成分分析####
#6.1主成分分析
WT<-RunPCA(WT, features = VariableFeatures(object=WT))
print(WT[["pca"]],dims=1:5,nfeatures=5)         
VizDimLoadings(WT,dims = 1:2,reduction ="pca")   
#6.2PCA可视化
#(1)dimplot
DimPlot(WT,reduction = "pca")
#(2)DimHeatmap
DimHeatmap(WT,dims = 1:30,cells=500,balanced = TRUE)

####7.determine the dimensionality' of the dataset计算数据维度
WT<-JackStraw(WT,num.replicate = 100)
WT<-ScoreJackStraw(WT,dims = 1:20)
#(1)
ElbowPlot(WT)
#(2)
JackStrawPlot(WT,dims = 1:15)#慎用，耗时非常长，大数据集不建议

####8.cluster the cells细胞聚类
WT<-FindNeighbors(WT,dims = 1:30)
WT<-FindClusters(WT,resolution = 0.6)#resolution值越高越精细
head(Idents(WT),3)

####9.UMAP/tSNE非线性降维####
#(1)一般认为，umap适合大型数据
WT<-RunUMAP(WT,dims = 1:30)
DimPlot(WT,reduction = "umap")
#(2)tsne适合小型数据，因为计算起来比较慢，但分离效果优于umap
WT<-RunTSNE(WT,dims = 1:30)
DimPlot(WT,reduction = "tsne")
