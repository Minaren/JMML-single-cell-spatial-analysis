if(!require(multtest))install.packages("multtest")
if(!require(R.utils))install.packages("R.utils")
if(!require(Seurat))install.packages("Seurat")
if(!require(R.utils))install.packages("dplyr")
if(!require(mindr))install.packages("tidyverse")
if(!require(mindr))install.packages("mindr")
####1.准备工作####
#1.1读取数据
rm(list = ls())
counts1 <- Read10X(data.dir = "Kras") 
#load the PBMC datasets
#1.2创建Seurat对象
Kras<-CreateSeuratObject(counts = counts1, project = "Kras", min.cells = 3,min.features=200
)
Kras#查看数据
Krasdataframe<-as.data.frame(Kras[["RNA"]]@counts)#保存数据
write.table(Krasdataframe,'krascounts.txt',sep = '\t')
####2.QC####
#2.1细胞过滤即QC条件设置
Kras[["percent.mt"]]<-PercentageFeatureSet(Kras,pattern = "mt")#计算与线粒体基因对应的转录物百分比,鼠源的pattern需要换成mt
head(Kras@meta.data,50)#查看数据前5行
Kras<-subset(Kras,subset = nFeature_RNA>200&nFeature_RNA<6000&percent.mt<10)

#2.2可视化QC结果
#(1)小提琴图

VlnPlot(Kras,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)#nfeature_RNA是指每个barcode中检测到的基因的数目，nCount-RNA是指每个barcode中识别的RNA数量
#(2)散点图查看特征之间的关系
plot1<-FeatureScatter(Kras,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(Kras,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1,plot2))

####3.normalizing the data标准化数据####
ncol(as.data.frame(Kras[["RNA"]]@counts))
Kras<-NormalizeData(Kras,normalization.method = "LogNormalize",scale.factor = 2000)
#默认情况下，使用全局缩放规范化方法LogNormalize，该方法通过总表达式对每个单元格的特征表达式度量进行标准化，并将其乘以一个缩放因子(默认为10,000)，然后对结果进行log转换。标准化值存储在pbmc[["RNA"]]@data中。

####4.Identification of highly variable features寻找高差异基因####
Kras<-FindVariableFeatures(Kras,selection.method = "vst",nfeatures = 2000)
#for PCA DoHeatmap
top10<-head(VariableFeatures(Kras),10)
#identity the 10 most highly variable genes
top10

plot1<-VariableFeaturePlot(Kras)
plot2<-LabelPoints(plot=plot1,points = top10,repel=TRUE)
plot1+plot2


####5.缩放数据####
#5.1缩放的标准
all.genes<-rownames(Kras)
Kras<-ScaleData(Kras,features = all.genes)
#pbmc<-ScaleData(pbmc)##如果觉得scale所有的基因比较缓慢的话，可以只选择高变基因

#5.2删除不需要的数据(校正协变量，矫正线粒体基因比例的影响)
#pbmc<-SclaeData(pbmc,vars.to.regress="percent.mt")
#非常耗时

####6.PCA主成分分析####
#6.1主成分分析
Kras<-RunPCA(Kras, features = VariableFeatures(object=Kras))
print(Kras[["pca"]],dims=1:5,nfeatures=5)         
VizDimLoadings(Kras,dims = 1:2,reduction ="pca")   
#6.2PCA可视化
#(1)dimplot
DimPlot(Kras,reduction = "pca")
#(2)DimHeatmap
DimHeatmap(Kras,dims = 1:30,cells=500,balanced = TRUE)

####7.determine the dimensionality' of the dataset计算数据维度
Kras<-JackStraw(Kras,num.replicate = 100)
Kras<-ScoreJackStraw(Kras,dims = 1:20)
#(1)
ElbowPlot(Kras)
#(2)
JackStrawPlot(Kras,dims = 1:15)#慎用，耗时非常长，大数据集不建议

####8.cluster the cells细胞聚类
Kras<-FindNeighbors(Kras,dims = 1:30)
Kras<-FindClusters(Kras,resolution = 0.6)#resolution值越高越精细
head(Idents(Kras),3)

####9.UMAP/tSNE非线性降维####
#(1)一般认为，umap适合大型数据
Kras<-RunUMAP(Kras,dims = 1:30)
DimPlot(Kras,reduction = "umap")
#(2)tsne适合小型数据，因为计算起来比较慢，但分离效果优于umap
Kras<-RunTSNE(Kras,dims = 1:30)
DimPlot(Kras,reduction = "tsne")

####10.细胞注释####
#10.1寻找marker基因
#（1）方法1：找cluster5的marker
cluster5.markers<-FindMarkers(CB35,ident.1 = 5,ident.2 = c(0,3),min.pct = 0.25)
#（2）方法3：找某个cluster与其他cluster的差异基因（marker）
cluster5.markers <- FindMarkers(seurat.obj, 
                                ident.1 = 5, # 找cluster5的marker
                                ident.2 = c(0, 3), # 和cluster 0-3对比
                                min.pct = 0.25)
head(cluster5.markers, n = 5)
#（3）方法3：找所有的cluster的marker
CB35.markers<-FindAllMarkers(CB35,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
if(!require(dplyr))install.packages("dplyr")#这个包是为了%<%的使用
CB35.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
VlnPlot(CB35,features = c("NKG7","PF4"),slot = "counts",log = TRUE)#小提琴图可视化cluster marker
VlnPlot(CB35,features = c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A","LY2"))
FeaturePlot(CB35,features = c("MS4A1","GNLY","CD3E","CD14","FCER1A"))#降维图可视化cluster marker

top10<-CB35.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
DoHeatmap(CB35,features = top10$gene)+NoLegend()
#http://biocc.hrbmu.edu.cn/CellMarker/
#查看某一基因在该分群中的表达情况
VlnPlot(Kras, features = c("ITPR3"))
#10.2标记cluster
new.cluster.ids<-c("Naive CD4 T","CD14+ Mono","Memory CD4 T","8","CD8 T","FCGR3A+ Mono","NK","DC","Platelet")
names(new.cluster.ids)<-levels(pbmc)
pbmc<-RenameIdents(pbmc,new.cluster.ids)
DimPlot(pbmc,reduction="umap",label=TRUE,pt.size=0.5)+NoLegend()
sessionInfo()#检查学习所用的包

saveRDS(pbmc,'pbmc.rds')#保存文件
pbmc<-readRDS('pbmc.rds')
