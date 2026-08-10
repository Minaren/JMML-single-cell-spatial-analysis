sce.harm <- readRDS('sce.harmrenamed.rds')
sce.harm
DimPlot(sce.harm)
names(sce.harm@meta.data)
unique(sce.harm$group)
DimPlot(sce.harm,split.by = 'group')
sce.harm$celltype.group <- paste(sce.harm$celltype, sce.harm$group, sep = "_")
sce.harm$celltype <- Idents(sce.harm)
Idents(sce.harm) <- "celltype.group"
mydeg <- FindMarkers(sce.harm,ident.1 = 'VSMC_AS1',ident.2 = 'VSMC_C57', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
cellfordeg<-levels(sce.harm$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(sce.harm, ident.1 = paste0(cellfordeg[i],"_P3"), ident.2 = paste0(cellfordeg[i],"_AS1"), verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()
library(dplyr)
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10
sce.harm <- ScaleData(sce.harm, features =  rownames(sce.harm))
DoHeatmap(sce.harm,features = top10,size=3)
Idents(sce.harm) <- "celltype"
VlnPlot(sce.harm,features = top10,split.by = 'group',idents = 'EC')
FeaturePlot(sce.harm,features = top10,split.by = 'group')
#DotPlot(sce.harm,features = top10,split.by ='group')#默认只有两种颜色
DotPlot(sce.harm,features = top10,split.by ='group',cols = c('blue','yellow','pink'))
####提取表达量#######
mymatrix <- as.data.frame(sce.harm@assays$RNA@data)
mymatrix2<-t(mymatrix)%>%as.data.frame()
mymatrix2[,1]<-sce.harm$celltype
colnames(mymatrix2)[1] <- "celltype"

mymatrix2[,ncol(mymatrix2)+1]<-sce.harm$group
colnames(mymatrix2)[ncol(mymatrix2)] <- "group"

#绘图
library(ggplot2)
p1<- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=Thbs1,fill=group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p1
