#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
load("D:/aJMML/scRNA1/sce_har_anno.RData")
WT <- subset(sce.harm, orig.ident == "WT")
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(WT@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- WT@meta.data 
p_data$celltype <- WT@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(WT),row.names = row.names(WT))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds1<-cds
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。

#expressed_genes <- VariableFeatures(WT) 
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype",cores=1) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
##差异基因的结果文件保存
write.table(deg,file="train.monocle_DEG_WT.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
deg<-read.table("train.monocle_DEG_WT.xls",sep="\t")
deg<-deg[-1,]
ordergene <- deg$V5
#ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
table(cds@featureData@data[["use_for_ordering"]])
pdf("train_ordergenes_WT.pdf")
plot_ordering_genes(cds)
dev.off()
trace('project2MST', edit = T, where = asNamespace("monocle"))
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
#cds <- orderCells(cds)
#cds <- orderCells(cds,root_state = 4)
#以pseudotime值上色 
pdf("train.monocle.pseudotime_WT.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="Pseudotime", size=1,show_backbone=TRUE) 
dev.off()
#以细胞类型上色
pdf("train.monocle.celltype_WT.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by="celltype", size=1,show_backbone=TRUE)
dev.off()
#以细胞状态上色
pdf("train.monocle.state_WT.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "orig.ident",size=0.05,show_backbone=TRUE)
dev.off()

#按照seurat分群排序细胞
pdf("seurat.clusters_WT.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
dev.off()

#以细胞状态上色（拆分）“分面”轨迹图，以便更容易地查看每个状态的位置。
pdf("train_monocle_celltype_WT.pdf",width =30 ,height = 10)
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~celltype", nrow = 1)
dev.off()

library(ggsci)
p1=plot_cell_trajectory(cds, color_by = "celltype")  + scale_color_futurama() 
p2=plot_cell_trajectory(cds, color_by = "orig.ident")  + scale_color_futurama()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
p3=plot_cell_trajectory(cds, color_by = "orig.ident")  + scale_color_manual(values = colour)
p=p1|p2|p3
ggsave("p换色.pdf", plot = p, width = 30, height = 10)

p1 <- plot_cell_trajectory(cds, x = 1, y = 2, color_by = "celltype") + 
  theme(legend.position='none',panel.border = element_blank()) + #去掉第一个的legend
  scale_color_futurama()
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "celltype")+
  scale_color_futurama()+
  theme(legend.title = element_blank()) 
p=p1|p2
ggsave("p树形.pdf", plot = p2, width = 10, height = 10)
library(ggpubr)
df <- pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
View(df)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggplot(df, aes(Pseudotime, colour = orig.ident, fill=orig.ident)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
ggsave("pseudotime_orig.ident.pdf",width = 15, height = 10)

pdf("train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
save(cds,file = 'cds_WT.RData')
write.csv(pData(cds), "pseudotime_WT.csv")
