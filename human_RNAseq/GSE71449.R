library("limma")
library("tidyverse")
library("ggplot2")
library("dplyr")
library("ggVolcano")
library("magrittr")
library("ggrepel")
library("pheatmap")
library("GEOquery")
if(!require(ggsci))install.packages("ggsci")
if(!require(ggpubr))install.packages("ggpubr")

####整理表达数据####FF
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE9476', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
#提取子集
gset[[1]]

###pdata获得疾病表型数据
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
table(pdata$title)

###exprs函数获取表达数据
exp <- exprs(gset[[1]])
#exp<-as.data.frame(exp)

write.table(exp,file = "exp.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp,file = "exp.csv")
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)


###limma包进行数据校正
library(limma)
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)#一般绝对值在20以内的就是经过转化的，所以不需要后面两个步骤
#exp <- log2(exp+1)
#range(exp)
write.csv(exp,file = "exp_校正后.csv")

####id注释####
GPL=getGEO(filename = 'GSE9476_family.soft.gz')
gpl=GPL@gpls[[1]]@dataTable@table
colnames(gpl)
ids=gpl[,c(1,11)]
write.table(ids,file = "ids.txt",sep = "\t",row.names=F,col.names = T)
write.csv(ids,file = "ids.csv")

ids<-read.csv(file = "ids.csv")
colnames(ids) = c("X","probe_id","symbol")#注意！csv里面有一列是X，所以增加了一列X
exp=as.data.frame(exp)
exp$probe_id=rownames(exp)  #将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id")  # 合并数据
# 按照symbol列去重
exp2=exp2[!duplicated(exp2$symbol),]
# 数据框probe_exp的行名变成symbol
rownames(exp2)=exp2$symbol
exp2=exp2[,c(-1,-ncol(exp2))]
exp2<-exp2[,-66]#这里也是因为前面增加了一列X所以要去掉
#输出ids转化后的表达量文件
write.table(exp2,file = "ids_exprs.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp2,file = "ids_exprs.csv")
write.csv(pdata,file = "pdata.csv")


#观察某个基因在实验组和对照组中表达量
#这里是把基因表达量和分组信息整合到excel里了

library(ggpubr)

# Read data
ITPR3 <- read.csv(file = "ITPR3.csv")
colnames(ITPR3) <- c("group", "gene", "expression")

# Create a ggboxplot
P <- ggboxplot(ITPR3, x = "group", y = "expression", color = "group", 
               ylab = "Gene expression",
               xlab = "",
               legend.title = "Your Legend Title",  # Replace with your actual legend title
               palette = c("#66c2a5", "#fc8d62"),  # Custom colors
               width = 0.6, add = "none")

# Rotate x-axis labels for better readability
p <- P + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add statistical significance annotations centered in the middle
p1 <- p + stat_compare_means(aes(group = group),
                             method = "wilcox.test",
                             symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                symbols = c("", "*", "*", " ")),
                             label = "p.signif", 
                             size = 4, vjust = 0.5)  # Adjust size and vjust for better placement

# Customize the plot theme
p2 <- p1 + theme_minimal(base_size = 14) +
  theme(legend.position = "top", legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold"))

# Save the plot
ggsave("p1_updated.pdf", p2, width = 10, height = 8)







ITPR3<-read.csv(file = "ITPR3.csv")
colnames(ITPR3) = c("group","gene","expression")
head(ITPR3)
ITPR3<-as.data.frame(ITPR3)
P<-ggboxplot(ITPR3,x="group", y="expression", color = "group", 
            ylab="Gene expression",
            xlab="",
            legend.title="x",
            palette = c("blue","red"),
            width=0.6, add = "none")

p=P+rotate_x_text(60)
p1<-p+stat_compare_means(aes(group=group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

ggsave("p1.pdf",width =10, height =10)

ANGPTL2<-read.csv(file = "ANGPTL2.csv")
colnames(ANGPTL2) = c("group","gene","expression")
head(ANGPTL2)
ANGPTL2<-as.data.frame(ANGPTL2)
P<-ggboxplot(ANGPTL2,x="group", y="expression", color = "group", 
             ylab="Gene expression",
             xlab="",
             legend.title=x,
             palette = c("blue","red"),
             width=0.6, add = "none")

p=P+rotate_x_text(60)
p1<-p+stat_compare_means(aes(group=group),
                         method="wilcox.test",
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                         label = "p.signif")

ggsave("ANGPTL2.pdf",width =10, height =10)




PEIZO1<-read.csv(file = "PIEZO1.csv")
colnames(PEIZO1) = c("group","gene","expression")
head(PEIZO1)
PEIZO1<-as.data.frame(PEIZO1)
P<-ggboxplot(PEIZO1,x="gene", y="expression", color = "group", 
             ylab="Gene expression",
             xlab="",
             legend.title=x,
             palette = c("blue","red"),
             width=0.6, add = "none")

p=P+rotate_x_text(60)
p2<-p+stat_compare_means(aes(group=group),
                         method="wilcox.test",
                         symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                         label = "p.signif")

ggsave("p2.pdf",width =10, height =10)


#输出图片
pdf(width=6, height=5)
print(p1)
dev.off()


####差异分析####
#读取整理好的表达数据
exp <- read.table("ids_exprs.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#读取样本信息表
group<-read.csv("sampleinfo.csv",header = T,row.names = 1,sep = ",")
#差异分析
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(exp)

contrast.matrix <- makeContrasts(Control-JMML,levels = design)

fit <- lmFit(exp,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
colnames(DEG)

##标记上下调基因
DEG$regulate <- ifelse(DEG$P.Value > 0.05, "unchanged",
                       ifelse(DEG$logFC > 1, "up-regulated",
                              ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))

foldChange = 1 # 自定义修改筛选参数
padj = 0.05 # 自定义修改筛选参数

All_diffSig <- DEG[(DEG$adj.P.Val < padj & (DEG$logFC > foldChange | DEG$logFC < (-foldChange))),]
#dim(All_diffSig)
write.csv(All_diffSig, "all_diffsig_filtered.csv") 


diffup <-  All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC > foldChange)),]
write.csv(diffup, "diffup_filtered.csv")
diffdown <- All_diffSig[(All_diffSig$P.Value < padj & (All_diffSig$logFC < -foldChange)),]
write.csv(diffdown,"diffdown_filtered.csv")

DEG1<-DEG
DEG1$Genes<-rownames(DEG1)

pdf("volcano1.pdf",width = 15,height = 15)
ggvolcano(data = DEG1,x = "logFC",y = "P.Value",output = FALSE,label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") #标签位置为up right
dev.off()




DEG_genes <- DEG1[DEG1$adj.P.Val<0.05&abs(DEG1$logFC)>1,]
DEG_gene_expr <- exp[DEG1$Genes,]

expr_data1<-exp
expr_data1$Genes<-rownames(expr_data1)
DEG_gene_expr<-inner_join(DEG1, expr_data1, by = "Genes")

rownames(DEG_gene_expr)<-DEG_gene_expr$Genes
DEG_gene_expr<-DEG_gene_expr[,-1]#循环几次直到只剩表达量
DEG_gene_expr <- DEG_gene_expr[apply(DEG_gene_expr, 1, function(x) sd(x)!=0),]



pheatmap(DEG_gene_expr,
            color = colorRampPalette(c("blue","white","red"))(100), #颜色
            scale = "row", #归一化的方式
            border_color = NA, #线的颜色
            fontsize = 10, #文字大小
            show_rownames = F,
            filename ="DEG_pheatmap2.pdf")

p2<-pheatmap(DEG_gene_expr,
             color = colorRampPalette(c("blue","white","red"))(100), #颜色
             scale = "row", #归一化的方式
             border_color = NA, #线的颜色
             fontsize = 10, #文字大小
             show_rownames = T,
             filename ="DEG_pheatmap2.pdf")



























####热图####
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()

####火山图(完整版)#####

head(deg)
deg$'-log10(pvalue)' <- -log10(deg$P.Value)
library(ggpubr)
library(ggthemes)
ggscatter(deg, x="logFC", y="-log10(pvalue)")+theme_base()

#为了美观，新加一列Group
deg$Group= "not-significant"
#将adj.P.Val小于0.05，logFC大于2的基因设为显著上调基因
#将adj.P.Val小于0. 05，logFC小于2的基因设为显著下调基因
deg$Group[which( (deg$adj.P.Val<0.01)&(deg$logFC > 1.5))]= "up-regulated"
deg$Group[which( (deg$adj.P.Val<0.01)&(deg$logFC < -1.5))] = "down-regulated"
#查看上调和下调基因数目
table(deg$Group)

deg$Symbol<-rownames(deg)
#新加一列Label
deg$Label=""
#对差异表达基因的p值进行从小到大排序
deg<-deg[order(deg$adj.P.Val),]
#高表达的基因中，选择adj.P .Val最小的10个
up.genes<-head(deg$Symbol[which(deg$Group == "up-regulated")],10)
#低表达的基因中，选择adj.P.Val最小的10个
down.genes<- head(deg$Symbol[which(deg$Group == "down-regulated")],10)
#将up. genes和down. genes合并， 并加入到Label中
deg.top10.genes<-c(as.character(up.genes),as.character(down.genes))
deg$Label[match(deg.top10.genes, deg$Symbol)] <- deg.top10.genes

#绘制新的火山图
ggscatter(deg, x="logFC", y="-log10(pvalue)",
          color = "Group",
          palette=c("#2f5688","#BBBBBB","#CC0000"),
          size=1,
          label=deg$Label,
          font.label = 8,
          repel=T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",)+theme_base()+
  geom_hline(yintercept = 1.30,linetype="dashed")+
  geom_vline(xintercept = c(-2,2),linetype="dashed")

####GEO三大富集分析####
setwd("JMML-Control")
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library("clusterProfiler")
deg <- read.table("deg_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
deg <- deg %>% filter(change!="stable")

DEG <- deg
DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO分析
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
barplot(ego,showCategory=20,drop=T)
dotplot(ego,showCategory=50)

#KEGG
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
dotplot(kk, showCategory=30)

####GSEA
###安装并导入所需要的R包
BiocManager::install("clusterProfiler") #感谢Y叔的clusterprofiler包
BiocManager::install("enrichplot")  #画图需要
BiocManager::install("org.Hs.eg.db") #基因注释需要
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
###导入数据
#df = read.table("gene_diff.txt",header = T) #读入txt
df = read.csv("gene_diff.csv",header = T) #读入csv
head(df)#查看前面几行
dim(df)#数据总共几行几列

###转换基因id
#基因名是symbol，需要将基因ID转换为Entrez ID格式。Entrez ID实际上是指的Entrez gene ID，是对应于染色体上一个gene location的。每一个发现的基因都会被编制一个统一的编号，而Entrez ID是指的来自于NCBI旗下的Entrez gene数据库所使用的编号。因为Entrez ID具有特异性，所以后续分析更适合用Entrez ID。
df_id<-bitr(df$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db

#把两个数据框df 和 df_id根据SYMBOL列合并
df_all_gsea<-merge(df,df_id,by="SYMBOL",all=F)#使用merge合并
head(df_all_gsea) #再看看数据
dim(df_all_gsea) #因为有一部分没转换成功，所以数量就少了。
head(df_all_gsea)

###GAEA
df_all_sort <- df_all_gsea[order(df_all_gsea$logFC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$logFC #把foldchange按照从大到小提取出来
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(gene_fc)
KEGG <- gseKEGG(gene_fc, organism = "hsa") #具体参数在下面
head(KEGG)#看一下这个文件
#ID 代表KEGG中的信号通路:ID 代表KEGG中的信号通路;Description 对信号通路的描述;setSize该信号通路的基因个数;enrichmentScore 富集分数，也就是ES;NES 标准化以后的ES，全称normalized enrichment score、;qvalues ，或者说FDR q-val（false discovery rate）错误发现率;rank 排名;core_enrichment，富集该目的通路的基因列表。

#GO富集
GO <- gseGO(
  gene_fc, #gene_fc
  ont = "BP",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,#人类注释基因
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",#p值校正方法
)
#KEGG富集
gseKEGG(
  geneList,
  organism = "hsa",
  keyType = "kegg",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea",
  ...
)
sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]#按照enrichment score从高到低排序
head(sortKEGG)
dim(sortKEGG)
write.table(sortKEGG,"gsea_sortKEGG.txt") #保存结果

#gseaplot2用法
gseaplot2(
  x, #gseaResult object，即GSEA结果
  geneSetID,#富集的ID编号
  title = "", #标题
  color = "green",#GSEA线条颜色
  base_size = 11,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = FALSE, #是否添加 pvalue table
  ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
)
#基础操作
paths <- c("hsa03010", "hsa05152", "hsa05171", "hsa04512")#选取你需要展示的通路ID
gseaplot2(KEGG,paths, pvalue_table = TRUE)

gseaplot2(KEGG, "hsa05134", color = "firebrick", rel_heights=c(1, .2, .6))

#pathsGSEA
gseaplot2(KEGG,paths,color = colorspace::rainbow_hcl(4),subplots=c(1,2), pvalue_table = TRUE)
#换个颜色，只显示上面两个副图

