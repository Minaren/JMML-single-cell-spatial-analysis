剪视频不容易  求关注三连啦

#xena.R 文件有更新，看我上传的最新版
#设置工作目录
setwd("xena")
#install.packages("tidyverse")
library(tidyverse)#每次重新打开R都要library一下
#xena官网
https://xenabrowser.net/datapages/

##文件的读取
#读取tsv文件
counts1 = read.table(file = 'TCGA-LIHC.htseq_counts.tsv', sep = '\t', header = TRUE) 
rownames(counts1) <- counts1[,1] #Alt <- 
x <- counts1[,1:3]
counts1 = counts1[,-1]
#substr函数
substr("wanglihong",1,4)
#table函数
table(substr(colnames(counts1),14,16))
#c("01A","11A")
#%in%符号用于判断是否属于
counts1 <- counts1[,substr(colnames(counts1),14,16)%in% c("01A","11A")]

table(substr(colnames(counts1),14,16))

#保留行名前15位
rownames(counts1) <- substr(rownames(counts1),1,15)
ceiling(1.2)
ceiling(3.8)
counts <- ceiling(2^(counts1)-1)

##文件的输出
#输出为文本
write.table(counts,"counts.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#输出为表格
write.csv(counts, file = "counts.csv")

#9.11
####差异分析####
#设置工作目录
setwd("xena")
##读取文本文件
counts <- read.table("counts.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#加载基因注释文件
Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#美元符号代表提取列
#取行名交集
comgene <- intersect(rownames(counts),rownames(Ginfo))
counts <- counts[comgene,]
class(counts)#判断数据类型
class(comgene)
Ginfo <- Ginfo[comgene,]
a <- rownames(counts)
b <- rownames(Ginfo)
identical(a,b)

counts$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
counts <- counts[!duplicated(counts$Gene),]   #去重复
rownames(counts) <- counts$Gene   #将行名变为Gene Symbol
ncol(Ginfo)
nrow
counts <- counts[,-ncol(counts)]   #去除最后一列
write.table(counts, file = "LIHC_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的counts
tumor <- colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A <- counts[,tumor]
write.table(counts_01A, file = "LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#差异分析
library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)

counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"),levels = c("N","T"))) %>% 
  column_to_rownames("sample")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = conditions,
  design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)
save(res,file = "LIHC_DEG.rda")#一定要保存！
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需要
#DEG:differentially expressed genes
#读取差异基因文件
#直接在文件夹双击

####整理fpkm文件####
#与counts几乎相同，fpkm不需进行log转换
#读取tsv文件
library(tidyverse)
setwd("xena")
fpkm1 = read.table(file = 'TCGA-LIHC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1]  
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1

Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] #只要编码RNA
#取行名交集
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
fpkm <- fpkm[comgene,]
Ginfo <- Ginfo[comgene,]
fpkm$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
fpkm <- fpkm[!duplicated(fpkm$Gene),]   #去重复
rownames(fpkm) <- fpkm$Gene   #将行名变为Gene Symbol
fpkm <- fpkm[,-ncol(fpkm)]   #去除最后一列
#保存所以患者的fpkm文件
write.table(fpkm, file = "LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的fpkm文件
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#

####9.14####
setwd("xena")
library(tidyverse)
fpkm_01A <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#读取之前的差异分析结果，还记得怎么读取吗#
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)

gene <- c("LIN28B","CTAG2","REG3A")
a <- fpkm_01A[gene,]
b <- fpkm_11A[gene,]
a <- t(a)
b <- t(b)
class(a)
a <- as.data.frame(a)
b <- as.data.frame(b)
##运用传导符%>%  cltrl+shift+M 
a <- a %>% t() %>% as.data.frame()
b <- b %>% t() %>% as.data.frame()
write.csv(a, file = "01A.csv")
write.csv(b, file = "11A.csv")
#Graphpad等等


####GEO数据库的使用####
####代表什么
#pubmed插件 文献检索
#GEO网站：https://www.ncbi.nlm.nih.gov/geo/
####GSE84402####
setwd("GSE84402")
###加载R包
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
chooseBioCmirror()
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]

#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))
##2.2 通过exprs函数获取表达矩阵并校正
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###数据校正
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
#使用R包转换id
index = gset[[1]]@annotation
if(!require("hgu133plus2.db"))
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(29:30)]
write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####中秋快乐####
setwd("GSE84402")
###加载R包
library(tidyverse)
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))

##读取上节课整理好的表达数据exp##
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#差异分析
library(limma)
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

##热图##
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




####cox回归分析####
#设置工作目录
setwd("cox")
#安装加载R包
install.packages("survival")
install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#表达数据整理完毕
#读取tcga差异分析结果
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)
#整合
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)

#Cox分析
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}


write.table(Coxoutput, file = "cox results.txt",sep = "\t",row.names = F,col.names = T,quote = F)
###筛选top基因
pcutoff <- 0.001
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] # 取出p值小于阈值的基因
topgene <- topgene[1:10,]

#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "12" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()


####计算患者免疫评分与肿瘤纯度#####
setwd("TCGA ESTIMATE")  #设置工作目录
#安装包
library(utils) #这个包应该不用下载，自带的
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "LIHC_fpkm_mRNA_01A.txt",   #输入文件名
                  output.f = "LIHC_fpkm_mRNA_01A.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("LIHC_fpkm_mRNA_01A.gct",   #刚才的输出文件名
              "LIHC_fpkm_mRNA_01A_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打分
result <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分



####ROC####
#读取生存信息tsv文件
setwd("ROC")
library(tidyverse)
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#保存整理好的生存信息
write.table(surv, file = "survival.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#提取上次cox作图的10个基因
gene <- c("SPP1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()
#整合表达谱与生存信息
exp_sur <- cbind(exp10,surv)
#准备R包
install.packages("ROCR")
install.packages("rms")
library(ROCR)
library(rms)
#构建ROC预测模型
ROC1 <- prediction(exp_sur$SPP1,exp_sur$OS)   #构建ROC预测模型 
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR值
AUC <- performance(ROC1,"auc")   #计算曲线下面积(AUC)值

AUC<- 0.5604839 #改 根据结果对AUC进行赋值

#1.4 绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)   #绘制对角线
dev.off()
