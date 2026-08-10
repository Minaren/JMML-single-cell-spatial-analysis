####1.读取目标genneset文件####
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets<-human
genesets <- read.csv("AUCell_human.csv",header=F)
genesets <- subset(genesets, select = c("V1","V2")) %>% as.data.frame()
genesets <- split(genesets$V2, genesets$V1)

####2.提取分组平均表达矩阵####
#Idents设置按什么来取组的表达均值（计算case control之间的均值也可以）
expr <- AverageExpression(HSC, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

####3.GSVA富集分析####
# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesets, method="ssgsea") 
saveRDS(gsva.res, "gsva.res.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res.csv", row.names = F)
pheatmap::pheatmap(gsva.res, show_colnames = T, scale = "row")
