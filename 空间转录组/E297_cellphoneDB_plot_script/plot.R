library("ktplots")
library("CellChat")
library(tidyr)
setwd("/data/project/E297/WT/")
df.net <- read.table("count_network.txt",header = T,sep = "\t",stringsAsFactors = F)
df.net <- spread(df.net, TARGET, count)
rownames(df.net) <- df.net$SOURCE
df.net <- df.net[, -1]
df.net <- as.matrix(df.net)

pvals_stat <- read.delim("pvalues.txt", check.names = FALSE)
mymeans <- read.delim("means.txt", check.names = FALSE)
#pdf("CellPhoneDB.heatmap.pdf")
#p<-plot_cpdb_heatmap(pvals = pvals_stat, cellheight = 10, cellwidth = 10,return_tables =T)
#dev.off()
#write.table(p$count_network,"count_network.txt",row.names=T,sep="\t",quote=F)

#df.net <-p$count_network
#meta.data <- read.table("/data/project/20240902_Chenglong_4SingleCell/8Samples/rename/CellphoneDB/data/LN2_meta.txt",
#                        header = T,sep = "\t",stringsAsFactors = F)
#groupSize <- as.numeric(table(meta.data$cell_type))

pdf("CellPhoneDB.net.circle.All.pdf")
netVisual_circle(df.net,
#                 vertex.weight = groupSize,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

pdf(paste0("CellPhoneDB.net.circle.split.pdf"),12,10)
par(mar=c(1,1,1,1), mfrow=c(4,4))
for (i in 1:nrow(df.net)) {
  mat2 <- matrix(0, nrow = nrow(df.net), ncol = ncol(df.net), dimnames = dimnames(df.net))
  mat2[i, ] <- df.net[i, ]
  netVisual_circle(mat2,
#                   vertex.weight = groupSize,
                   weight.scale = T,
                   edge.weight.max = max(df.net),
                   title.name = "",
                   arrow.size=0.5)
}

dev.off()

source("/data/project/E297/cellphoneDB_Dotplot.mod.R")
pdf("cellphoneDB_Dotplot.HSC.pdf")
cellphoneDB_Dotplot(pvals.data = pvals_stat,means.data = mymeans,key="HSC",
                    target.cells_1 = c("HSC"),target.cells_2=c("HSC"))
dev.off()
