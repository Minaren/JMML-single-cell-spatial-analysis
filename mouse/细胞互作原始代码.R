####安装和加载包####
# 安装必要的包
#devtools::install_github("satijalab/seurat-data")
devtools::install_github("sqjin/CellChat")
BiocManager::install("BiocNeighbors")
# 加载所需的 R 包
library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(future)

# 设置选项
options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 4)  # 支持并行处理


####数据准备####
# 加载 PBMC 数据
pbmc3k.final <- readRDS("pbmc.rds")

# 查看分析流程记录
pbmc3k.final@commands$FindClusters

# 使用 mindr 包构建思维导图（可选）
if (!requireNamespace("mindr", quietly = TRUE)) install.packages("mindr")
library(mindr)
(out <- capture.output(str(pbmc3k.final)))
out2 <- paste(out, collapse = "\n")
mm(gsub("\\.\\.@", "# ", gsub("\\.\\. ", "#", out2)), type = "text", root = "Seurat")

####3. 创建 CellChat 对象####
cellchat <- createCellChat(object = pbmc3k.final, group.by = "cell_type")
summary(cellchat)
levels(cellchat@idents)  # 查看细胞类型
groupSize <- as.numeric(table(cellchat@idents))  # 计算每种细胞群的数量

####4. 导入配体-受体数据库####
CellChatDB <- CellChatDB.human  # 对人类数据使用
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use  # 设置用于分析的数据库


####5. 数据预处理####
cellchat <- subsetData(cellchat)  # 筛选信号基因
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)  # 基于蛋白相互作用校正信号基因表达值


####推断细胞通讯网络####
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

# 推断信号通路层面的通讯
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

####7. 可视化细胞互作####
# 全局网络可视化
cellchat <- aggregateNet(cellchat)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weights/strength")
#逐一绘制每种细胞类型的通讯
mat <- cellchat@net$count
par(mfrow = c(3, 3), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                   arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

####8. 分析信号通路贡献####
# 查看显著的信号通路
pathways.show <- c("TGFb")
vertex.receiver <- c(1, 2, 4, 6)

# 不同布局的可视化
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# 提取对信号通路贡献最大的配体-受体对
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1, ]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver)

####9. 批量绘制所有信号通路####
pathways.show.all <- cellchat@netP$pathways
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"), layout = "circle")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot = gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}
setwd("../")
#泡泡图展示
p <- netVisual_bubble(cellchat, sources.use = c(3, 5, 7, 8, 9), 
                      targets.use = c(1, 2, 4, 6), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.pdf", p, width = 8, height = 12)



