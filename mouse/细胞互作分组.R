# 加载必要的包
library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)

# 清理环境
rm(list = ls())
options(stringsAsFactors = FALSE)

# ----------------------------
# 数据准备部分
# ----------------------------
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")  # 设置Python环境

# 加载数据
load("D:/aJMML/scRNA1/sce_har_anno.RData")
table(sce.harm$celltype)  # 查看细胞类型分布
Idents(sce.harm) <- 'celltype'

# 筛选指定的细胞类型
#celltypes_to_include <- c('B cells', 'CD4+ T cells', 'CD8+ T cells', 'Dendritic cells', 'Monocytes', 'NK cells')
#scRNA <- subset(scRNA, idents = celltypes_to_include)
#scRNA$celltype <- as.factor(as.character(scRNA$celltype))


# 根据样本标识筛选
table(sce.harm$orig.ident)
Idents(scRNA) <- 'orig.ident'
scRNA1 <- subset(sce.harm, subset = orig.ident %in% c('WT'))
scRNA1 <- subset(sce.harm, subset = orig.ident %in% c('Kras'))
scRNA1 <- subset(sce.harm, subset = orig.ident %in% c('Remission'))
scRNA1 <- subset(sce.harm, subset = orig.ident %in% c('Resistance'))


# 创建CellChat对象
c_scRNA1 <- createCellChat(scRNA1@assays$RNA@data, meta = scRNA1@meta.data, group.by = "celltype")
c_scRNA2 <- createCellChat(scRNA2@assays$RNA@data, meta = scRNA2@meta.data, group.by = "celltype")
save(c_scRNA1, c_scRNA2, file = "cco.rda")

# ----------------------------
# 分析部分
# ----------------------------
# 设置分析输出路径
dir.create("./Compare", showWarnings = FALSE)
setwd("./Compare")

# 分析WT样本的细胞通讯网络
cellchat <- c_scRNA1
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ce_scRNA1 <- cellchat
saveRDS(ce_scRNA1, "ce_scRNA1.rds")

# 分析Kras样本的细胞通讯网络
cellchat <- c_scRNA2
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ce_scRNA2 <- cellchat
saveRDS(ce_scRNA2, "ce_scRNA2.rds")

# 合并CellChat对象
cco.list <- list(WT = ce_scRNA1, Kras= ce_scRNA2)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

# ----------------------------
# 可视化部分
# ----------------------------
# 可视化通讯数量和强度对比
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

# 差异网络图
par(mfrow = c(1, 2))
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
# save as Counts_Compare_net.pdf
# 保存为PDF

# 热图可视化
par(mfrow = c(1, 1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1 + h2
# save as Counts_Compare_net_heatmap.pdf
# 保存为PDF

# 指定细胞互作数量对比网络图
par(mfrow = c(1, 2))
s.cell <- c("HSC", "T")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = TRUE, label.edge = TRUE, edge.weight.max = weight.max, 
                 title.name = paste0("Number of interactions - ", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = TRUE, label.edge = TRUE, edge.weight.max = weight.max, 
                 title.name = paste0("Number of interactions - ", names(cco.list)[2]))
# save as Counts_Compare_select.pdf 10*6.5

## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 6)

#特定信号通路的对比
pathways.show <- c("IL16") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

#热图
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

####3.7 配体-受体对比分析
#气泡图展示所有配体受体对的差异

levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1,11), targets.use = c(1,11),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 12, height = 8)

#气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use = c(1,11), targets.use = c(1,11), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in TIL", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(1,11), targets.use = c(1,11), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", pc, width = 12, height = 5.5)

# 保存分析结果
saveRDS(cellchat, "cellchat.rds")


