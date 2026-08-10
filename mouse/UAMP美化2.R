#载入相关R包：
rm(list = ls())
install.packages('tidydr')
library(tidydr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cols4all)

#读取scRNA-seq测试数据：
##已完成分析的Seurat对象
load("D:/aJMML/scRNA1/sce_har_anno.RData")
sce<-sce.harm
sce[['celltype']] <- Idents(sce)
dim(sce)
celltype <- Idents(sce)
table(celltype)
#UMAP提取方法相同：
UMAP <- as.data.frame(sce@reductions$umap@cell.embeddings)
UMAP <- cbind(UMAP,celltype)
head(UMAP)
#建立自定义主题：
mytheme <- theme_void() + #空白主题，便于我们后期添加UMAP箭头
  theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整

#建立映射，添加散点：
p <- ggplot(data = UMAP, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = celltype),
             size = 0.4,
             alpha = 0.8)
p

#添加实线椭圆置信区间：
p1 <- p +
  stat_ellipse(aes(color = celltype),
               level = 0.95, linetype = 1, show.legend = F) +
  mytheme
p1

#添加虚线椭圆置信区间：
p2 <- p +
  stat_ellipse(aes(color = celltype),
               level = 0.95, linetype = 2, show.legend = F) +
  mytheme
p2


#添加填充型置信区间：
p3 <- p +
  stat_ellipse(aes(color = celltype, fill = celltype),
               level = 0.95, linetype = 1, show.legend = F,
               geom = 'polygon', alpha = 0.1) +
  mytheme
p3

#添加UMAP坐标轴箭头：
p4 <- p3 +
  theme_dr(xlength = 0.2, #x轴长度
           ylength = 0.2, #y轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), #箭头大小/长度
                               ends = 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())
p4

#在图中增加亚群标签:
##计算每个亚群散点的中位数，作为图中标签的坐标
label <- UMAP %>%
  group_by(celltype)%>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
head(label)

p5 <- p4 +
  geom_text(data = label,
            aes(x = UMAP_1, y = UMAP_2, label = celltype),
            fontface = "bold", #粗体强调
            color = 'black', size = 4)
p5

p6 <- p5 +
  guides(color = guide_legend(override.aes = list(size = 8))) #放大图例中的散点
p6


#自定义颜色：
#c4a_gui() #查看本次使用配色包的色板
mycol <-c(
  "#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da",
  "#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84")
p7 <- p6 +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol)
p7




