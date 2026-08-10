load("D:/aJMML/scRNA1/sce_har_anno.RData")
HSC <- subset(sce.harm, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC,resolution =0.7)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
DimPlot(HSC, label = T,pt.size = 1)

g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(HSC)) #提取HSC矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(HSC)) #提取HSC矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对HSC进行细胞周期评分
HSC <- CellCycleScoring(HSC, g2m.features=g2m_genes, s.features=s_genes)
colnames(HSC@meta.data)
table(HSC$Phase)
table(HSC$orig.ident)#查看各组细胞数
prop.table(table(Idents(HSC)))
table(Idents(HSC), HSC$orig.ident)#各组不同细胞群细胞数
table(Idents(HSC), HSC$Phase)
Cellratio <- prop.table(table(HSC$Phase,Idents(HSC)), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_cellcycle.pdf',width = 7,height = 6,dpi = 500)

#MPP亚群再分析
MPP <- subset(sce.harm, celltype=="MPP")
MPP <- ScaleData(MPP, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
MPP <- FindVariableFeatures(MPP, nfeatures = 4000)
MPP <- RunPCA(MPP, npcs = 50, verbose = FALSE)
MPP <- FindNeighbors(MPP, reduction = "pca", dims = 1:50)
MPP <- FindClusters(MPP, resolution =0.7 )
MPP <- RunUMAP(MPP, reduction = "pca", dims = 1:50)
MPP$seurat_clusters <- MPP@active.ident
DimPlot(MPP, label = T,pt.size = 1)
MPP <- CellCycleScoring(MPP, g2m.features=g2m_genes, s.features=s_genes)
table(MPP$orig.ident)#查看各组细胞数
prop.table(table(Idents(MPP)))
table(Idents(MPP), MPP$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(MPP$Phase,Idents(MPP)), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('percentage_MPP_cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('percentage_MPP_cellcycle.pdf',width = 7,height = 6,dpi = 500)

#tHSC1
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tHSC1)) #提取HSC矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tHSC1)) #提取HSC矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对HSC进行细胞周期评分
tHSC1 <- CellCycleScoring(tHSC1, g2m.features=g2m_genes, s.features=s_genes)
colnames(tHSC1@meta.data)
table(tHSC1$Phase)
table(tHSC1$orig.ident)#查看各组细胞数
prop.table(table(Idents(tHSC1)))
table(Idents(tHSC1), tHSC1$orig.ident)#各组不同细胞群细胞数
table(tHSC1$orig.ident, tHSC1$Phase)
Cellratio <- prop.table(table(tHSC1$Phase,tHSC1$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tHSC1',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tHSC1cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tHSC1cellcycle.pdf',width = 7,height = 6,dpi = 500)

#tHSC2
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tHSC2)) #提取HSC矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tHSC2)) #提取HSC矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对HSC进行细胞周期评分
tHSC2 <- CellCycleScoring(tHSC2, g2m.features=g2m_genes, s.features=s_genes)
colnames(tHSC2@meta.data)
table(tHSC2$Phase)
table(tHSC2$orig.ident)#查看各组细胞数
prop.table(table(Idents(tHSC2)))
table(Idents(tHSC2), tHSC2$orig.ident)#各组不同细胞群细胞数
table(tHSC2$orig.ident, tHSC2$Phase)
Cellratio <- prop.table(table(tHSC2$Phase,tHSC2$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tHSC2',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tHSC2cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tHSC2cellcycle.pdf',width = 7,height = 6,dpi = 500)


g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tHSC0)) #提取HSC矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tHSC0)) #提取HSC矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对HSC进行细胞周期评分
tHSC0 <- CellCycleScoring(tHSC0, g2m.features=g2m_genes, s.features=s_genes)
colnames(tHSC0@meta.data)
table(tHSC0$Phase)
table(tHSC0$orig.ident)#查看各组细胞数
prop.table(table(Idents(tHSC0)))
table(Idents(tHSC0), tHSC0$orig.ident)#各组不同细胞群细胞数
table(tHSC0$orig.ident, tHSC0$Phase)
Cellratio <- prop.table(table(tHSC0$Phase,tHSC0$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tHSC0',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tHSC0cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tHSC0cellcycle.pdf',width = 7,height = 6,dpi = 500)

tMPP0<-subset(MPP, idents = 0)
tMPP1<-subset(MPP, idents = 1)
tMPP2<-subset(MPP, idents = 2)
tMPP3<-subset(MPP, idents = 3)
tMPP4<-subset(MPP, idents = 4)
#tMPP1
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tMPP1)) #提取MPP矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tMPP1)) #提取MPP矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对MPP进行细胞周期评分
tMPP1 <- CellCycleScoring(tMPP1, g2m.features=g2m_genes, s.features=s_genes)
colnames(tMPP1@meta.data)
table(tMPP1$Phase)
table(tMPP1$orig.ident)#查看各组细胞数
prop.table(table(Idents(tMPP1)))
table(Idents(tMPP1), tMPP1$orig.ident)#各组不同细胞群细胞数
table(tMPP1$orig.ident, tMPP1$Phase)
Cellratio <- prop.table(table(tMPP1$Phase,tMPP1$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tMPP1',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tMPP1cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tMPP1cellcycle.pdf',width = 7,height = 6,dpi = 500)

#tMPP2
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tMPP2)) #提取MPP矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tMPP2)) #提取MPP矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对MPP进行细胞周期评分
tMPP2 <- CellCycleScoring(tMPP2, g2m.features=g2m_genes, s.features=s_genes)
colnames(tMPP2@meta.data)
table(tMPP2$Phase)
table(tMPP2$orig.ident)#查看各组细胞数
prop.table(table(Idents(tMPP2)))
table(Idents(tMPP2), tMPP2$orig.ident)#各组不同细胞群细胞数
table(tMPP2$orig.ident, tMPP2$Phase)
Cellratio <- prop.table(table(tMPP2$Phase,tMPP2$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tMPP2',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tMPP2cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tMPP2cellcycle.pdf',width = 7,height = 6,dpi = 500)


g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tMPP0)) #提取MPP矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tMPP0)) #提取MPP矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对MPP进行细胞周期评分
tMPP0 <- CellCycleScoring(tMPP0, g2m.features=g2m_genes, s.features=s_genes)
colnames(tMPP0@meta.data)
table(tMPP0$Phase)
table(tMPP0$orig.ident)#查看各组细胞数
prop.table(table(Idents(tMPP0)))
table(Idents(tMPP0), tMPP0$orig.ident)#各组不同细胞群细胞数
table(tMPP0$orig.ident, tMPP0$Phase)
Cellratio <- prop.table(table(tMPP0$Phase,tMPP0$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tMPP0',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tMPP0cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tMPP0cellcycle.pdf',width = 7,height = 6,dpi = 500)


#tMPP3
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tMPP3)) #提取MPP矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tMPP3)) #提取MPP矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对MPP进行细胞周期评分
tMPP3 <- CellCycleScoring(tMPP3, g2m.features=g2m_genes, s.features=s_genes)
colnames(tMPP3@meta.data)
table(tMPP3$Phase)
table(tMPP3$orig.ident)#查看各组细胞数
prop.table(table(Idents(tMPP3)))
table(Idents(tMPP3), tMPP3$orig.ident)#各组不同细胞群细胞数
table(tMPP3$orig.ident, tMPP3$Phase)
Cellratio <- prop.table(table(tMPP3$Phase,tMPP3$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tMPP3',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tMPP3cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tMPP3cellcycle.pdf',width = 7,height = 6,dpi = 500)

#tMPP4
g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(tMPP4)) #提取MPP矩阵中的G2M期marker基因
s_genes <- cc.genes$s.genes   #获取S期marker基因 
s_genes <- CaseMatch(search=s_genes, match=rownames(tMPP4)) #提取MPP矩阵中的S期marker基因 
#通过提取到的g2m期基因和s期基因，使用CellCycleScoring函数，对MPP进行细胞周期评分
tMPP4 <- CellCycleScoring(tMPP4, g2m.features=g2m_genes, s.features=s_genes)
colnames(tMPP4@meta.data)
table(tMPP4$Phase)
table(tMPP4$orig.ident)#查看各组细胞数
prop.table(table(Idents(tMPP4)))
table(Idents(tMPP4), tMPP4$orig.ident)#各组不同细胞群细胞数
table(tMPP4$orig.ident, tMPP4$Phase)
Cellratio <- prop.table(table(tMPP4$Phase,tMPP4$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.2,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='tMPP4',y = 'Ratio')+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('tMPP4cellcycle.png',width = 7,height = 6,dpi = 500)
ggsave('tMPP4cellcycle.pdf',width = 7,height = 6,dpi = 500)

