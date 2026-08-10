####素有类型细胞比例####
library(Seurat) 
load("D:/aJMML/scRNA1/sc_seurat_integr.RData")

celltype = data.frame(ClusterID = 0:13,celltype = 'NA')
celltype[celltype$ClusterID %in% c(0),2] = "HSC"
celltype[celltype$ClusterID %in% c(10,1,5),2] = "MPP"
celltype[celltype$ClusterID %in% c(8,9),2] = 'CMP'
celltype[celltype$ClusterID %in% c(11),2] = 'CLP'
celltype[celltype$ClusterID %in% c(3),2] = 'MEP'
celltype[celltype$ClusterID %in% c(4,12),2] = 'Ery' 
celltype[celltype$ClusterID %in% c(13),2] = "B" 
celltype[celltype$ClusterID %in% c(2,6,7),2] = 'Granulocyte'


head(celltype)
celltype 
table(celltype$celltype)
sc_integr$celltype = "NA"
for(i in 1:nrow(celltype)){
  sc_integr@meta.data[which(sc_integr@meta.data$integrated_snn_res.0.6 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sc_integr@meta.data$celltype)

Idents(sc_integr) = sc_integr$celltype
Idents(sc_integr) = factor(Idents(sc_integr),levels = c('HSC','MPP','CLP','CMP','MEP','Ery','B','Granulocyte'))
sc_integr$celltype <- Idents(sc_integr)
Idents(sc_integr) <- factor(Idents(sc_integr),levels = rev(levels(Idents(sc_integr))))

HSC <- subset(sc_integr, celltype=="HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.4)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
DimPlot(HSC, label = T,pt.size = 1)

save(HSC,file = "sc_seurat_integr_HSC.RData")

DefaultAssay(HSC) <- "RNA"
