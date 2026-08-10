library(clustree)
library(patchwork)
load("D:/aJMML/scRNA/sc_filt.RData")
sc_n= NormalizeData(sc_filt)
sc_n= FindVariableFeatures(sc_n,selection.method = "vst",nfeatures = 2000, verbose = FALSE)
sc_n = ScaleData(sc_n)
sc_n

sc_n <- RunPCA(sc_n, npcs=50, verbose=FALSE)
ElbowPlot(sc_n, ndims = 50)
sce_ha = sc_n
pc.num=1:30
sce_har<-RunHarmony(sce_ha,group.by.vars = "orig.ident",project.dim = F,plot_convergence = T)

sce.har = FindNeighbors(sce.har, dims = pc.num,reduction = "harmony")
sce.har = FindClusters(sce.har, resolution = c(seq(0,1.6,.2)))
table(sce.har@active.ident)

####clustree观察分辨率####
p1 <- clustree(sce.har@meta.data, prefix = 'RNA_snn_res.') + coord_flip()
#p2 <- DimPlot(sce.harm, group.by = 'RNA_snn_res.0.6', label = T)
#p<-p1 + p2 + plot_layout(widths = c(3, 1))
ggsave("RNA_snn_res.png", plot = p1, width = 17, height = 17)
ggsave("RNA_snn_res.pdf", plot = p1, width = 17, height = 17)

####AUC####
##res = 0.6
sce.har %>% 
  SetIdent(value = 'RNA_snn_res.0.6') %>% 
  FindAllMarkers(test.use = 'roc') %>% 
  filter(myAUC > 0.6) %>% 
  count(cluster, name = 'number')
##res = 0.2
sce.har %>% 
  SetIdent(value = 'RNA_snn_res.0.2') %>% 
  FindAllMarkers(test.use = 'roc') %>% 
  filter(myAUC > 0.6) %>% 
  count(cluster, name = 'number')
