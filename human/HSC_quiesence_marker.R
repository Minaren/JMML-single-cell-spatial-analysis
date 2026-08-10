genes_to_check =c("CD41","CD229","CD244","CD63","CD82")
genes_to_check =c("ITGA2B","LY9","CD244","CD63","CD82","ITGA6","PTEN","KLF4","HIF1A","FBW7","PML","NFATC1","ATM","GFI1",
                  "NECDIN","CDK6","CDKN1A","CDKN1B")
DefaultAssay(HSC) = "RNA"
DotPlot(HSC,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
FeaturePlot(object = HSC, features = c('PTEN',"KLF4","ITGA2B","LY9","CD244","CD63","CD82"),min.cutoff = "q10", max.cutoff = "q90")
