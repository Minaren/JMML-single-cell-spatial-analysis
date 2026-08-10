Formal class 'Seurat' [package "SeuratObject"] with 13 slots
  #  assays      :List of 1
  #..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
  ####  counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ######  i       : int [1:44400707] 3 9 11 12 34 55 65 87 92 96 ...
  ######  p       : int [1:14889] 0 2257 5513 10449 13865 16689 18046 22462 25008 28687 ...
  ######  Dim     : int [1:2] 22495 14888
  ######  Dimnames:List of 2
  ######..$ : chr [1:22495] "Gm22307" "Gm6085" "Gm6119" "Mrpl15" ...
  ######..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ######  x       : num [1:44400707] 1 1 3 1 1 1 9 1 2 1 ...
  ######  factors : list()
  ####  data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ######  i       : int [1:44400707] 3 9 11 12 34 55 65 87 92 96 ...
  ######  p       : int [1:14889] 0 2257 5513 10449 13865 16689 18046 22462 25008 28687 ...
  ######  Dim     : int [1:2] 22495 14888
  ######  Dimnames:List of 2
  ######..$ : chr [1:22495] "Gm22307" "Gm6085" "Gm6119" "Mrpl15" ...
  ######..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ######  x       : num [1:44400707] 1.02 1.02 1.85 1.02 1.02 ...
  ######  factors : list()
  ####  scale.data   : num [1:2000, 1:14888] -0.0136 -0.0389 -0.0315 6.4792 -0.8658 ...
  ####..- attr(*, "dimnames")=List of 2
  #####..$ : chr [1:2000] "St18" "Prex2" "Sulf1" "Eya1" ...
  #####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ####  assay.orig   : NULL
  ####  var.features : chr [1:2000] "Hbb-bs" "Hba-a2" "Hba-a1" "Hbb-bt" ...
  ####  meta.features:'data.frame':	22495 obs. of  5 variables:
  ####..$ vst.mean                 : num [1:22495] 0.000537 0.002955 0.00047 1.460707 0.004433 ...
  ####..$ vst.variance             : num [1:22495] 0.000537 0.002947 0.000604 3.839525 0.004548 ...
  ####..$ vst.variance.expected    : num [1:22495] 0.00057 0.003477 0.000496 3.8714 0.005306 ...
  ####..$ vst.variance.standardized: num [1:22495] 0.943 0.848 1.218 0.992 0.857 ...
  ####..$ vst.variable             : logi [1:22495] FALSE FALSE FALSE FALSE FALSE FALSE ...
  ####  misc         : list()
  ####  key          : chr "rna_"
  #  meta.data   :'data.frame':	14888 obs. of  7 variables:
  #..$ orig.ident     : chr [1:14888] "Kras" "Kras" "Kras" "Kras" ...
  #..$ nCount_RNA     : num [1:14888] 5619 22243 28865 16048 16664 ...
  #..$ nFeature_RNA   : int [1:14888] 2257 3256 4936 3416 2824 1357 4416 2546 3679 2251 ...
  #..$ percent.mt     : num [1:14888] 2.011 0.0854 0.4989 0.4798 0.18 ...
  #..$ RNA_snn_res.0.6: Factor w/ 23 levels "0","1","2","3",..: 8 1 2 5 1 10 19 1 2 1 ...
  #..$ seurat_clusters: Factor w/ 23 levels "0","1","2","3",..: 8 1 2 5 1 10 19 1 2 1 ...
  #..$ celltype       : Factor w/ 12 levels "HSC","MPP","GMP",..: 12 6 2 3 6 6 8 6 2 6 ...
  #  active.assay: chr "RNA"
  #  active.ident: Factor w/ 12 levels "cDC","T","B",..: 1 7 11 10 7 7 5 7 11 7 ...
  #..- attr(*, "names")= chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  #  graphs      :List of 2
  #..$ RNA_nn :Formal class 'Graph' [package "SeuratObject"] with 7 slots
  ####  assay.used: chr "RNA"
  ####  i         : int [1:297760] 0 4641 7134 7182 8095 8333 8648 9912 10513 10710 ...
  ####  p         : int [1:14889] 0 16 62 79 106 122 148 155 162 192 ...
  ####  Dim       : int [1:2] 14888 14888
  ####  Dimnames  :List of 2
  ####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ####  x         : num [1:297760] 1 1 1 1 1 1 1 1 1 1 ...
  ####  factors   : list()
  #..$ RNA_snn:Formal class 'Graph' [package "SeuratObject"] with 7 slots
  ####  assay.used: chr "RNA"
  ####  i         : int [1:1136260] 0 1248 1315 1421 1846 1895 2067 2779 3496 3563 ...
  ####  p         : int [1:14889] 0 59 161 250 331 430 493 530 625 688 ...
  ####  Dim       : int [1:2] 14888 14888
  ####  Dimnames  :List of 2
  ####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  ####  x         : num [1:1136260] 1 0.1111 0.1111 0.1429 0.0811 ...
  ####  factors   : list()
  #  neighbors   : list()
  #  reductions  :List of 3
  #..$ pca    :Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  ####  cell.embeddings           : num [1:14888, 1:50] 9.89 -11.04 10.27 -1.52 -11.46 ...
  ####..- attr(*, "dimnames")=List of 2
  #####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  #####..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
  ####  feature.loadings          : num [1:2000, 1:50] 0.01341 0.01331 0.01648 0.01724 0.00112 ...
  ####..- attr(*, "dimnames")=List of 2
  #####..$ : chr [1:2000] "Hbb-bs" "Hba-a2" "Hba-a1" "Hbb-bt" ...
  #####..$ : chr [1:50] "PC_1" "PC_2" "PC_3" "PC_4" ...
  ####  feature.loadings.projected: num[0 , 0 ] 
  ####  assay.used                : chr "RNA"
  ####  global                    : logi FALSE
  ####  stdev                     : num [1:50] 11.91 8.57 7.72 6.5 5.73 ...
  ####  jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  ######  empirical.p.values     : num[0 , 0 ] 
  ######  fake.reduction.scores  : num[0 , 0 ] 
  ######  empirical.p.values.full: num[0 , 0 ] 
  ######  overall.p.values       : num[0 , 0 ] 
  ####  misc                      :List of 1
  ####..$ total.variance: num 1328
  ####  key                       : chr "PC_"
  #..$ harmony:Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  ####  cell.embeddings           : num [1:14888, 1:50] 11.882 -11.02 11.657 -0.837 -11.44 ...
  ####..- attr(*, "dimnames")=List of 2
  #####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  #####..$ : chr [1:50] "harmony_1" "harmony_2" "harmony_3" "harmony_4" ...
  ####  feature.loadings          : num[0 , 0 ] 
  ####  feature.loadings.projected: num[0 , 0 ] 
  ####  assay.used                : chr "RNA"
  ####  global                    : logi FALSE
  ####  stdev                     : num [1:50] 11.89 8.42 7.49 5.88 5.55 ...
  ####  jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  ######  empirical.p.values     : num[0 , 0 ] 
  ######  fake.reduction.scores  : num[0 , 0 ] 
  ######  empirical.p.values.full: num[0 , 0 ] 
  ######  overall.p.values       : num[0 , 0 ] 
  ####  misc                      : list()
  ####  key                       : chr "harmony_"
  #..$ umap   :Formal class 'DimReduc' [package "SeuratObject"] with 9 slots
  ####  cell.embeddings           : num [1:14888, 1:2] 6.02 -2.9 10.67 1.61 -2.61 ...
  ####..- attr(*, "scaled:center")= num [1:2] -0.426 -0.17
  ####..- attr(*, "dimnames")=List of 2
  #####..$ : chr [1:14888] "Kras_AAACCCACAGAAACCG-1" "Kras_AAACCCACAGTATACC-1" "Kras_AAACCCAGTACGTGAG-1" "Kras_AAACCCATCGTTAGTG-1" ...
  #####..$ : chr [1:2] "UMAP_1" "UMAP_2"
  ####  feature.loadings          : num[0 , 0 ] 
  ####  feature.loadings.projected: num[0 , 0 ] 
  ####  assay.used                : chr "RNA"
  ####  global                    : logi TRUE
  ####  stdev                     : num(0) 
  ####  jackstraw                 :Formal class 'JackStrawData' [package "SeuratObject"] with 4 slots
  ######  empirical.p.values     : num[0 , 0 ] 
  ######  fake.reduction.scores  : num[0 , 0 ] 
  ######  empirical.p.values.full: num[0 , 0 ] 
  ######  overall.p.values       : num[0 , 0 ] 
  ####  misc                      : list()
  ####  key                       : chr "UMAP_"
  #  images      : list()
  #  project.name: chr "SeuratProject"
  #  misc        : list()
  #  version     :Classes 'package_version', 'numeric_version'  hidden list of 1
  #..$ : int [1:3] 4 1 0
  #  commands    :List of 7
  #..$ NormalizeData.RNA        :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "NormalizeData.RNA"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:29:09"
  ####  assay.used : chr "RNA"
  ####  call.string: chr "NormalizeData(sc_filt)"
  ####  params     :List of 5
  ####..$ assay               : chr "RNA"
  ####..$ normalization.method: chr "LogNormalize"
  ####..$ scale.factor        : num 10000
  ####..$ margin              : num 1
  ####..$ verbose             : logi TRUE
  #..$ FindVariableFeatures.RNA :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "FindVariableFeatures.RNA"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:29:17"
  ####  assay.used : chr "RNA"
  ####  call.string: chr [1:2] "FindVariableFeatures(sc_n, selection.method = \"vst\", nfeatures = 2000, " "    verbose = FALSE)"
  ####  params     :List of 12
  ####..$ assay              : chr "RNA"
  ####..$ selection.method   : chr "vst"
  ####..$ loess.span         : num 0.3
  ####..$ clip.max           : chr "auto"
  ####..$ mean.function      :function (mat, display_progress)  
  ####..$ dispersion.function:function (mat, display_progress)  
  ####..$ num.bin            : num 20
  ####..$ binning.method     : chr "equal_width"
  ####..$ nfeatures          : num 2000
  ####..$ mean.cutoff        : num [1:2] 0.1 8
  ####..$ dispersion.cutoff  : num [1:2] 1 Inf
  ####..$ verbose            : logi FALSE
  #..$ ScaleData.RNA            :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "ScaleData.RNA"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:29:19"
  ####  assay.used : chr "RNA"
  ####  call.string: chr "ScaleData(sc_n)"
  ####  params     :List of 10
  ####..$ features          : chr [1:2000] "Hbb-bs" "Hba-a2" "Hba-a1" "Hbb-bt" ...
  ####..$ assay             : chr "RNA"
  ####..$ model.use         : chr "linear"
  ####..$ use.umi           : logi FALSE
  ####..$ do.scale          : logi TRUE
  ####..$ do.center         : logi TRUE
  ####..$ scale.max         : num 10
  ####..$ block.size        : num 1000
  ####..$ min.cells.to.block: num 3000
  ####..$ verbose           : logi TRUE
  #..$ RunPCA.RNA               :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "RunPCA.RNA"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:29:31"
  ####  assay.used : chr "RNA"
  ####  call.string: chr "RunPCA(sc_n, npcs = 50, verbose = FALSE)"
  ####  params     :List of 10
  ####..$ assay          : chr "RNA"
  ####..$ npcs           : num 50
  ####..$ rev.pca        : logi FALSE
  ####..$ weight.by.var  : logi TRUE
  ####..$ verbose        : logi FALSE
  ####..$ ndims.print    : int [1:5] 1 2 3 4 5
  ####..$ nfeatures.print: num 30
  ####..$ reduction.name : chr "pca"
  ####..$ reduction.key  : chr "PC_"
  ####..$ seed.use       : num 42
  #..$ FindNeighbors.RNA.harmony:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "FindNeighbors.RNA.harmony"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:30:47"
  ####  assay.used : chr "RNA"
  ####  call.string: chr "FindNeighbors(sce.har, dims = pc.num, reduction = \"harmony\")"
  ####  params     :List of 17
  ####..$ reduction      : chr "harmony"
  ####..$ dims           : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
  ####..$ assay          : chr "RNA"
  ####..$ k.param        : num 20
  ####..$ return.neighbor: logi FALSE
  ####..$ compute.SNN    : logi TRUE
  ####..$ prune.SNN      : num 0.0667
  ####..$ nn.method      : chr "annoy"
  ####..$ n.trees        : num 50
  ####..$ annoy.metric   : chr "euclidean"
  ####..$ nn.eps         : num 0
  ####..$ verbose        : logi TRUE
  ####..$ force.recalc   : logi FALSE
  ####..$ do.plot        : logi FALSE
  ####..$ graph.name     : chr [1:2] "RNA_nn" "RNA_snn"
  ####..$ l2.norm        : logi FALSE
  ####..$ cache.index    : logi FALSE
  #..$ FindClusters             :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "FindClusters"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:30:50"
  ####  assay.used : chr "RNA"
  ####  call.string: chr [1:2] "FindClusters(sce.har, graph.name = \"RNA_snn\", resolution = 0.6, " "    algorithm = 1)"
  ####  params     :List of 10
  ####..$ graph.name      : chr "RNA_snn"
  ####..$ modularity.fxn  : num 1
  ####..$ resolution      : num 0.6
  ####..$ method          : chr "matrix"
  ####..$ algorithm       : num 1
  ####..$ n.start         : num 10
  ####..$ n.iter          : num 10
  ####..$ random.seed     : num 0
  ####..$ group.singletons: logi TRUE
  ####..$ verbose         : logi TRUE
  #..$ RunUMAP.RNA.harmony      :Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
  ####  name       : chr "RunUMAP.RNA.harmony"
  ####  time.stamp : POSIXct[1:1], format: "2022-08-17 18:31:35"
  ####  assay.used : chr "RNA"
  ####  call.string: chr "RunUMAP(sce.har, dims = pc.num, reduction = \"harmony\")"
  ####  params     :List of 26
  ####..$ dims                : int [1:30] 1 2 3 4 5 6 7 8 9 10 ...
  ####..$ reduction           : chr "harmony"
  ####..$ assay               : chr "RNA"
  ####..$ slot                : chr "data"
  ####..$ umap.method         : chr "uwot"
  ####..$ return.model        : logi FALSE
  ####..$ n.neighbors         : int 30
  ####..$ n.components        : int 2
  ####..$ metric              : chr "cosine"
  ####..$ learning.rate       : num 1
  ####..$ min.dist            : num 0.3
  ####..$ spread              : num 1
  ####..$ set.op.mix.ratio    : num 1
  ####..$ local.connectivity  : int 1
  ####..$ repulsion.strength  : num 1
  ####..$ negative.sample.rate: int 5
  ####..$ uwot.sgd            : logi FALSE
  ####..$ seed.use            : int 42
  ####..$ angular.rp.forest   : logi FALSE
  ####..$ densmap             : logi FALSE
  ####..$ dens.lambda         : num 2
  ####..$ dens.frac           : num 0.3
  ####..$ dens.var.shift      : num 0.1
  ####..$ verbose             : logi TRUE
  ####..$ reduction.name      : chr "umap"
  ####..$ reduction.key       : chr "UMAP_"
  #  tools       : list()
