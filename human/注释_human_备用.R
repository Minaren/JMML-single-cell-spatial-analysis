library(Seurat) 
library(tidyverse) 
load("D:/aJMML/scRNA1/sc_seurat_integr.RData")
DefaultAssay(sc_integr) = "RNA"

#HSC
genes_to_check =c("ITGAM","IL3RA","CD19","MS4A1","CD3D","CD33","CD34",
                  "CD38", "CD4","PTPRC","CD8A","THY1","ADMP","ERC2","CYP2C8","KIT",
                  "NECTIN3","ITGA6","CD135","CD110")
genes_to_check =c("CD2","MME","ITGAM","CD14","FUT4","CD19","NCAM1","GYPA","PTPRC","CD34",
                 "CD38","THY1","ITGA6","IL3RA","CD7","FLT3","MPL","CEACAM1","IL5RA","CCR3")

genes_to_check =c("CD2","ITGAM","ITGAL","ANPEP","CD14","FUT4","CD19","CD33","IL6R","IL7R","MME","CD7","CD127",
                  "TFRC","EPOR","ITGA2B","GP1BA","PTPRC","CD34","CD38","THY1","ITGA6","KIT","FLT3","IL3RA")
DefaultAssay(sc_integr) = "RNA"
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_HSC5.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show_HSC5.png",width = 14,height = 12,dpi = 500)

#LMPP
genes_to_check =c("KIAA0087","PRSS1","RXFP1","PRSS3","LTB","COBLL1","KCNK17","HOXA9",
                  "LZTFL1","DCK","ARID5B")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_LMPP.pdf",width = 14,height = 12,dpi = 500)

#CMP
genes_to_check =c("HTRA3","C1QTNF4","PRAM1","IGLL1","CSF3R","MZB1","KBTBD11","KBTBD11",
                  "SPARC","PHGDH","RAPGEF1","EID2B","LAT2","IKBKE","SCMH1","SPTAN1",
                  "PKM","NCDN","P2RY8","PTPRCAP","DOK3","HCST","TRIM14","KCNAB2",
                  "RAP1GAP2","CDCA7","TNFSF13B","ADA","MAPKAPK3","LSP1","EDEM2","IMPA2
                  ","TMEM106C")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_CMP.pdf",width = 14,height = 12,dpi = 500)

#GMP
genes_to_check =c("CTSG","AZU1","ELANE","PRTN3","CST7","SRGN","RNASE3","RNASE2","SLPI","COL23A1")
DotPlot(sc_integr,features = unique(genes_to_check)) + coord_flip()
ggsave("har_gene_show_GMP.pdf",width = 14,height = 12,dpi = 500)

#MEP
genes_to_check =c("TGM2","PPP1R14A","FAM129B","TNFRSF25","SPTA1","CXADR","ECE1","PVT1",
"STEAP3","MTSS1","ATP7B","BLVRB","AKR1C1","PRKAR2B","TUBB2A","NECAB1","AMMECR1","SGPP1",
"GBGT1","ELOVL6")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_MEP.pdf",width = 14,height = 12,dpi = 500)

#CLP
genes_to_check =c("LTA","TBCD","ZCCHC7","ARHGEF7","ZADH2","PRPF39","TOPBP1","UFM1",
"PTEN","VPREB3","CALM1","BAZ1A","UBE2G1","PMAIP1","FAM76B","NPY","SMAD5")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show7.pdf",width = 14,height = 12,dpi = 500)

#CLP
genes_to_check =c("LTA","TBCD","ZCCHC7","ARHGEF7","ZADH2","PRPF39","TOPBP1","UFM1",
                  "PTEN","VPREB3","CALM1","BAZ1A","UBE2G1","PMAIP1","FAM76B","NPY","SMAD5")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_CLP.pdf",width = 14,height = 12,dpi = 500)


#B13
genes_to_check =c("IGKC","CD79A","CD37","CD79B","CD52","CXCR4","IGHD","SPIB","HERPUD1",
                  "VPREB3","FAM129C","SMIM14","CD83","CD24","MZB1","NCF1","RALGPS2","CLEC2D","LTB",
                  "BCL2","BCL6","LMO2","MME", "MS4A1", "CD34","CD19","CD24")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_B.pdf",width = 14,height = 12,dpi = 500)



#Macrophage
genes_to_check =c("SPN","CFD","PLK2","LGALS9","LGALS3","EREG","MYO10","PLAUR","OLR1","FABP2","CAPG","MPO",
                  "SRD5A1","BTG1","CD68","ADGRE1","IRF5","CD14","PECAM1","CD33","SIRPA","ITGAM","FUT4","CD40")


DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show12.pdf",width = 14,height = 12,dpi = 500)

#	Monocyte 2
genes_to_check =c("IL1B","LYZ","MNDA","CSTA","PLAUR","S100A6","TSPO","SOD2","BCL2A1"
                  ,"TYROBP","CTSS","LGALS1","NAMPT","CCL4","NFKBIA","PLBD1")
genes_to_check =c("ITGAM","FUT4","CD40","CD63","CD68","LILRA5","CSF1R","SIGLEC1"
                  ,"CCR5","TREM1","CXCL10","CXCL11","CCL4","CXCL9","EMR1","GPNMB","CXCL2"
                  ,"VSIG4","VSIR","CD14","CD33","SIRPA")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show13.pdf",width = 14,height = 12,dpi = 500)

#	Eosinophil 4,12
genes_to_check =c("ITGA5","CSF2RB","LGALS9","LTC4S","CCL20","STX4","EPX","TM9SF3",
                  "GAB3","CYBA","CR1","CD63","CD69","CD40","F2R","CD244","CAMK1D")
genes_to_check =c("ITGAM","ANPEP","FUT4","FCGR2A","CD33","MBP","RNASE2","EPX","CD9",
                  "CD24","CR1","SPN","FCGR1A","CSF2RA","IL3RA","IL5RA","IL6R","CCR3",
                  "GGT1")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show14.pdf",width = 14,height = 12,dpi = 500)


#Basophil 6,12
genes_to_check =c("CSF2RA","FCER1G","LTC4S","CD63","CD69","CD244","CDK9","IL3RA","NCF4","CD59","RNASE3","SELPLG")
genes_to_check =c("ITGAM","ANPEP","FUT4","FCGR2A","CD33","CD9","ITGAL","ANPEP","FCGR3A"
                  ,"IL2RA",'CD33','CD38','SPN','CD88','IL3RA','IL5RA','CD40LG','CCR2'
                  ,"ENPP3","IL18R1","TLR2","TLR4","TLR6","TLR1","TLR9")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show15.pdf",width = 14,height = 12,dpi = 500)


#	Neutrophil 4,12
genes_to_check =c("CHST2","PNP","CFP","KHDRBS1","CLIC1","RAC2","VCL","CSF2RB","CHST1",
                  "NAPRT1","CAP1","CAPG","MPO","RAP1B","ATP6V1B2","SYT2","GNAI2","ICMT")
genes_to_check =c("ITGAM","ANPEP","FUT4","FCGR2A","MME","CD24","CR1",
                  "SPN","CEACAM1","CEACAM8","CEACAM6","CEACAM3","FCAR","CD93","NECTIN2","CSF3R"
                  ,"CSF2RA","BTS1","CD177","CXCR1","TLR2","TLR4","TLR6","S100A9","S100A8",
                  "IL6","IL12","TNF","IL1")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show16.pdf",width = 14,height = 12,dpi = 500)


#	Erythroid 4,12
genes_to_check =c("EPOR","KLF1","TFR2","CSF2RB","APOE","APOC1","CNRIP1","CD35
                  ","CD36",'TFRC',"ENG","ALDH","SLC4A1",
                  "TFRC","EPOR","GPA","HBG1","HBG2","CD71")
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show17.pdf",width = 14,height = 12,dpi = 500)





genes_to_check =c("AVP","MLLT3","CAVIN1","HOPX","SPINK2","RAMP1","SELL",
                  "CD52", "CTHRC1","FAM43A","GIMAP7","MECOM","ABCB1","PTAFR","CALCRL","ADGRG1","CRHBP")
DefaultAssay(sc_integr) = "RNA"
DotPlot(sc_integr,features = unique(genes_to_check), cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_HSPC.pdf",width = 14,height = 12,dpi = 500)
ggsave("har_gene_show1_HSPC.png",width = 14,height = 12,dpi = 500)