install.packages("remotes")
remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)

DefaultAssay(sce.harm) = "RNA"
VlnPlot(sce.harm, features = c("Hpd"),pt.size=0)
DefaultAssay(HSC) = "RNA"
VlnPlot(HSC, features = c("Iigp1"),pt.size=0)
VlnPlot(HSC, features = c("Hlf","Ifitm1","Flt3","Dntt","Gata1","Car1","Mki67","Elane"
                          ,"Mpo","Pf4","Vwf","Itga2b",'Cebpe'),stacked=T,pt.size=0)
VlnPlot(HSC, features = c("Hlf","Ifitm1","Flt3","Dntt","Gata1","Car1","Mki67","Elane"
                          ,"Mpo","Pf4","Vwf","Itga2b",'Cebpe'),stacked=T,pt.size=0)
VlnPlot(HSC, features = c("Cd74","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                          "Cd34","Cd63",'Cd9',"Cd117"),
        stacked=T,pt.size=0)
VlnPlot(HSC, features = c("Egr1","Gata2","Junb"),stacked=T,pt.size=0)
ggsave("vlnHSC_共享.pdf",width = 20,height = 60,units = "cm")
VlnPlot(sce.harm, 
        features = c("Cd74","Cd27","Cd52","Cd93","Cd53","Cd79a","Cd55","Cd59a","Cd69",
                     "Cd34","Cd63",'Cd9',"Cd117"),
        pt.size = 0)

VlnPlot(HSC, 
        features = c("Cd74","Cd27","Cd52","Cd93","Cd53","Cd55","Cd59a","Cd69",
                     "Cd34","Cd63",'Cd9',"Cd117"),
        pt.size = 0,
        ncol = 2,
        split.by = "orig.ident")

VlnPlot(HSC, 
        features = c("Txnip","Mllt3","Socs2","Mpl", "Mycn","Cdkn1c","Ndn"),
        pt.size = 0)
ggsave("vln11.pdf",width = 20,height = 40,units = "cm")

VlnPlot(HSC, 
        features = c("Iigp1"),
        pt.size = 0,
        split.by = "orig.ident")
ggsave("vln_Iigp1.pdf",width = 40,height = 40,units = "cm")

VlnPlot(HSC, 
        features = c("Phlda1"),
        pt.size = 0,
        split.by = "orig.ident")
ggsave("vln_Phlda1.pdf",width = 40,height = 40,units = "cm")

VlnPlot(sce.harm, 
        features = c("Hpd"),
        pt.size = 0)
ggsave("vln10.pdf",width = 40,height = 40,units = "cm")



load("D:/aJMML/scRNA1/sc_seurat_integr.RData")
DefaultAssay(sc_integr) = "RNA"
VlnPlot(sc_integr, 
        features = c("HPD"),
        pt.size = 0)
ggsave("vln10.pdf",width = 40,height = 40,units = "cm")


VlnPlot(HSC, 
        features = c("KIf1","Gm15915","Tspo2","Mfsd2b", "Gata1","Il1rl1","Aqp1","Car1","Trib2"
                     ,"Apoe","Abcb4","Csrp3","Samd14","Atp1b2","Cd69","Cd55","Cd59a"),split.by = "orig.ident",pt.size = 0)
ggsave("vln11.pdf",width = 20,height = 40,units = "cm")