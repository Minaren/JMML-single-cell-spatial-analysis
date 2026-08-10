#### Preparation ####
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(220625)

# Suppress package startup messages and load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(tidyverse)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(cowplot)
  library(stringr)
})

#### I. Merging Multiple Single-Cell Samples ####
## 1. Read and Merge Data
### 1.1 Read Data
dir <- c('D:/aJMML/JMML/Kras', 'D:/aJMML/JMML/WT')
samples_name <- c('Kras', 'WT')

### 1.2 Batch Creation of Seurat Objects
scRNAlist <- lapply(1:length(dir), function(i) {
  counts <- Read10X(data.dir = dir[i])
  seurat_obj <- CreateSeuratObject(counts, project = samples_name[i], min.cells = 3, min.features = 200)
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = samples_name[i])
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj
})

### 1.3 Merge Seurat Objects
sc_merge <- do.call(merge, scRNAlist)
table(sc_merge$orig.ident)
head(sc_merge@meta.data, 5)
save(sc_merge, file = 'sc_merge.Rdata')

#### II. Data Quality Control ####
## 2.1 Calculate Gene Ratios
scRNA <- sc_merge
theme.set2 <- theme(axis.title.x = element_blank())
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
group <- "orig.ident"

plots <- lapply(plot.features, function(feature) {
  VlnPlot(scRNA, group.by = group, pt.size = 0, features = feature) + theme.set2 + NoLegend()
})

violin <- wrap_plots(plots = plots, nrow = 2)
ggsave("vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

## 2.2 Set Quality Control Standards
minGene <- 500
maxGene <- 6000
pctMT <- 10

## 2.3 Visualization After Filtering
sc_filt <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

plots <- lapply(plot.features, function(feature) {
  VlnPlot(sc_filt, group.by = group, pt.size = 0, features = feature) + theme.set2 + NoLegend()
})

violin <- wrap_plots(plots = plots, nrow = 2)
ggsave("vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8)
save(sc_filt, file = "sc_filt.RData")

#### III. Data Integration ####
## 3.1 Prepare Data
load("D:/aJMML/scRNA/sc_filt.RData")
sc_n <- NormalizeData(sc_filt)
sc_n <- FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_n <- ScaleData(sc_n)
sc_n

## 3.2 Integrate Data Using Harmony
sc_n <- RunPCA(sc_n, npcs = 50, verbose = FALSE)
ElbowPlot(sc_n, ndims = 50)
ggsave('elbowplot_n.pdf', width = 6.5, height = 5, dpi = 500)
ggsave('elbowplot_n.png', width = 6.5, height = 5, dpi = 500)

sce_ha <- sc_n
pc.num <- 1:30
sce_har <- RunHarmony(sce_ha, group.by.vars = "orig.ident", project.dim = FALSE, plot_convergence = TRUE)

## 3.3 Clustering
sce.har <- sce_har
sce.har <- FindNeighbors(sce.har, dims = pc.num, reduction = "harmony")
sce.har <- FindClusters(sce.har, graph.name = "RNA_snn", resolution = 0.6, algorithm = 1)
table(sce.har@active.ident)

## 3.4 Dimensionality Reduction and Visualization
sce.harm <- RunUMAP(sce.har, dims = pc.num, reduction = "harmony")
save(sce.harm, file = "sce_harm_tu.RData")

DimPlot(sce.harm, label = TRUE, repel = TRUE) + NoLegend()
ggsave('umap.pdf', width = 6, height = 6, dpi = 500)
ggsave('umap.png', width = 6, height = 6, dpi = 500)

DimPlot(sce.harm, label = FALSE, group.by = 'orig.ident')
ggsave('umap_indi.pdf', width = 7, height = 6, dpi = 500)
ggsave('umap_indi.png', width = 7, height = 6, dpi = 500)

#### IV. Cell Annotation ####
## 4.1 Gene Check
rm(list = ls())
load("D:/aJMML/scRNA1/sce_harm_tu.RData")

genes_to_check <- c("Prss34","Cd200r3","Lmo4","Vpreb3","Cd79a","Cd79b","Cd74","Cd19",
                    "Prg2","Itga2b","Plek","Cebpe","Lmo4","Prss34",
                    "Blvrb","Rhd","Hmbs","Gata1","Car1",
                    "Hba-a1","Hbb-bt","Hbb-bs","Blvrb","Hba-a2","Hmbs",
                    "Cd68","Cd74","Mpeg1","Adgre1","Csf1r","C1qc",
                    "Cd74","Mpeg1","Cd68","Csf1r","Adgre1",
                    "Csf1r","Cd74","Mpeg1","Cd68","S100a4",
                    "Pf4","Ctla2a","Itga2b","Gata2","Plek","Cd9",
                    "Prtn3","Elane","Mpo","Ctsg",
                    "Elane","Mpo","Prtn3",
                    "Ms4a6c","F13a1","Irf8","Csf1r",
                    "Cd34","Flt3","Ctla2a",
                    "Camp","Ngp","Lcn2","Cebpe","Fcnb",
                    "Mpo","Cd177",
                    "Lcn2",
                    "Ngp","Lcn2","Cd177","Cebpe","Ltf",
                    "Ngp","Lcn2","Cd177","Cebpe",
                    "Mmp8","Lcn2","Ngp","Camp","Retnlg",
                    "Lmo4","Runx1",
                    "Slamf1","Kit","Cd34","Sca1","Cd48","Procr","Ifitm1","Mecom", "Hoxa9", "Mycn", "Hlf","Cd201","Cd47",
                    "Ms4a1",
                    "Il7r,Ccr7,Cd8a")

DefaultAssay(sce.harm) <- "RNA"
DotPlot(sce.harm, features = unique(genes_to_check)) + coord_flip()
ggsave("har_gene_show.pdf", width = 14, height = 12, dpi = 500)
ggsave("har_gene_show.png", width = 14, height = 12, dpi = 500)

allmarkers <- FindAllMarkers(sce.harm, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)
write.csv(allmarkers, 'allmarkers_celltype_mouse.csv')

top10 <- allmarkers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
write.csv(top10, 'top10_celltype_mouse.csv')
top10 <- read.csv("top10.csv", header = TRUE)

p <- DoHeatmap(sce.harm, features = top10$gene) + NoLegend()
ggsave("markers.heatmap_celltype_mouse.png", plot = p, width = 17, height = 17)
ggsave("markers.heatmap_celltype_mouse_top10.pdf", plot = p, width = 30, height = 30)

p <- DotPlot(sce.harm, features = unique(top10$gene), assay = 'RNA') + coord_flip()
plotd <- p + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))

ggsave("markers.dotplot.png", plot = plotd, width = 17, height = 49)
ggsave("markers.dotplot.pdf", plot = plotd, width = 17, height = 49)

## 4.2 Enrichment Analysis
rm(list = ls())
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

allmarkers <- read.csv("D:/aJMML/scRNA/allmarkers.csv")
ids <- bitr(allmarkers$gene, 'SYMBOL', 'ENTREZID', 'org.Mm.eg.db')
allmarkers <- merge(allmarkers, ids, by.x = 'gene', by.y = 'SYMBOL')

entrezid <- allmarkers$ENTREZID
GOdata <- enrichGO(entrezid, OrgDb = 'org.Mm.eg.db', ont = 'BP', pvalueCutoff = 0.05, readable = TRUE)
KEGGdata <- enrichKEGG(entrezid, organism = 'mmu', pvalueCutoff = 0.05)

write.csv(GOdata, 'go_data.csv')
write.csv(KEGGdata, 'kegg_data.csv')
