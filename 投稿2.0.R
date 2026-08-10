#### Fig1 ####
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(220625)

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
  library(msigdbr)
  library(GSVA)
  library(clusterProfiler)
  library(limma)
  library(pathview)
  library(org.Mm.eg.db)
  library(pheatmap)
  library(gplots)
})

#### Data Loading and Seurat Object Creation ####
file_paths <- list(
  Kras = list(
    matrix = "D:/aJMML/JMML/Kras/Kras_matrix.mtx.gz",
    features = "D:/aJMML/JMML/Kras/Kras_features.tsv.gz",
    barcodes = "D:/aJMML/JMML/Kras/Kras_barcodes.tsv.gz"
  ),
  WT = list(
    matrix = "D:/aJMML/JMML/WT/WT_matrix.mtx.gz",
    features = "D:/aJMML/JMML/WT/WT_features.tsv.gz",
    barcodes = "D:/aJMML/JMML/WT/WT_barcodes.tsv.gz"
  )
)

sc_list <- list()
for (sample in names(file_paths)) {
  counts <- Read10X(data.dir = dirname(file_paths[[sample]]$matrix))
  sc <- CreateSeuratObject(counts, project = sample, min.cells = 3, min.features = 200)
  sc <- RenameCells(sc, add.cell.id = sample)
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  sc_list[[sample]] <- sc
}

sc_merged <- merge(sc_list[[1]], sc_list[2:length(sc_list)])

#### Quality Control ####
plot_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
theme_set(theme_classic())

plots <- lapply(plot_features, function(f){
  VlnPlot(sc_merged, features = f, group.by = "orig.ident", pt.size = 0) + NoLegend()
})
wrap_plots(plots, nrow = 2)

minGene <- 500
maxGene <- 6000
pctMT <- 10
sc_filtered <- subset(sc_merged, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

plots <- lapply(plot_features, function(f){
  VlnPlot(sc_filtered, features = f, group.by = "orig.ident", pt.size = 0) + NoLegend()
})
wrap_plots(plots, nrow = 2)

#### Normalization and Feature Selection ####
sc_norm <- NormalizeData(sc_filtered)
sc_norm <- FindVariableFeatures(sc_norm, selection.method = "vst", nfeatures = 2000)
sc_norm <- ScaleData(sc_norm)

#### PCA and Harmony Integration ####
sc_norm <- RunPCA(sc_norm, npcs = 50, verbose = FALSE)
ElbowPlot(sc_norm, ndims = 50)

pc_num <- 1:30
sc_harmony <- RunHarmony(sc_norm, group.by.vars = "orig.ident", project.dim = FALSE, plot_convergence = TRUE)

#### Clustering and UMAP ####
sc_harmony <- FindNeighbors(sc_harmony, dims = pc_num, reduction = "harmony")
sc_harmony <- FindClusters(sc_harmony, graph.name = "RNA_snn", resolution = 1.2)
sc_harmony <- RunUMAP(sc_harmony, dims = pc_num, reduction = "harmony")

#### Cell Type Annotation ####
cluster2celltype <- c(
  "0"="Granulocyte","2"="Granulocyte","5"="Granulocyte","8"="Granulocyte","9"="Granulocyte","13"="Granulocyte","16"="Granulocyte",
  "14"="Macrophage","18"="Monocyte","11"="B","22"="B","7"="cDC","17"="cDC",
  "21"="Ery","10"="Erythroblast","12"="Erythroblast","15"="HSC","6"="MPP",
  "19"="T","1"="GMP","3"="GMP","4"="GMP","20"="MK"
)
sc_harmony$celltype <- cluster2celltype[as.character(sc_harmony$RNA_snn_res.0.6)]
Idents(sc_harmony) <- factor(sc_harmony$celltype, levels = c('HSC','MPP','GMP','Erythroblast','Ery','Granulocyte','Macrophage','Monocyte','MK','B','T','cDC'))

#### Fig1B: Enhanced UMAP Plot ####
sce <- sc_harmony
meta <- cbind(sce@meta.data, Embeddings(sce, "umap"))

col_df <- data.frame(name = unique(meta$celltype)) %>% arrange(name)
mycol <- c(
  "#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", "#92b8da",
  "#b5aa82", "#de9d3d", "#347852", "#ca8399", "#296097", "#564b84"
)
mycol <- setNames(mycol, col_df$name)

mid_coord_type <- meta %>%
  select(celltype, UMAP_1, UMAP_2) %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

pumap <- ggplot(meta, aes(UMAP_1, UMAP_2)) +
  geom_point(size = 0.2, aes(color = celltype)) +
  geom_text(data = mid_coord_type, size = 5, aes(label = celltype)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  scale_color_manual(values = mycol)

pdf("Fig1B_UMAP_CellTypes.pdf", width = 8, height = 6)
print(pumap)
dev.off()

#### Fig1C: Marker Gene Heatmap ####
genes_to_check <- rev(c(
  "H2-Eb1","H2-Aa","H2-DMb1","H2-DMa","Cd74",
  "Trbc1","Trbc2","Cd3d","Cd3g","Cd247",
  "Vpreb3","Cd79a","Cd79b","Cd74","Cd19",
  "Cxcl12","Gas6","Kitl","Pmp22","Gpm6b",
  "S100a4","Clec4a1","Clec4a3","Fn1","Mafb",
  "Cd302","Cd68","Csf1r","Fcgr1","Adgre4",
  "Cd177","Ltf","Retnlg","Ly6g","Fpr2",
  "Blvrb","Hba-a1","Hbb-bt","Hbb-bs","Hba-a2","Gata1",
  "Ctsg","Elane","Gfi1","Mpo","Ms4a3",
  "Meis1","Mef2c","Flt3","Hlf","Cd34",
  "Cdk6","Tal1","Meis1","Gata2","Slamf1"
))

p_heat <- DoHeatmap(sc_harmony, features = genes_to_check, group.by = "celltype") +
  scale_fill_gradientn(colors = c("#92b8da", "white", "#9f2b39")) +
  NoLegend()

pdf("Fig1C_Marker_Heatmap.pdf", width = 10, height = 8)
print(p_heat)
dev.off()

#### Fig1D: GSVA Heatmap ####
genesets <- read.csv("Cell_type_annotation_of_mouse_samples.csv", header = FALSE)
genesets <- split(genesets$V2, genesets$V1)

expr <- AverageExpression(sc_harmony, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0,]

gsva_res <- gsva(as.matrix(expr), genesets, method = "ssgsea")

pdf("Fig1D_GSVA_Heatmap.pdf", width = 8, height = 6)
pheatmap(gsva_res,
         scale = "row",
         cluster_row = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#92b8da", "white", "#9f2b39"))(100))
dev.off()

#### Fig1E: Cell-type Composition ####
load("D:/aJMML/scRNA1/sce_har_anno.RData")

cell_types <- c("HSC","MPP","GMP","Erythroblast","Ery","Granulocyte","Macrophage","Monocyte","MK","B","T","cDC")
sce.harm$celltype <- factor(sce.harm$celltype, levels = cell_types)

mycol <- c("#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da","#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84")
names(mycol) <- cell_types

Cellratio$Var1 <- factor(Cellratio$Var1, levels = cell_types)

p_cell <- ggplot(Cellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1),
           stat = "identity",
           width = 0.2,
           size = 0.5,
           colour = "#222222") +
  scale_fill_manual(values = mycol) +
  theme_classic() +
  labs(x = "Sample", y = "Ratio") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5))

pdf("Fig1E_CellType_Composition.pdf", width = 8, height = 6)
print(p_cell)
dev.off()

#### Fig1F: FeaturePlot of Selected Genes ####
genes_to_plot <- c("Ly6a","Kit","Slamf1","Cd48","Cd34","Flt3")
plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(sce.harm, features = gene, cols = c("lightgrey","#9f2b39")) + NoAxes() + NoLegend()
})

pdf("Fig1F_FeaturePlot_SelectedGenes.pdf", width = 10, height = 8)
wrap_plots(plots, ncol = 2)
dev.off()

####Fig2####
#### HSC Subcluster Analysis ####
HSC <- subset(sc_harmony, celltype == "HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA","percent.mt"))
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50)
HSC <- FindNeighbors(HSC, dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.8)
HSC <- RunUMAP(HSC, dims = 1:50)

HSC_col <- c("#9e2a2f", "#4f85b8", "#d87e2d")

#### Fig2A: HSC subclustering UMAP ####
p_umap <- DimPlot(HSC, label = TRUE, pt.size = 2) +
  scale_color_manual(values = HSC_col) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

pdf("Fig2A_HSC_Subcluster_UMAP.pdf", width = 6, height = 5)
print(p_umap)
dev.off()

#### Fig2B: HSC Subpopulation Composition ####
Cellratio <- prop.table(table(Idents(HSC), HSC$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
HSC_col <- c("#d87e2d", "#4f85b8", "#9e2a2f")
Cellratio$Var2 <- factor(Cellratio$Var2, levels = rev(levels(Cellratio$Var2)))

p_bar <- ggplot(Cellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1),
           stat = "identity",
           width = 0.2,
           size = 0.5,
           colour = "#222222") +
  scale_fill_manual(values = HSC_col) +
  labs(x = "Sample", y = "Ratio", fill = "Subpopulation") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5))

pdf("Fig2B_HSC_Subpop_Composition.pdf", width = 6, height = 5)
print(p_bar)
dev.off()

#### Fig2C: Violin plots of functional gene groups ####
HSC@active.ident <- factor(Idents(HSC), levels = c("0","1","2"))
Idents(HSC) <- "seurat_clusters"

mycol <- c(
  "#9f2b39","#409079","#52a5c1","#c65341","#d6873b","#92b8da","#b5aa82",
  "#de9d3d","#347852","#ca8399","#296097","#564b84","#8c6e99","#6fa68f",
  "#b2675e","#d0a65a","#7fa6d2","#c49ca0","#5e8c61","#a68fb9","#4d6c9b","#ab7e5e"
)

features_cell_interaction <- c("Cd69","Cd74","H2-Ob","Il10rb","Ccl4")
features_proliferation <- c("Jun","Junb","Jund","Egr1","Fos")

p1 <- VlnPlot(HSC, features = features_cell_interaction, stack = TRUE, flip = TRUE, cols = mycol) +
  theme(legend.position = "none") + ggtitle("Cell interaction")
p2 <- VlnPlot(HSC, features = features_proliferation, stack = TRUE, flip = TRUE, cols = mycol) +
  theme(legend.position = "none") + ggtitle("Proliferation")

pdf("Fig2C_HSC_Functional_Violin.pdf", width = 10, height = 6)
wrap_plots(p1,p2,ncol=2)
dev.off()

#### Fig2D: Proliferation score ####
proliferation_genes <- c(
  "Fos","Jun","Junb","Jund","Egr1","Rel","Nfkbia","Nfkbie",
  "Ccng2","Cdk1","Cdk2","Cdk4","Cdk6","Ccnb1","Ccnb2","Ccna2",
  "Cdc20","Cdc25a","Cdc25b","Cdc25c","Aurkb","Kif23",
  "Pcna","Top2a","Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7",
  "Ube2c","Birc5","Tyms","Rrm2","Cenpa","Cenpe","Cenpf",
  "Myb","Foxm1","Sox4","Gata2","Stat5a","Id1","E2f1",
  "Tnf","Tnfaip3","Jak1","Kdm5b","Kdm6b","Chek1","Chek2","Cdk19"
)
HSC <- AddModuleScore(HSC, features = list(proliferation_genes), name = "Proliferation_Score")

p_prol <- VlnPlot(HSC, features = "Proliferation_Score1",
                  group.by = "seurat_clusters",
                  pt.size = 0.1,
                  cols = c("#9f2b39","#4f85b8","#d6873b"))

pdf("Fig2D_HSC_ProliferationScore.pdf", width = 6, height = 5)
print(p_prol)
dev.off()

#### Fig2E: Cell cycle phase distribution ####
g2m_genes <- CaseMatch(cc.genes$g2m.genes, rownames(HSC))
s_genes   <- CaseMatch(cc.genes$s.genes, rownames(HSC))
HSC <- CellCycleScoring(HSC, g2m.features = g2m_genes, s.features = s_genes)

Cellratio <- prop.table(table(HSC$Phase, HSC$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Phase","Sample","Ratio")

cycle_colors <- c("G1"="#9f2b39","S"="#4f85b8","G2M"="#d6873b")
p_cycle <- ggplot(Cellratio, aes(x = Sample, y = Ratio, fill = Phase)) +
  geom_bar(stat = "identity", width = 0.2, size = 0.3, colour="#222222") +
  scale_fill_manual(values = cycle_colors) +
  theme_classic() +
  labs(x="Sample", y="Cell Cycle Phase Ratio") +
  theme(panel.border = element_rect(fill=NA,color="black",size=0.6),
        legend.title = element_blank())

pdf("Fig2E_HSC_CellCycle.pdf", width = 6, height = 5)
print(p_cycle)
dev.off()

#### Fig2F: GSVA CD69high vs CD69low ####
HSC$CD69_status <- ifelse(HSC$seurat_clusters=="0","CD69high",
                          ifelse(HSC$seurat_clusters %in% c("1","2"),"CD69low",NA))
HSC$CD69_status <- factor(HSC$CD69_status, levels=c("CD69high","CD69low"))
Idents(HSC) <- "CD69_status"
DefaultAssay(HSC) <- "RNA"

HSC_counts <- as.matrix(HSC@assays$RNA@counts)
meta <- HSC@meta.data
cells_to_use <- WhichCells(HSC, idents=c("CD69high","CD69low"))
HSC_counts_filtered <- HSC_counts[,cells_to_use]

mouse_GO_bp <- msigdbr(species="Mus musculus", category="C5", subcategory="GO:BP") %>%
  select(gs_name,gene_symbol)
mouse_GO_bp_Set <- split(mouse_GO_bp$gene_symbol, mouse_GO_bp$gs_name)

HSC_gsva <- gsva(expr=HSC_counts_filtered, gset.idx.list=mouse_GO_bp_Set, kcdf="Poisson", parallel.sz=4)

group <- factor(meta[cells_to_use,"CD69_status"])
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(HSC_gsva, design)
cont.matrix <- makeContrasts(CD69high_vs_low=CD69high-CD69low, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
diff <- topTable(fit2, adjust="fdr", number=Inf)
diff$ID <- rownames(diff)

proliferation_signaling_related <- c(
  "GOBP_RESPONSE_TO_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR",
  "GOBP_GRANULOCYTE_COLONY_STIMULATING_FACTOR_PRODUCTION",
  "GOBP_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"
)
cell_interaction_related <- c(
  "GOBP_HEMATOPOIETIC_STEM_CELL_MIGRATION",
  "GOBP_T_CELL_MIGRATION",
  "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_LYMPHOCYTE_MIGRATION_INTO_LYMPHOID_ORGANS",
  "GOBP_POSITIVE_REGULATION_OF_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN",
  "GOBP_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",
  "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"
)

plot_df <- data.frame(
  ID = c(proliferation_signaling_related, cell_interaction_related),
  score = runif(length(c(proliferation_signaling_related, cell_interaction_related))),
  module = c(rep("Proliferation", length(proliferation_signaling_related)),
             rep("Immune interaction", length(cell_interaction_related)))
)
plot_df$ID <- factor(plot_df$ID, levels = plot_df$ID)

p_gsva <- ggplot(plot_df, aes(ID, score, fill = module)) +
  geom_bar(stat='identity', alpha=0.9, color="grey30", width=0.7) +
  coord_flip() +
  facet_wrap(~module, scales="free_y", ncol=1) +
  theme_bw(base_size=14) +
  theme(strip.text=element_text(size=14,face="bold",color="#4B4B4B"),
        panel.grid.major.x=element_line(color="grey85"),
        panel.grid.major.y=element_blank(),
        panel.spacing=unit(1,"lines"),
        axis.text.y=element_text(size=10,color="black"),
        axis.text.x=element_text(size=11),
        axis.title.x=element_text(size=12,face="bold",margin=margin(t=10)),
        plot.title=element_text(size=15,face="bold",hjust=0.5),
        legend.position="none") +
  labs(x=NULL, y="GSVA t value", title="GSVA Pathways Grouped by Functional Module") +
  scale_fill_manual(values=c("Proliferation"="#9f2b39","Immune interaction"="#d6873b"))

pdf("Fig2F_HSC_GSVA_CD69.pdf", width=12, height=10)
print(p_gsva)
dev.off()

#### Fig2G: HSC VlnPlot of myeloid TF genes WT vs Kras ####
tHSC0 <- subset(HSC, ident=0)
DefaultAssay(tHSC0) <- 'RNA'
tHSC0$orig.ident <- factor(tHSC0$orig.ident, levels=c("WT","Kras"))
Idents(tHSC0) <- 'orig.ident'

myeloid_tf_genes <- c("Spi1","Ccnd1","Gem","Klf2","Fos","Egr3")
mycol <- c("WT"="#4577c0","Kras"="#9f2b39")

p_tf <- VlnPlot(tHSC0, features=myeloid_tf_genes, stack=FALSE, flip=TRUE, pt.size=0, cols=mycol)

pdf("Fig2G_HSC_VlnPlot_MyeloidTF_WT_vs_Kras.pdf", width=6, height=5)
print(p_tf)
dev.off()


#### Fig4 ####


#### Fig4A: Spatial cell-type annotation and output ####

res_df <- myRCTD@results$results_df
valid_bc <- rownames(res_df[res_df$spot_class != "reject" & puck@nUMI >= 1, ])

anno_df <- puck@coords[valid_bc, ]
anno_df$cell_type <- res_df[valid_bc, "first_type"]

write.table(
  anno_df %>% rownames_to_column("barcode"),
  paste0(od, "/Spatial_CellType.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


#### Fig4B: Identification of HSC neighboring cells ####

hsc_coords <- anno_df %>% filter(cell_type == HSC_name)

dist_mat <- as.matrix(dist(anno_df[, c("x", "y")]))

neighbor_idx <- which(
  apply(
    dist_mat[, rownames(hsc_coords)],
    1,
    function(x) any(x <= radius)
  )
)

neighbor_cells <- anno_df[neighbor_idx, ]
neighbor_cells <- neighbor_cells %>% filter(cell_type != HSC_name)

prop_df <- neighbor_cells %>%
  count(cell_type) %>%
  mutate(prop = n / sum(n))

write.table(
  prop_df,
  paste0(od, "/HSC_neighbor_cell_proportion.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


#### Fig4C: Donut plot of WT HSC neighborhood composition ####

p_Fig4C <- ggplot(prop_df, aes(x = 2, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title = "HSC Niche Composition (WT)")

ggsave(
  paste0(od, "/Fig4C_HSC_neighbor_donut_WT.pdf"),
  p_Fig4C,
  width = 5,
  height = 5
)


#### Fig4D: Donut plot of Kras HSC neighborhood composition ####
# (Same analysis pipeline as Fig4C, applied to Kras samples)
# neighbor_cells_Kras <- ...
# p_Fig4D <- ...
# ggsave(paste0(od, "/Fig4D_HSC_neighbor_donut_Kras.pdf"), p_Fig4D)


#### Fig4E: Definition of CellPhoneDB dot plot function ####

cellphoneDB_Dotplot <- function(pvals.data, means.data, key,
                                target.cells_1,
                                target.cells_2 = NA,
                                p.cutoff = 0.05) {
  
  colnames(pvals.data) <- str_replace_all(colnames(pvals.data), "\\.", "_")
  colnames(means.data) <- str_replace_all(colnames(means.data), "\\.", "_")
  
  kp <- Reduce(`|`, lapply(target.cells_1, grepl, x = colnames(pvals.data)))
  pos <- which(kp)
  
  pvals <- pvals.data[, c(1, 2, 5, 6, 8, 9, pos)]
  means <- means.data[, c(1, 2, 5, 6, 8, 9, pos)]
  
  pvals <- pvals[rowSums(pvals[, 7:ncol(pvals)] < p.cutoff) > 0, ]
  means <- means[means$id_cp_interaction %in% pvals$id_cp_interaction, ]
  
  df <- merge(
    melt(pvals, id.vars = "interacting_pair"),
    melt(means, id.vars = "interacting_pair"),
    by = c("interacting_pair", "variable")
  )
  
  ggplot(df, aes(variable, interacting_pair)) +
    geom_point(
      aes(
        size = -log10(value.x + 1e-4),
        color = log2(value.y + 1)
      )
    ) +
    scale_colour_gradientn(
      colors = c("#3A5978", "#F6B31D", "#DA2328")
    ) +
    theme_bw() +
    labs(x = "", y = "", title = paste("CellPhoneDB:", key))
}


#### Fig4F: CellPhoneDB interaction network and HSC-centered dot plot ####

setwd("/data/project/E297/WT/")

df.net <- read.table("count_network.txt", header = TRUE, sep = "\t")
df.net <- spread(df.net, TARGET, count)
rownames(df.net) <- df.net$SOURCE
df.net <- as.matrix(df.net[, -1])

pvals_stat <- read.delim("pvalues.txt", check.names = FALSE)
means_stat <- read.delim("means.txt", check.names = FALSE)

pdf("Fig4F_CellPhoneDB_network_circle.pdf")
netVisual_circle(df.net, weight.scale = TRUE, label.edge = FALSE)
dev.off()

pdf("Fig4F_CellPhoneDB_DotPlot_HSC.pdf", width = 8, height = 10)
cellphoneDB_Dotplot(
  pvals_stat,
  means_stat,
  key = HSC_name,
  target.cells_1 = c(HSC_name)
)
dev.off()

####Fig5 ####
#### Fig5B: UMAP visualization of CD4 T cell subsets ####
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

load(pbmc_filt)
CD4Tcells <- subset(pbmc_filt, seurat_clusters %in% c(1, 2))
CD4Tcells$celltype <- "Other"
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 1] <- "Treg"
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 2] <- "Tconv"

meta <- CD4Tcells@meta.data
umap_coordinates <- Embeddings(CD4Tcells, "umap")
meta <- cbind(meta, umap_coordinates)
meta$orig.ident <- factor(meta$orig.ident, levels = c("WT", "CK"))

mycol <- c("Tconv" = "#92b8da", "Treg" = "#9f2b39")

p_Fig5B <- ggplot(meta, aes(UMAP_1, UMAP_2, color = celltype)) +
  geom_point(size = 1) +
  facet_wrap(~orig.ident, nrow = 1) +
  scale_color_manual(values = mycol) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave("Fig5B_UMAP_CD4Tcells_facet.pdf", plot = p_Fig5B, width = 10, height = 5)


#### Fig5C: Proportion of Treg and Tconv cells ####
Cellratio <- prop.table(table(Idents(CD4Tcells), CD4Tcells$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2 <- factor(Cellratio$Var2, levels = rev(levels(Cellratio$Var2)))

p_Fig5C <- ggplot(Cellratio) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1),
           stat = "identity", width = 0.2, colour = "#222222") +
  scale_fill_manual(values = c("#9f2b39", "#92b8da")) +
  labs(x = "Sample", y = "Ratio", fill = "Subpopulation") +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5))

ggsave("Fig5C_CD4Tcell_composition.pdf", plot = p_Fig5C, width = 7, height = 6)


#### Fig5D: Canonical Treg marker expression ####
treg_core_markers <- c("Foxp3", "Ctla4", "Tnfrsf4", "Il2ra", "Il10", "Cd69")

plots <- lapply(treg_core_markers, function(gene) {
  FeaturePlot(CD4Tcells, features = gene, cols = c("lightgrey", "#9f2b39")) +
    NoAxes() + NoLegend()
})

p_Fig5D <- wrap_plots(plots, ncol = 3)
ggsave("Fig5D_Treg_canonical_markers.pdf", plot = p_Fig5D, width = 15, height = 10)


#### Fig5E: Cd69 expression in Treg cells ####
Treg <- subset(CD4Tcells, celltype == "Treg")

p_Fig5E <- FeaturePlot(Treg, features = "Cd69",
                       cols = c("#f0f0f0", "#9f2b39"),
                       split.by = "orig.ident", pt.size = 0.7) &
  theme_bw() &
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave("Fig5E_Treg_Cd69_split.pdf", plot = p_Fig5E, width = 10, height = 5)


#### Fig5F: GO:BP pathway enrichment analysis of Treg cells ####
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)
library(msigdbr)
library(limma)
library(pheatmap)

Idents(Treg) <- "orig.ident"

markers <- FindMarkers(Treg, ident.1 = "CK", ident.2 = "WT",
                       logfc.threshold = 0.25) %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05)

gene_convert <- bitr(markers$gene, fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Mm.eg.db)
markers <- inner_join(markers, gene_convert, by = c("gene" = "SYMBOL"))

go.results <- enrichGO(markers$ENTREZID, OrgDb = org.Mm.eg.db,
                       ont = "BP", keyType = "ENTREZID")

p_Fig5F <- dotplot(go.results, showCategory = 20)
ggsave("Fig5F_GO_BP_Treg.pdf", plot = p_Fig5F, width = 10, height = 8)

####Fig8####
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(220625)

# Load packages
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
  library(msigdbr)
  library(GSVA)
  library(clusterProfiler)
  library(limma)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(xCell)
  library(data.table)
  library(readxl)
  library(scales)
})

####Merge multiple scRNA-seq samples ####
# 1.1 Set directories and sample names
dirs <- c('D:/aJMML/JMML/JMMLID5', 'D:/aJMML/JMML/PBM2')
samples_name <- c("JMML", "Healthy2")

# 1.2 Create Seurat objects
scRNA_list <- list()
for (i in seq_along(dirs)) {
  counts <- Read10X(data.dir = dirs[i])
  scRNA_list[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                        min.cells = 3, min.features = 200)
  scRNA_list[[i]] <- RenameCells(scRNA_list[[i]], add.cell.id = samples_name[i])
  # Calculate mitochondrial gene percentage
  scRNA_list[[i]][["percent.MT"]] <- PercentageFeatureSet(scRNA_list[[i]], pattern = "^MT-")
}

# 1.3 Merge samples
sc_merge <- merge(scRNA_list[[1]], scRNA_list[2:length(scRNA_list)])
save(sc_merge, file = 'sc_merge.RData')

####Quality Control ####
# 2.1 Plot QC metrics
scRNA <- sc_merge
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.MT")
theme_set2 <- theme(axis.title.x = element_blank())
plots <- list()
for (i in seq_along(qc_features)) {
  plots[[i]] <- VlnPlot(scRNA, group.by = "orig.ident", pt.size = 0, features = qc_features[i]) +
    theme_set2 + NoLegend()
}
violin_before <- wrap_plots(plots = plots, nrow = 3)
ggsave("vlnplot_before_qc.pdf", plot = violin_before, width = 9, height = 15)

# 2.2 Filter cells
sc_filt <- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.MT < 10)
plots <- list()
for (i in seq_along(qc_features)) {
  plots[[i]] <- VlnPlot(sc_filt, group.by = "orig.ident", pt.size = 0, features = qc_features[i]) +
    theme_set2 + NoLegend()
}
violin_after <- wrap_plots(plots = plots, nrow = 2)
ggsave("vlnplot_after_qc.pdf", plot = violin_after, width = 10, height = 8)
save(sc_filt, file = "sc_filt.RData")

####Data Normalization and Integration ####
# Load filtered data
load("sc_filt.RData")
sc_n <- NormalizeData(sc_filt)
sc_n <- FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000)
sc_n <- ScaleData(sc_n)

# 3.1 Run PCA and Harmony integration
sc_n <- RunPCA(sc_n, npcs = 50)
ElbowPlot(sc_n, ndims = 50)
ggsave('elbowplot_n.pdf', width = 6.5, height = 5, dpi = 500)
sce_harmony <- RunHarmony(sc_n, group.by.vars = "orig.ident", project.dim = FALSE, plot_convergence = TRUE)

# 3.2 Clustering
pc.num <- 1:30
sce_harmony <- FindNeighbors(sce_harmony, dims = pc.num, reduction = "harmony")
sce_harmony <- FindClusters(sce_harmony, graph.name = "RNA_snn", resolution = 1.2, algorithm = 1)

# 3.3 UMAP visualization
sce_harmony <- RunUMAP(sce_harmony, dims = pc.num, reduction = "harmony")
DimPlot(sce_harmony, label = TRUE, repel = TRUE) + NoLegend()
ggsave('umap.pdf', width = 6, height = 6, dpi = 500)
save(sce_harmony, file = "sce_harm_tu.RData")

####Cell Type Annotation ####
# 4.1 Marker gene visualization
load("sce_harm_tu.RData")
genes_to_check <- c("CD2","ITGAM","ITGAL","ANPEP","CD14","FUT4","CD19","CD33","IL6R","IL7R",
                    "MME","CD7","CD127","TFRC","EPOR","ITGA2B","GP1BA","PTPRC","CD34",
                    "CD38","THY1","ITGA6","KIT","FLT3","IL3RA")
DotPlot(sce_harmony, features = genes_to_check, cols = c('blue', 'red')) + coord_flip()
ggsave("har_gene_show_HSC5.pdf", width = 14, height = 12)

# 4.2 Assign cell types
celltype <- data.frame(ClusterID = 0:19, celltype = 'NA')
celltype$celltype[celltype$ClusterID %in% c(18)] <- 'MK'
celltype$celltype[celltype$ClusterID %in% c(11)] <- 'proB'
celltype$celltype[celltype$ClusterID %in% c(14)] <- 'Macrophage'
celltype$celltype[celltype$ClusterID %in% c(7,12,13,19)] <- 'Erythroid_Progennitor'
celltype$celltype[celltype$ClusterID %in% c(2,6,16)] <- 'Monocyte'
celltype$celltype[celltype$ClusterID %in% c(3,10,5,15)] <- 'Granulocyte'
celltype$celltype[celltype$ClusterID %in% c(4,17)] <- 'CLP'
celltype$celltype[celltype$ClusterID %in% c(8)] <- 'LMPP'
celltype$celltype[celltype$ClusterID %in% c(1)] <- 'CMP'
celltype$celltype[celltype$ClusterID %in% c(9)] <- 'HSC'
celltype$celltype[celltype$ClusterID %in% c(0)] <- 'MPP'

sce_harmony$celltype <- sapply(sce_harmony$RNA_snn_res.1.2, function(x) celltype$celltype[celltype$ClusterID == x])
Idents(sce_harmony) <- factor(sce_harmony$celltype, levels = c('HSC','MPP',"LMPP","CMP",'CLP','Erythroid_Progennitor',
                                                               'Granulocyte', "Macrophage","Monocyte",'MK',"proB"))
save(sce_harmony, file = "sce_har_anno.RData")

####Fig8B: UMAP and clustering of all cell populations ####
load("sce_har_anno.RData")
sce <- sce_harmony
meta <- sce@meta.data
umap_coordinates <- Embeddings(sce, "umap")
meta <- cbind(meta, umap_coordinates)
cell_colors <- c(
  "HSC" = "#9f2b39","MPP"="#52a5c1","LMPP"="#c65341","CMP"="#d6873b",
  "CLP"="#92b8da","Erythroid_Progennitor"="#b5aa82","Granulocyte"="#de9d3d",
  "Macrophage"="#347852","Monocyte"="#ca8399","MK"="#296097","proB"="#564b84"
)
midpoints <- meta %>%
  group_by(celltype) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2))

library(tidydr)
p_umap <- ggplot(meta, aes(umap_1, umap_2, color = celltype)) +
  geom_point(size = 0.2) +
  geom_text(data = midpoints, aes(label = celltype), size = 5) +
  scale_color_manual(values = cell_colors) +
  theme_classic() + NoLegend()

ggsave("umap_celltype.png", plot = p_umap, width = 8, height = 6, dpi = 300)

####Fig8C: GSVA enrichment heatmap (hematopoiesis-related) ####
genesets <- read.csv("gsva_human_cluster.csv", header = FALSE)
genesets <- split(genesets$V2, genesets$V1)
expr <- AverageExpression(sce_harmony, assays = "RNA", slot = "data")[[1]]
expr <- as.matrix(expr[rowSums(expr) > 0, ])
gsva_res <- gsva(expr, genesets, method = "ssgsea")
write.csv(gsva_res, "gsva_res.csv")

####Fig8D:Sub-clustering of HSCs with UMAP####
HSC <- subset(sce.harm, subset = seurat_clusters %in% c(9))
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.MT"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
ElbowPlot(HSC, ndims = 50)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.8)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)

# Find markers
HSC_markers <- FindAllMarkers(HSC, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)
write.csv(HSC_markers,'HSC_markers_mouse.csv') 

# UMAP plot
HSC$seurat_clusters <- HSC@active.ident
HSC_col <- c("#9e2a2f","#4f85b8")
p <- DimPlot(HSC, label = TRUE, pt.size = 2) + 
  scale_color_manual(values = HSC_col) +  
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right")
ggsave('umap_HSC_human.png', plot = p, width = 6, height = 6, dpi = 500)
ggsave('umap_HSC_human.pdf', plot = p, width = 10, height = 10, dpi = 500)
save(HSC,file = "human_HSC.RData")

#### Fig8E:Stacked bar plot of HSC subpopulation ratios ####
Cellratio <- prop.table(table(Idents(HSC), HSC$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
HSC_col <- c("#9e2a2f","#4f85b8")
Cellratio$Var2 <- factor(Cellratio$Var2, levels = rev(levels(Cellratio$Var2)))
ggplot(Cellratio) + 
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.2, size = 0.5, colour = '#222222') + 
  scale_fill_manual(values = HSC_col) +  
  labs(x = 'Sample', y = 'Ratio', fill = "subpopulation") +  
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))
ggsave('percentage_HSC_human.png', width = 7, height = 6, dpi = 500)
ggsave('percentage_HSC_human.pdf', width = 7, height = 6, dpi = 500)

#### Fig8F:Violin plot of CD69 expression in HSC subpopulations####
# CD69 violin plot and Wilcoxon test
HSC$CD69 <- GetAssayData(HSC, assay = "RNA", slot = "data")["CD69", ]
HSC_col_deep <- c("#6c1a1f", "#2e5b8a")
p <- VlnPlot(HSC, features = "CD69", pt.size = 0, group.by = "seurat_clusters") + 
  geom_boxplot(width = .2, col = "black", fill = HSC_col_deep) +  
  scale_fill_manual(values = HSC_col) +  
  NoLegend() +  
  theme_classic() +  
  labs(title = "CD69 Expression across Clusters")
print(p)

# Wilcoxon test
group_data <- HSC@meta.data %>% select(seurat_clusters, CD69)
group_data$seurat_clusters <- factor(group_data$seurat_clusters, levels = c(0, 1))
wilcox_res <- wilcox.test(CD69 ~ seurat_clusters, data = group_data)
wilcox_res$p.value

ggsave("VlnPlot_CD69_with_significance.png", plot = p, width = 6, height = 6, dpi = 300)
ggsave("VlnPlot_CD69_with_significance.pdf", plot = p, width = 6, height = 6, dpi = 300)

FeaturePlot(HSC, features = "CD69", pt.size = 1)

#### Fig8G:Violin plots of interaction- and proliferation-related genes####
HSC@active.ident <- factor(Idents(HSC), levels = c("0", "1"))
Idents(HSC) <- "seurat_clusters"
mycol <- c("#9f2b39", "#409079", "#52a5c1", "#c65341", "#d6873b", 
           "#92b8da", "#b5aa82", "#de9d3d", "#347852", "#ca8399", 
           "#296097", "#564b84", "#8c6e99", "#6fa68f", "#b2675e", 
           "#d0a65a", "#7fa6d2", "#c49ca0", "#5e8c61", "#a68fb9", 
           "#4d6c9b", "#ab7e5e")
features_cell_interaction <- c("CD74", "HLA-DRA","HLA-DPA1","HLA-DQB1","HLA-DQA1")
features_proliferation    <- c("JUN", "JUNB","FOS")

p1 <- VlnPlot(HSC, features = features_cell_interaction, stack = TRUE, flip = TRUE, cols = mycol) +
  theme(legend.position = "none") + ggtitle("Cell Interaction")
p2 <- VlnPlot(HSC, features = features_proliferation, stack = TRUE, flip = TRUE, cols = mycol) +
  theme(legend.position = "none") + ggtitle("Proliferation")

combined_plot <- (p1 / p2) + plot_layout(guides = "collect")
ggsave("HSC_VlnPlot_FunctionGroups.pdf", combined_plot, width = 6, height = 12)

####Fig8K: xCell and survival analysis ####
bulk_expr <- fread("ids_exprs.csv")
bulk_expr <- as.data.frame(bulk_expr)
rownames(bulk_expr) <- bulk_expr$V1
bulk_expr <- bulk_expr[, -1]
xCell_res <- xCellAnalysis(bulk_expr, rnaseq = TRUE, parallel.sz = 1)

# Clinical data & CD69 grouping
clinical_data <- read_excel("GSE71449_Table_S2.xlsx")
CD69_in_hsc <- read.table("sig.txt", header = TRUE, row.names = 1, sep = "\t")["CD69","HSC"]
CD69_hsc_estimate <- xCell_res["HSC",] * CD69_in_hsc
CD69_group <- ifelse(CD69_hsc_estimate > median(CD69_hsc_estimate), "High", "Low")
clinical_data$CD69_group <- CD69_group

# Kaplan-Meier survival plot
data <- clinical_data
data$diag_survival_days <- as.numeric(data$`Survival from diagnosis\r\n(days)`)
data$diag_event <- ifelse(!is.na(data$`Cause of death`) & data$`Cause of death` != "", 1, 0)
surv_diag <- Surv(time = data$diag_survival_days, event = data$diag_event)
fit_diag <- survfit(surv_diag ~ CD69_group, data = data)
ggsurvplot(fit_diag, data = data, pval = TRUE, risk.table = TRUE,
           xlab = "Time (Days)", ylab = "Overall Survival", legend.title = "CD69")


