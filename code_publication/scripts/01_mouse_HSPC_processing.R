# ============================================================================
# Script: 01_mouse_HSPC_processing.R
# Purpose: Mouse bone-marrow HSPC single-cell RNA-seq processing and analysis:
#          read 10x matrices (Kras and WT), QC filtering, Harmony integration,
#          clustering (Louvain, resolution 1.2), UMAP, cluster marker detection,
#          marker-based cell-type annotation, GSVA pathway scoring, and
#          cell-type proportion comparison between genotypes.
# Inputs:  data/raw/mouse_HSPC/Kras/ and data/raw/mouse_HSPC/WT/
#          (10x Cell Ranger output directories; barcodes/features/matrix)
#          data/gene_sets/gsva_mouse_cluster.csv  (curated module gene sets)
# Outputs: output/mouse_HSPC_annotated.rds  (Seurat object with $celltype)
#          output/gsva_res_mouse.csv, figures/umap_celltype.pdf, etc.
# Run order: after 00_setup.R (run this script first)
# Prerequisites: raw 10x output directories present (see data/README.md)
# NOTE: The final cell-type annotation is marker-based (manual cluster
#       assignment). An optional cross-validation by label transfer from
#       GSE122465 is provided separately in 02_mouse_annotation_crossvalidation.R.
# ============================================================================

source("scripts/00_setup.R")

# --- 1. read and merge samples ---------------------------------------------
dirs <- c(file.path(raw_dir, "mouse_HSPC", "Kras"),
          file.path(raw_dir, "mouse_HSPC", "WT"))
samples_name <- c("Kras", "WT")

scRNAlist <- list()
for (i in seq_along(dirs)) {
  counts <- Read10X(data.dir = dirs[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                       min.cells = 3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
}
sc_merge <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

# --- 2. QC filtering --------------------------------------------------------
# Thresholds from the original analysis: 500-6,000 detected genes, <10% mt reads
sc_filt <- subset(sc_merge,
                  subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# --- 3. normalization, variable features, scaling --------------------------
sc_n <- NormalizeData(sc_filt)
sc_n <- FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sc_n <- ScaleData(sc_n)

# --- 4. PCA, Harmony integration, clustering, UMAP -------------------------
sc_n <- RunPCA(sc_n, npcs = 50, verbose = FALSE)
pc.num <- 1:30
sce_har <- RunHarmony(sc_n, group.by.vars = "orig.ident", project.dim = FALSE,
                      plot_convergence = FALSE)
sce.har <- FindNeighbors(sce_har, dims = pc.num, reduction = "harmony")
sce.har <- FindClusters(sce.har, graph.name = "RNA_snn", resolution = 1.2, algorithm = 1)
sce.harm <- RunUMAP(sce.har, dims = pc.num, reduction = "harmony")

# --- 5. cluster markers -----------------------------------------------------
# NOTE: meta column referenced for manual annotation below is RNA_snn_res.0.6;
#       confirm the correspondence between the stored clustering column and
#       the resolution 1.2 clustering used for the UMAP.
allmarkers <- FindAllMarkers(sce.harm, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)
write.csv(allmarkers, file.path(output_dir, "allmarkers_celltype_mouse.csv"))

# --- 6. marker-based cell-type annotation -----------------------------------
# Manual cluster->celltype mapping from the original analysis.
celltype <- data.frame(ClusterID = 0:22, celltype = NA_character_)
celltype[celltype$ClusterID %in% c(0,2,5,8,9,13,16), "celltype"] <- "Granulocyte"
celltype[celltype$ClusterID %in% 14, "celltype"] <- "Macrophage"
celltype[celltype$ClusterID %in% 18, "celltype"] <- "Monocyte"
celltype[celltype$ClusterID %in% c(11,22), "celltype"] <- "B"
celltype[celltype$ClusterID %in% c(7,17), "celltype"] <- "cDC"
celltype[celltype$ClusterID %in% 21, "celltype"] <- "Ery"
celltype[celltype$ClusterID %in% c(10,12), "celltype"] <- "Erythroblast"
celltype[celltype$ClusterID %in% 15, "celltype"] <- "HSC"
celltype[celltype$ClusterID %in% 6, "celltype"] <- "MPP"
celltype[celltype$ClusterID %in% 19, "celltype"] <- "T"
celltype[celltype$ClusterID %in% c(1,3,4), "celltype"] <- "GMP"
celltype[celltype$ClusterID %in% 20, "celltype"] <- "MK"

sce.harm$celltype <- NA_character_
for (i in seq_len(nrow(celltype))) {
  sce.harm@meta.data[which(sce.harm@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),
                     "celltype"] <- celltype$celltype[i]
}
# NOTE: the original code referenced the meta column RNA_snn_res.0.6; if the
#       object stores a different clustering column, adjust the column name.

celltype_order <- c("HSC","MPP","GMP","Erythroblast","Ery","Granulocyte",
                    "Macrophage","Monocyte","MK","B","T","cDC")
Idents(sce.harm) <- factor(sce.harm$celltype, levels = celltype_order)
sce.harm$celltype <- Idents(sce.harm)

saveRDS(sce.harm, file.path(output_dir, "mouse_HSPC_annotated.rds"))

# --- 7. GSVA pathway scoring (curated modules) ------------------------------
# Module gene sets are provided as a two-column CSV (V1 = module, V2 = gene).
genesets <- read.csv(file.path(gene_dir, "gsva_mouse_cluster.csv"), header = FALSE)
genesets <- split(genesets$V2, genesets$V1)

expr <- AverageExpression(sce.harm, assays = "RNA", layer = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ]
expr <- as.matrix(expr)

if (packageVersion("GSVA") >= "1.50.0") {
  param <- gsvaParam(expr, genesets, kcdf = "Gaussian")
  gsva.res <- gsva(param, verbose = FALSE)
} else {
  gsva.res <- gsva(expr, genesets, method = "ssgsea")
}
write.csv(data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = FALSE),
          file.path(output_dir, "gsva_res_mouse.csv"), row.names = FALSE)

# --- 8. cell-type proportions by genotype -----------------------------------
Cellratio <- prop.table(table(Idents(sce.harm), sce.harm$orig.ident), margin = 2)
Cellratio_df <- as.data.frame(Cellratio)
colnames(Cellratio_df) <- c("celltype", "sample", "freq")
write.csv(Cellratio_df, file.path(output_dir, "mouse_celltype_proportions.csv"), row.names = FALSE)

p_prop <- ggplot(Cellratio_df, aes(x = sample, y = freq, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.2, colour = "#222222") +
  theme_classic() +
  labs(x = "Sample", y = "Cell proportion") +
  theme(panel.border = element_rect(fill = NA, colour = "black"))
ggsave(file.path(fig_dir, "mouse_celltype_proportion.pdf"), p_prop, width = 7, height = 6)