# ============================================================================
# Script: 04_mouse_Tcell_analysis.R
# Purpose: Mouse bone-marrow CD4+ T-cell analysis: sub-clustering of CD4+
#          T cells, Treg/Tconv identification, differential expression between
#          genotypes, GO enrichment, and GSVA pathway activity of Tregs.
# Inputs:  data/raw/mouse_Tcell/            (10x Cell Ranger output, CD45+CD3+
#          sorted bone-marrow T cells; three mice pooled per genotype)
# Outputs: output/mouse_CD4Tcells.rds, output/Treg_DE.csv,
#          output/Treg_GO_BP.csv, output/Treg_GSVA_descriptive.csv
# Run order: after 00_setup.R
# NOTE 1: Three mice per genotype were pooled into one 10x library per
#         genotype; the WT vs Kras comparison is therefore based on one
#         library per genotype (exploratory, no library-level replication).
# NOTE 2: QC/integration of the T-cell object used the same parameters as the
#         HSPC pipeline (per manuscript). If the raw 10x output is available
#         in data/raw/mouse_Tcell/, it is rebuilt below; otherwise load a
#         previously prepared CD4+ T-cell object from output/.
# ============================================================================

source("scripts/00_setup.R")

# --- 1. (re)build the CD4+ T-cell object from raw data ----------------------
# If the pre-processed object already exists, load it and skip the build.
tcell_obj_file <- file.path(output_dir, "mouse_Tcell_CD4.rds")
if (file.exists(tcell_obj_file)) {
  pbmc_filt <- readRDS(tcell_obj_file)
} else {
  # NOTE: the original T-cell QC/integration was performed outside the shared
  #       scripts; the block below reproduces it with the HSPC parameters.
    samples_name <- c("Kras", "WT")
  dirs <- c(file.path(raw_dir, "mouse_Tcell", "Kras"),
            file.path(raw_dir, "mouse_Tcell", "WT"))
  scRNAlist <- list()
  for (i in seq_along(dirs)) {
    counts <- Read10X(data.dir = dirs[i])
    scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                         min.cells = 3, min.features = 200)
    scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-")
  }
  sc_merge <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
  sc_filt <- subset(sc_merge,
                    subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)
  sc_n <- NormalizeData(sc_filt)
  sc_n <- FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sc_n <- ScaleData(sc_n)
  sc_n <- RunPCA(sc_n, npcs = 50, verbose = FALSE)
  sc_n <- RunUMAP(sc_n, dims = 1:30, reduction = "pca",
                  reduction.name = "umap_pca", reduction.key = "UMAPpca_")
  p_no_harmony <- DimPlot(sc_n, reduction = "umap_pca", group.by = "orig.ident") +
    ggtitle("Mouse T cells before Harmony (one pooled library per genotype)")
  ggsave(file.path(fig_dir, "sensitivity_mouse_Tcell_preHarmony.pdf"),
         p_no_harmony, width = 7, height = 6)
  sce_har <- RunHarmony(sc_n, group.by.vars = "orig.ident")
  pbmc_filt <- FindNeighbors(sce_har, dims = 1:30, reduction = "harmony")
  pbmc_filt <- FindClusters(pbmc_filt, resolution = 0.6, algorithm = 1)
  pbmc_filt <- RunUMAP(pbmc_filt, dims = 1:30, reduction = "harmony")
  pbmc_filt$orig.ident <- factor(pbmc_filt$orig.ident, levels = c("WT", "Kras"))
  saveRDS(pbmc_filt, tcell_obj_file)
}
check_group_counts(pbmc_filt, "orig.ident", c(Kras = 233L, WT = 759L))

# --- 2. CD4+ T-cell subsetting and Treg/Tconv assignment --------------------
# The CD4+ T-cell compartment corresponds to clusters 1 and 2 of the
# CD45+CD3+ object (from the original analysis).
CD4Tcells <- subset(pbmc_filt, seurat_clusters %in% c(1, 2))
CD4Tcells$celltype <- "Other"
CD4Tcells$celltype[as.character(CD4Tcells$seurat_clusters) == "1"] <- "Treg"
CD4Tcells$celltype[as.character(CD4Tcells$seurat_clusters) == "2"] <- "Tconv"
Idents(CD4Tcells) <- "celltype"
check_group_counts(CD4Tcells, "orig.ident", c(Kras = 190L, WT = 261L))
saveRDS(CD4Tcells, file.path(output_dir, "mouse_CD4Tcells.rds"))

# --- 3. differential expression: Tregs, Kras vs WT --------------------------
Treg <- subset(CD4Tcells, celltype == "Treg")
Idents(Treg) <- "orig.ident"
markers <- FindMarkers(Treg, ident.1 = "Kras", ident.2 = "WT",
                       only.pos = FALSE, logfc.threshold = 0.25)
markers$gene <- rownames(markers)
markers$analysis_note <- paste(
  "Exploratory cell-level comparison; one pooled 10x library per genotype,",
  "so cell-level P values are not genotype-level biological-replicate inference."
)
write.csv(markers, file.path(output_dir, "Treg_DE.csv"))

# --- 4. GO enrichment (biological process) ---------------------------------
library(clusterProfiler)
go_genes <- markers$gene[markers$avg_log2FC > 0.25]
entrez <- bitr(go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
go.results <- enrichGO(entrez$ENTREZID, keyType = "ENTREZID", ont = "BP", OrgDb = org.Mm.eg.db)
write.csv(go.results@result, file.path(output_dir, "Treg_GO_BP.csv"))

# --- 5. descriptive GSVA pathway activity of Tregs (Kras vs WT) ------------
# One pooled library was generated per genotype. Therefore pathway scores are
# computed from genotype-level average log-normalized expression and reported
# descriptively; no cell-level limma P values are produced.
library(msigdbr)
Treg_avg <- AverageExpression(Treg, assays = "RNA", group.by = "orig.ident",
                              slot = "data", verbose = FALSE)[[1]]
Treg_avg <- as.matrix(Treg_avg[rowSums(Treg_avg) > 0, , drop = FALSE])

mouse_GO_bp_Set <- get_msigdb_sets("Mus musculus", "C5", "GO:BP")

Treg_gsva <- run_gsva(Treg_avg, mouse_GO_bp_Set, kcdf = "Gaussian")
write.table(Treg_gsva, file.path(output_dir, "Treg_GSVA_matrix.tsv"),
            sep = "\t", quote = FALSE)

required_groups <- c("Kras", "WT")
if (!all(required_groups %in% colnames(Treg_gsva))) {
  stop("Treg GSVA output does not contain both Kras and WT columns.")
}
diff <- data.frame(
  pathway = rownames(Treg_gsva),
  Kras = Treg_gsva[, "Kras"],
  WT = Treg_gsva[, "WT"],
  Kras_minus_WT = Treg_gsva[, "Kras"] - Treg_gsva[, "WT"],
  inference = "descriptive_only_one_pooled_library_per_genotype",
  row.names = NULL
)
write.csv(diff, file.path(output_dir, "Treg_GSVA_descriptive.csv"), row.names = FALSE)
