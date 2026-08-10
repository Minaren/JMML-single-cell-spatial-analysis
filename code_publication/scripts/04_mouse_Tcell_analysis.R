# ============================================================================
# Script: 04_mouse_Tcell_analysis.R
# Purpose: Mouse bone-marrow CD4+ T-cell analysis: sub-clustering of CD4+
#          T cells, Treg/Tconv identification, differential expression between
#          genotypes, GO enrichment, and GSVA pathway activity of Tregs.
# Inputs:  data/raw/mouse_Tcell/            (10x Cell Ranger output, CD45+CD3+
#          sorted bone-marrow T cells; three mice pooled per genotype)
# Outputs: output/mouse_CD4Tcells.rds, output/Treg_DE.csv,
#          output/Treg_GO_BP.csv, output/Treg_GSVA_diff.csv
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
  sce_har <- RunHarmony(sc_n, group.by.vars = "orig.ident")
  pbmc_filt <- FindNeighbors(sce_har, dims = 1:30, reduction = "harmony")
  pbmc_filt <- FindClusters(pbmc_filt, resolution = 0.6, algorithm = 1)
  pbmc_filt <- RunUMAP(pbmc_filt, dims = 1:30, reduction = "harmony")
  pbmc_filt$orig.ident <- factor(pbmc_filt$orig.ident, levels = c("WT", "Kras"))
  saveRDS(pbmc_filt, tcell_obj_file)
}

# --- 2. CD4+ T-cell subsetting and Treg/Tconv assignment --------------------
# The CD4+ T-cell compartment corresponds to clusters 1 and 2 of the
# CD45+CD3+ object (from the original analysis).
CD4Tcells <- subset(pbmc_filt, seurat_clusters %in% c(1, 2))
CD4Tcells$celltype <- "Other"
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 1] <- "Treg"
CD4Tcells$celltype[CD4Tcells$seurat_clusters == 2] <- "Tconv"
Idents(CD4Tcells) <- "celltype"

saveRDS(CD4Tcells, file.path(output_dir, "mouse_CD4Tcells.rds"))

# --- 3. differential expression: Tregs, Kras vs WT --------------------------
Treg <- subset(CD4Tcells, celltype == "Treg")
Idents(Treg) <- "orig.ident"
markers <- FindMarkers(Treg, ident.1 = "Kras", ident.2 = "WT",
                       only.pos = FALSE, logfc.threshold = 0.25)
markers$gene <- rownames(markers)
markers <- markers[markers$p_val_adj < 0.05, ]
write.csv(markers, file.path(output_dir, "Treg_DE.csv"))

# --- 4. GO enrichment (biological process) ---------------------------------
library(clusterProfiler)
entrez <- bitr(markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
go.results <- enrichGO(entrez$ENTREZID, keyType = "ENTREZID", ont = "BP", OrgDb = org.Mm.eg.db)
write.csv(go.results@result, file.path(output_dir, "Treg_GO_BP.csv"))

# --- 5. GSVA pathway activity of Tregs (Kras vs WT) -------------------------
# MSigDB C5 GO:BP gene sets; kcdf = Poisson for count data.
library(msigdbr)
Treg_counts <- as.matrix(GetAssayData(Treg, assay = "RNA", layer = "counts"))
cells_to_use <- WhichCells(Treg, idents = c("Kras", "WT"))
Treg_counts <- Treg_counts[, cells_to_use]

mouse_GO_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
mouse_GO_bp_Set <- split(mouse_GO_bp$gene_symbol, mouse_GO_bp$gs_name)

Treg_gsva <- gsva(as.matrix(Treg_counts), mouse_GO_bp_Set, kcdf = "Poisson", parallel.sz = 4)
write.table(Treg_gsva, file.path(output_dir, "Treg_GSVA_matrix.xls"), sep = "\t")

group <- factor(Treg@meta.data[cells_to_use, "orig.ident"])
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
fit <- lmFit(Treg_gsva, design)
cont <- makeContrasts(Kras_vs_WT = Kras - WT, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
diff <- topTable(fit2, adjust = "fdr", number = Inf)
diff$ID <- rownames(diff)
write.csv(diff, file.path(output_dir, "Treg_GSVA_diff.csv"))
