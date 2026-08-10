# ============================================================================
# Script: 06_human_HSC_survival_analysis.R
# Purpose: Human HSC sub-clustering (resolution 0.8), CD69high/CD69low subset
#          classification, proliferation scoring, GSVA pathway activity,
#          CD69high HSC gene-signature construction, bulk-tumor ssGSEA scoring
#          (GSE71449), and survival analysis (KM, log-rank, univariate and
#          multivariate Cox).
# Inputs:  output/human_bone_marrow.rds
#          data/bulk/GSE71449/ids_exprs.csv        (bulk expression matrix)
#          data/bulk/GSE71449/Table_S1.xlsx        (clinical data)
# Outputs: output/human_HSC.rds
#          output/HSC_CD69high_vs_low_markers.csv
#          output/HSC_CD69high_survival_table.csv
#          figures/KM_OS_HSC_CD69high.pdf, figures/Forest_plot_multivariate.pdf
# Run order: after 05_human_data_integration.R
# NOTE: cluster 0 = CD69high and cluster 1 = CD69low in the human HSC
#       sub-clustering, as in the original analysis.
# ============================================================================

source("scripts/00_setup.R")

sce.harm <- readRDS(file.path(output_dir, "human_bone_marrow.rds"))

# --- 1. HSC sub-clustering ---------------------------------------------------
HSC <- subset(sce.harm, idents = "HSC")
HSC <- ScaleData(HSC)
HSC <- RunPCA(HSC)
HSC <- FindNeighbors(HSC, dims = 1:30)
HSC <- FindClusters(HSC, resolution = 0.8)
HSC <- RunUMAP(HSC, dims = 1:30)
Idents(HSC) <- "seurat_clusters"

# CD69high/CD69low assignment (cluster 0 = high, cluster 1 = low)
HSC$CD69_status <- ifelse(HSC$seurat_clusters == "0", "CD69high",
                   ifelse(HSC$seurat_clusters == "1", "CD69low", "Other"))

# --- 2. proliferation score ---------------------------------------------------
proliferation_gene_sets <- list(
  IEGs = c("FOS", "JUN", "JUNB", "JUND", "EGR1", "REL", "NFKBIA", "NFKBIE"),
  Cell_Cycle_Regulators = c("CCNG2", "CDK1", "CDK2", "CDK4", "CDK6", "CCNB1", "CCNB2",
                            "CCNA2", "CDC20", "CDC25A", "CDC25B", "CDC25C", "AURKB", "KIF23"),
  DNA_Replication_Mitosis = c("PCNA", "TOP2A", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
                              "UBE2C", "BIRC5", "TYMS", "RRM2", "CENPA", "CENPE", "CENPF"),
  Transcription_Factors = c("MYB", "FOXM1", "SOX4", "GATA2", "STAT5A", "ID1", "E2F1"),
  Signaling_Checkpoint_Others = c("TNF", "TNFAIP3", "JAK1", "KDM5B", "KDM6B",
                                  "CHEK1", "CHEK2", "CDK19")
)
HSC <- AddModuleScore(HSC, features = list(unique(unlist(proliferation_gene_sets))),
                      name = "Proliferation_Score")

# --- 3. GSVA pathway activity (CD69high vs CD69low) --------------------------
# NOTE: the original human HSC GSVA block referenced workspace variables
#       (HSC_counts_filtered, mouse GO:BP sets, CD69_status). These are
#       rebuilt here deterministically; MSigDB human C5 GO:BP sets are used
#       (the original code referenced the mouse set by mistake).
library(msigdbr)
HSC_counts <- as.matrix(GetAssayData(HSC, assay = "RNA", layer = "counts"))
cells_use <- WhichCells(HSC, idents = c("0", "1"))
human_GO_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)
human_GO_bp_Set <- split(human_GO_bp$gene_symbol, human_GO_bp$gs_name)
HSC_gsva <- gsva(HSC_counts[, cells_use], human_GO_bp_Set, kcdf = "Poisson", parallel.sz = 4)

group <- factor(HSC@meta.data[cells_use, "CD69_status"])
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
fit <- lmFit(HSC_gsva, design)
cont <- makeContrasts(CD69high_vs_low = CD69high - CD69low, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
diff <- topTable(fit2, adjust = "fdr", number = Inf)
diff$ID <- rownames(diff)
write.csv(diff, file.path(output_dir, "human_HSC_GSVA_diff.csv"))

# --- 4. CD69high HSC gene signature ------------------------------------------
Idents(HSC) <- "seurat_clusters"
HSC <- JoinLayers(HSC)
hsc_markers <- FindAllMarkers(HSC, only.pos = TRUE, logfc.threshold = 0.25,
                              min.pct = 0.1, test.use = "wilcox")
write.csv(hsc_markers, file.path(output_dir, "HSC_CD69high_vs_low_markers.csv"), row.names = FALSE)
genes_high <- hsc_markers$gene[hsc_markers$cluster == "0"]
genes_low  <- hsc_markers$gene[hsc_markers$cluster == "1"]

saveRDS(HSC, file.path(output_dir, "human_HSC.rds"))

# --- 5. bulk ssGSEA scoring and survival -------------------------------------
bulk_raw <- read.delim(file.path(bulk_dir, "GSE71449", "ids_exprs.csv"),
                       sep = ",", row.names = 1, check.names = FALSE)
bulk_expr_mat <- as.matrix(bulk_raw); mode(bulk_expr_mat) <- "numeric"
rownames(bulk_expr_mat) <- toupper(rownames(bulk_expr_mat))
if (any(duplicated(rownames(bulk_expr_mat)))) {
  bulk_df <- as.data.frame(bulk_expr_mat)
  bulk_df$Gene <- rownames(bulk_df)
  bulk_df <- aggregate(. ~ Gene, data = bulk_df, FUN = mean)
  rownames(bulk_df) <- bulk_df$Gene; bulk_df$Gene <- NULL
  bulk_expr_mat <- as.matrix(bulk_df)
}
gene_sets <- list(HSC_CD69high = intersect(toupper(genes_high), rownames(bulk_expr_mat)),
                  HSC_CD69low  = intersect(toupper(genes_low),  rownames(bulk_expr_mat)))
if (packageVersion("GSVA") >= "1.50.0") {
  gsva_res <- gsva(gsvaParam(bulk_expr_mat, gene_sets, kcdf = "Gaussian"), verbose = FALSE)
} else {
  gsva_res <- gsva(bulk_expr_mat, gene_sets, method = "ssgsea", kcdf = "Gaussian")
}
hsc_high_score <- gsva_res["HSC_CD69high", ]

clinical <- readxl::read_excel(file.path(bulk_dir, "GSE71449", "Table_S1.xlsx"))
clinical$SampleID <- as.character(clinical$ID)
clinical$diag_survival_days <- as.numeric(clinical[["Survival from diagnosis\r\n(days)"]])
clinical$diag_event <- ifelse(!is.na(clinical[["Cause of death"]]) &
                                clinical[["Cause of death"]] != "", 1, 0)
clinical$HSC_CD69high_score <- hsc_high_score[match(clinical$SampleID, names(hsc_high_score))]
clinical_valid <- clinical[!is.na(clinical$HSC_CD69high_score) &
                             !is.na(clinical$diag_survival_days) &
                             !is.na(clinical$diag_event), ]

# median split and optimal cutoff
clinical_valid$group_high <- ifelse(clinical_valid$HSC_CD69high_score >
                                      median(clinical_valid$HSC_CD69high_score), "High", "Low")
cut <- surv_cutpoint(clinical_valid, time = "diag_survival_days", event = "diag_event",
                     variables = "HSC_CD69high_score")
clinical_valid$group_optimal <- ifelse(clinical_valid$HSC_CD69high_score > cut$cutpoint$cutpoint,
                                       "High", "Low")

# KM + log-rank (median split)
surv_obj <- Surv(time = clinical_valid$diag_survival_days, event = clinical_valid$diag_event)
fit <- survfit(surv_obj ~ group_high, data = clinical_valid)
km <- ggsurvplot(fit, data = clinical_valid, pval = TRUE, risk.table = TRUE,
                 palette = "jco", xlab = "Days from diagnosis", ylab = "Overall Survival",
                 legend.title = "HSC_CD69high score")
pdf(file.path(fig_dir, "KM_OS_HSC_CD69high.pdf"), width = 8, height = 7)
print(km)
dev.off()

# Cox models
clinical_valid <- clinical_valid %>%
  mutate(age_years = as.numeric(.data[["age at diagnosis\r\n(years)"]]),
         monosomy7 = ifelse(Karyotype == "mono7", 1, 0),
         mut_group = factor(case_when(Mutation == "PTPN11" ~ "PTPN11",
                                      Mutation == "quad neg" ~ "Quad_neg",
                                      TRUE ~ "Other"),
                            levels = c("PTPN11", "Other", "Quad_neg")))
cox_multi <- coxph(Surv(diag_survival_days, diag_event) ~ HSC_CD69high_score + age_years +
                     monosomy7 + mut_group, data = clinical_valid)
write.csv(as.data.frame(summary(cox_multi)$conf.int),
          file.path(output_dir, "HSC_CD69high_multivariate_Cox.csv"))
pdf(file.path(fig_dir, "Forest_plot_multivariate.pdf"), width = 9, height = 5)
ggforest(cox_multi, data = clinical_valid)
dev.off()