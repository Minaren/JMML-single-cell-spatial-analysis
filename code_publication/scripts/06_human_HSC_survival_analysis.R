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
HSC <- FindVariableFeatures(HSC, selection.method = "vst", nfeatures = 2000)
HSC <- ScaleData(HSC, features = VariableFeatures(HSC))
HSC <- RunPCA(HSC, features = VariableFeatures(HSC))
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

# --- 3. descriptive GSVA pathway activity (CD69high vs CD69low) ------------
# Scores are calculated from average log-normalized expression for the two
# HSC states. The underlying public analysis contains no independent library
# replicates for cell-level pathway inference, so no cell-as-replicate limma
# P values are reported.
library(msigdbr)
cells_use <- rownames(HSC@meta.data)[HSC$CD69_status %in% c("CD69high", "CD69low")]
HSC_pair <- subset(HSC, cells = cells_use)
HSC_avg <- AverageExpression(HSC_pair, assays = "RNA", group.by = "CD69_status",
                             slot = "data", verbose = FALSE)[[1]]
HSC_avg <- as.matrix(HSC_avg[rowSums(HSC_avg) > 0, , drop = FALSE])
human_GO_bp_Set <- get_msigdb_sets("Homo sapiens", "C5", "GO:BP")
HSC_gsva <- run_gsva(HSC_avg, human_GO_bp_Set, kcdf = "Gaussian")
if (!all(c("CD69high", "CD69low") %in% colnames(HSC_gsva))) {
  stop("Human HSC GSVA output lacks CD69high or CD69low.")
}
diff <- data.frame(
  pathway = rownames(HSC_gsva),
  CD69high = HSC_gsva[, "CD69high"],
  CD69low = HSC_gsva[, "CD69low"],
  CD69high_minus_CD69low = HSC_gsva[, "CD69high"] - HSC_gsva[, "CD69low"],
  inference = "descriptive_only_no_independent_library_replicates",
  row.names = NULL
)
write.csv(diff, file.path(output_dir, "human_HSC_GSVA_descriptive.csv"), row.names = FALSE)

# --- 4. CD69high HSC gene signature ------------------------------------------
Idents(HSC) <- "seurat_clusters"
hsc_markers <- FindMarkers(HSC, ident.1 = "0", ident.2 = "1",
                           only.pos = FALSE, logfc.threshold = 0.25,
                           min.pct = 0.1, test.use = "wilcox")
hsc_markers$gene <- rownames(hsc_markers)
write.csv(hsc_markers, file.path(output_dir, "HSC_CD69high_vs_low_markers.csv"), row.names = FALSE)
genes_high <- hsc_markers$gene[hsc_markers$avg_log2FC > 0.25 & hsc_markers$p_val_adj < 0.05]
genes_low  <- hsc_markers$gene[hsc_markers$avg_log2FC < -0.25 & hsc_markers$p_val_adj < 0.05]
if (length(genes_high) < 10L || length(genes_low) < 10L) {
  stop("Fewer than 10 significant genes in the CD69high or CD69low signature.")
}

saveRDS(HSC, file.path(output_dir, "human_HSC.rds"))

# --- 5. bulk ssGSEA scoring and survival -------------------------------------
bulk_expr_file <- file.path(bulk_dir, "GSE71449", "ids_exprs.csv")
clinical_file <- file.path(bulk_dir, "GSE71449", "Table_S1.xlsx")
require_files(c(bulk_expr_file, clinical_file))
bulk_raw <- read.delim(bulk_expr_file,
                       sep = ",", row.names = 1, check.names = FALSE)
bulk_expr_mat <- as.matrix(bulk_raw); mode(bulk_expr_mat) <- "numeric"
if (anyNA(bulk_expr_mat) || any(!is.finite(bulk_expr_mat))) {
  stop("GSE71449 expression matrix contains missing or non-finite values.")
}
if (anyDuplicated(colnames(bulk_expr_mat))) {
  stop("GSE71449 expression matrix contains duplicated sample IDs.")
}
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
if (any(lengths(gene_sets) < 10L)) {
  stop("Insufficient signature overlap with GSE71449 expression matrix: ",
       paste(names(gene_sets), lengths(gene_sets), sep = "=", collapse = ", "))
}
write.csv(data.frame(signature = names(gene_sets), n_genes = lengths(gene_sets)),
          file.path(output_dir, "HSC_signature_overlap_GSE71449.csv"), row.names = FALSE)
gsva_res <- run_ssgsea(bulk_expr_mat, gene_sets)
hsc_high_score <- gsva_res["HSC_CD69high", ]

clinical <- readxl::read_excel(clinical_file)
required_clinical <- c("ID", "Survival from diagnosis\r\n(days)",
                       "age at diagnosis\r\n(years)",
                       "Karyotype", "Mutation")
missing_clinical <- setdiff(required_clinical, names(clinical))
if (length(missing_clinical)) {
  stop("Missing required GSE71449 clinical column(s): ",
       paste(missing_clinical, collapse = ", "))
}
clinical$SampleID <- as.character(clinical$ID)
clinical$diag_survival_days <- as.numeric(clinical[["Survival from diagnosis\r\n(days)"]])
# Prefer an explicit vital/event-status field. If the public table has none,
# retain the original cause-of-death rule but record that fallback in the
# exported audit table so that it can be checked against the source dataset.
status_candidates <- c("diag_event", "OS_event", "Event", "event",
                       "Vital status", "Vital Status", "Status")
status_col <- status_candidates[status_candidates %in% names(clinical)][1]
if (!is.na(status_col)) {
  raw_status <- tolower(trimws(as.character(clinical[[status_col]])))
  clinical$diag_event <- ifelse(raw_status %in% c("1", "dead", "deceased", "death", "event"), 1,
                                ifelse(raw_status %in% c("0", "alive", "censored", "no event"), 0,
                                       NA_real_))
  event_definition <- paste0("explicit_status_column:", status_col)
} else {
  if (!"Cause of death" %in% names(clinical)) {
    stop("No explicit survival-status column or Cause of death column was found.")
  }
  clinical$diag_event <- ifelse(!is.na(clinical[["Cause of death"]]) &
                                  trimws(as.character(clinical[["Cause of death"]])) != "", 1, 0)
  event_definition <- "fallback_nonempty_cause_of_death"
  warning("No explicit survival-status column found; event status was inferred from non-empty Cause of death. Review before publication.")
}
clinical$HSC_CD69high_score <- hsc_high_score[match(clinical$SampleID, names(hsc_high_score))]
clinical_valid <- clinical[!is.na(clinical$HSC_CD69high_score) &
                             !is.na(clinical$diag_survival_days) &
                             !is.na(clinical$diag_event), ]

# Prespecified median split used in the manuscript.
clinical_valid$group_high <- ifelse(clinical_valid$HSC_CD69high_score >
                                      median(clinical_valid$HSC_CD69high_score), "High", "Low")
if (nrow(clinical_valid) != 44L) {
  warning("Expected 44 evaluable GSE71449 patients but found ", nrow(clinical_valid), ".")
}
clinical_valid$event_definition <- event_definition
write.csv(clinical_valid, file.path(output_dir, "HSC_CD69high_survival_table.csv"),
          row.names = FALSE)

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
         monosomy7 = ifelse(grepl("mono7|monosomy[ -]?7|(^|[^0-9])-7([^0-9]|$)",
                                  Karyotype, ignore.case = TRUE), 1, 0),
         mutation_normalized = tolower(trimws(as.character(Mutation))),
         mut_group = factor(case_when(grepl("^ptpn11", mutation_normalized) ~ "PTPN11",
                                      grepl("quad.*neg", mutation_normalized) ~ "Quad_neg",
                                      TRUE ~ "Other"),
                            levels = c("PTPN11", "Other", "Quad_neg")))

tidy_cox <- function(model, model_name) {
  sm <- summary(model)
  data.frame(
    model = model_name,
    term = rownames(sm$coefficients),
    HR = sm$conf.int[, "exp(coef)"],
    CI95_low = sm$conf.int[, "lower .95"],
    CI95_high = sm$conf.int[, "upper .95"],
    P_value = sm$coefficients[, "Pr(>|z|)"],
    n = model$n,
    events = model$nevent,
    row.names = NULL,
    check.names = FALSE
  )
}

# Supplementary Table S6: univariable Cox models.
univariable_formulas <- list(
  HSC_CD69high_score = Surv(diag_survival_days, diag_event) ~ HSC_CD69high_score,
  age_years = Surv(diag_survival_days, diag_event) ~ age_years,
  monosomy7 = Surv(diag_survival_days, diag_event) ~ monosomy7,
  mut_group = Surv(diag_survival_days, diag_event) ~ mut_group
)
cox_uni <- bind_rows(lapply(names(univariable_formulas), function(nm) {
  model <- coxph(univariable_formulas[[nm]], data = clinical_valid,
                 na.action = na.omit, x = TRUE)
  tidy_cox(model, paste0("univariable_", nm))
}))
write.csv(cox_uni, file.path(output_dir, "Table_S6_univariable_Cox.csv"), row.names = FALSE)

# Figure 8H: multivariable Cox model, fit on one explicit complete-case set.
model_vars <- c("diag_survival_days", "diag_event", "HSC_CD69high_score",
                "age_years", "monosomy7", "mut_group")
cox_data <- clinical_valid[complete.cases(clinical_valid[, model_vars]), ]
cox_multi <- coxph(Surv(diag_survival_days, diag_event) ~ HSC_CD69high_score +
                       age_years + monosomy7 + mut_group,
                     data = cox_data, x = TRUE)
cox_multi_tidy <- tidy_cox(cox_multi, "multivariable_Figure_8H")
write.csv(cox_multi_tidy,
          file.path(output_dir, "Figure_8H_multivariable_Cox.csv"), row.names = FALSE)

n_parameters <- sum(!is.na(coef(cox_multi)))
events_per_parameter <- cox_multi$nevent / n_parameters
model_audit <- data.frame(
  n_complete_cases = cox_multi$n,
  n_events = cox_multi$nevent,
  n_fitted_parameters = n_parameters,
  events_per_parameter = events_per_parameter
)
write.csv(model_audit, file.path(output_dir, "Figure_8H_model_audit.csv"),
          row.names = FALSE)
if (events_per_parameter < 10) {
  warning("Figure 8H has fewer than 10 events per fitted Cox parameter; interpret as exploratory.")
}

ph_test <- as.data.frame(cox.zph(cox_multi)$table)
ph_test$term <- rownames(ph_test)
write.csv(ph_test, file.path(output_dir, "Figure_8H_proportional_hazards_test.csv"),
          row.names = FALSE)
pdf(file.path(fig_dir, "Forest_plot_multivariate.pdf"), width = 9, height = 5)
ggforest(cox_multi, data = cox_data)
dev.off()
