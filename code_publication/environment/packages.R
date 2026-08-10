# ============================================================================
# environment/packages.R
# Lists all R packages required by the analysis scripts and reports the
# installed versions (adapted from the original analysis code).
# Recommended: install with renv (renv::init(); renv::install(pkgs);
# renv::snapshot()) to lock a reproducible environment. R >= 4.2 is required.
# ============================================================================

pkgs <- c(
  # core single-cell / spatial
  "Seurat", "harmony", "spacexr", "CellChat", "slingshot", "mgcv",
  # data wrangling / plotting
  "tidyverse", "dplyr", "tidyr", "ggplot2", "patchwork", "ggpubr",
  "cowplot", "stringr", "reshape2", "pheatmap", "viridis", "ggrepel",
  # enrichment / pathway
  "GSVA", "msigdbr", "clusterProfiler", "limma", "org.Mm.eg.db", "AnnotationDbi",
  # survival
  "survival", "survminer",
  # misc
  "gplots", "corrplot", "readxl", "data.table", "here"
)

# Report installed versions
pkg_info <- do.call(rbind, lapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    data.frame(Package = p, Version = as.character(packageVersion(p)), stringsAsFactors = FALSE)
  } else {
    data.frame(Package = p, Version = "NOT INSTALLED", stringsAsFactors = FALSE)
  }
}))
print(pkg_info, row.names = FALSE)
write.csv(pkg_info, file.path("output", "package_versions.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# External software used outside R (not covered by renv):
#   - Cell Ranger v6.1.1 (10x Genomics)      : scRNA-seq preprocessing
#   - BSTMatrix v1.0 (Biomarker Technologies): spatial transcriptomics upstream
#     processing and read mapping (mouse reference genome mm10; version to confirm)
#   - CellPhoneDB (v4+; version to confirm)  : ligand-receptor analysis
#   - 10x Genomics Chromium Single Cell 3' v3 : library chemistry
#   - Illumina NovaSeq 6000 (PE150)          : sequencing platform
# ----------------------------------------------------------------------------