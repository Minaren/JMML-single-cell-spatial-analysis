# ============================================================================
# environment/packages.R
# Lists all R packages required by the analysis scripts and reports the
# installed versions (adapted from the original analysis code).
# Recommended: install with renv (renv::init(); renv::install(pkgs);
# renv::snapshot()) to lock a reproducible environment. R >= 4.5 and
# Bioconductor 3.21 are required for GSVA v2.2.0.
# ============================================================================

pkgs <- c(
  # core single-cell / spatial
  "Seurat", "harmony", "spacexr", "slingshot", "mgcv", "RANN", "igraph",
  # data wrangling / plotting
  "tidyverse", "dplyr", "tidyr", "ggplot2", "patchwork", "ggpubr",
  "cowplot", "stringr", "pheatmap", "viridis", "ggrepel",
  # enrichment / pathway
  "GSVA", "msigdbr", "clusterProfiler", "DESeq2", "org.Mm.eg.db", "AnnotationDbi",
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
#   - bcl2fastq2 v2.20.0 (Illumina)          : BCL-to-FASTQ conversion
#   - Cell Ranger v7.0.0 (10x Genomics)      : scRNA-seq preprocessing
#   - BSTMatrix v1.0 (Biomarker Technologies): spatial transcriptomics upstream
#     processing and read mapping (mouse reference genome mm10; version to confirm)
#   - CellPhoneDB v5: exact installed version and database checksum are captured
#     by scripts/08_run_cellphonedb.py
#   - 10x Genomics Chromium Single Cell 3' v3 : library chemistry
#   - Illumina NovaSeq 6000 (PE150)          : sequencing platform
# ----------------------------------------------------------------------------
