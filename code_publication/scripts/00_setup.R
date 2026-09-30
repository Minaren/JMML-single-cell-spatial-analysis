# ============================================================================
# Script: 00_setup.R
# Purpose: Project-wide configuration shared by all analysis scripts.
#          Defines the project root and standard sub-directories, loads the
#          R packages used across the pipeline, and fixes a global seed.
# Usage:   source("scripts/00_setup.R") at the start of each analysis script.
#          The project root is detected from the .here file located in the
#          package root; all input/output paths are relative to it.
# Prerequisites: R >= 4.5 and the packages listed below (see README.md and
#          environment/packages.R for the full list and version guidance).
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)          # single-cell and spatial analysis (v4.3.0 used in the study)
  library(harmony)         # batch integration
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(cowplot)
  library(stringr)
  library(pheatmap)
  library(GSVA)            # gene set variation analysis / ssGSEA (v2.2.0 used)
  library(slingshot)       # trajectory inference (HSC subsets)
  library(mgcv)            # GAM smoothing of pseudotime trends
  library(org.Mm.eg.db)    # mouse gene annotations
  library(AnnotationDbi)
  library(survival)        # KM and Cox regression (human survival analysis)
  library(survminer)
# spacexr (RCTD) is loaded inside script 07
})

# The archived analysis is documented against Seurat 4.3.0 and GSVA 2.2.0.
# Fail early instead of silently switching between incompatible Seurat v4/v5
# assay-layer semantics or between the retired and current GSVA APIs.
if (packageVersion("Seurat") < "4.3.0" || packageVersion("Seurat") >= "5.0.0") {
  stop("This release requires Seurat >= 4.3.0 and < 5.0.0.")
}
if (packageVersion("GSVA") < "2.2.0") {
  stop("This release requires GSVA >= 2.2.0 and the parameter-object API.")
}

# --- project paths ---------------------------------------------------------
project_root <- here::here()          # requires the .here marker file in the package root
data_dir     <- file.path(project_root, "data")
raw_dir      <- file.path(data_dir, "raw")
gene_dir     <- file.path(data_dir, "gene_sets")
ref_dir      <- file.path(data_dir, "reference")
spatial_dir  <- file.path(data_dir, "spatial")
bulk_dir     <- file.path(data_dir, "bulk")
output_dir   <- file.path(project_root, "output")
fig_dir      <- file.path(project_root, "figures")

for (d in c(output_dir, fig_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(220625)          # seed used in the original analysis

# --- shared validation and compatibility helpers --------------------------
require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Required input file(s) not found:\n", paste(missing, collapse = "\n"))
  }
  invisible(paths)
}

get_assay_matrix <- function(object, assay = "RNA", slot = "data") {
  Seurat::GetAssayData(object = object, assay = assay, slot = slot)
}

read_gene_sets <- function(path, min_size = 5L) {
  require_files(path)
  x <- read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                col.names = c("set", "gene"))
  x$set <- trimws(x$set)
  x$gene <- trimws(x$gene)
  x <- x[nzchar(x$set) & nzchar(x$gene), , drop = FALSE]

  # Strings such as "1-Sep" and "7-Mar" are characteristic of spreadsheet
  # date conversion of gene symbols.  They must be restored from the source
  # gene-set file rather than guessed or silently discarded.
  date_like <- grepl("^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$",
                     x$gene, ignore.case = TRUE)
  if (any(date_like)) {
    bad <- unique(x$gene[date_like])
    stop("Possible spreadsheet-corrupted gene symbols in ", basename(path),
         ": ", paste(bad, collapse = ", "),
         ". Restore these symbols from the original source before analysis.")
  }

  sets <- lapply(split(x$gene, x$set), unique)
  sets <- sets[lengths(sets) >= min_size]
  if (!length(sets)) stop("No gene set with at least ", min_size, " genes in ", path)
  sets
}

run_gsva <- function(expr, gene_sets, kcdf = "Gaussian") {
  param <- GSVA::gsvaParam(exprData = expr, geneSets = gene_sets, kcdf = kcdf)
  GSVA::gsva(param, verbose = FALSE)
}

run_ssgsea <- function(expr, gene_sets) {
  param <- GSVA::ssgseaParam(exprData = expr, geneSets = gene_sets)
  GSVA::gsva(param, verbose = FALSE)
}

get_msigdb_sets <- function(species, collection = "C5", subcollection = "GO:BP") {
  api <- names(formals(msigdbr::msigdbr))
  if (all(c("collection", "subcollection") %in% api)) {
    tbl <- msigdbr::msigdbr(species = species, collection = collection,
                            subcollection = subcollection)
  } else {
    tbl <- msigdbr::msigdbr(species = species, category = collection,
                            subcategory = subcollection)
  }
  if (!all(c("gs_name", "gene_symbol") %in% colnames(tbl))) {
    stop("Unexpected msigdbr output: gs_name or gene_symbol is absent.")
  }
  split(tbl$gene_symbol, tbl$gs_name)
}

check_group_counts <- function(object, group_col, expected) {
  observed <- table(object@meta.data[[group_col]])
  absent <- setdiff(names(expected), names(observed))
  mismatch <- length(absent) || any(observed[names(expected)] != expected)
  if (mismatch) {
    stop("Unexpected retained-cell counts for ", group_col, ". Observed: ",
         paste(names(observed), as.integer(observed), sep = "=", collapse = ", "),
         "; expected: ",
         paste(names(expected), expected, sep = "=", collapse = ", "), ".")
  }
  invisible(observed)
}
