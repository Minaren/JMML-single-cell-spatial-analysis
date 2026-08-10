# ============================================================================
# Script: 00_setup.R
# Purpose: Project-wide configuration shared by all analysis scripts.
#          Defines the project root and standard sub-directories, loads the
#          R packages used across the pipeline, and fixes a global seed.
# Usage:   source("scripts/00_setup.R") at the start of each analysis script.
#          The project root is detected from the .here file located in the
#          package root; all input/output paths are relative to it.
# Prerequisites: R >= 4.2 and the packages listed below (see README.md and
#          environment/packages.R for the full list and version guidance).
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)          # single-cell and spatial analysis (v4.1.0 used in the study)
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
# spacexr (RCTD) and CellChat are loaded inside scripts 07 and 08
})

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
