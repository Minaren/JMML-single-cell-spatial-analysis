# ============================================================================
# Script: 08_cell_cell_communication.R
# Purpose: Visualisation and interpretation of cell-cell communication
#          inferred from the spatial transcriptomic data:
#          (1) CellPhoneDB interaction networks (circle plots) and
#              HSC-centred ligand-receptor dot plots;
#          (2) CellChat-based communication inference on the deconvoluted
#              spatial data.
# Inputs:  output/cellphonedb/ST_WT/ and output/cellphonedb/ST_Kras/
#          (CellPhoneDB output files: counts_network.txt, pvalues.txt,
#          means.txt, significant_means.txt)
#          output/spatial/ST_*_RCTD-annotated spot labels (Spatial_CellType.tsv)
# Outputs: figures/CellPhoneDB_net_circle_*.pdf,
#          figures/CellPhoneDB_DotPlot_HSC_*.pdf
# Run order: after 07_spatial_transcriptomics_analysis.R, after running
#          CellPhoneDB (see README.md "Cell-cell communication" section for
#          the required command)
# NOTE 1: CellPhoneDB (v4+) is a Python tool and must be run separately on the
#         RCTD-annotated spots (counts + meta files). Its output files are
#         then read here for visualisation. Version and database used by the
#         original analysis are to be confirmed by the authors.
# NOTE 2: HSC-centred dot plots are filtered at p.cutoff = 0.05.
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(CellChat))  # v1.6.1 used in the study
suppressPackageStartupMessages(library(reshape2))

cpdb_dir <- file.path(output_dir, "cellphonedb")

# --- 1. CellPhoneDB network circle plots and HSC dot plots -------------------
cellphoneDB_Dotplot <- function(pvals.data, means.data, key,
                                target.cells_1, p.cutoff = 0.05) {
  colnames(pvals.data) <- str_replace_all(colnames(pvals.data), "\\.", "_")
  colnames(means.data) <- str_replace_all(colnames(means.data), "\\.", "_")
  kp <- Reduce(`|`, lapply(target.cells_1, grepl, x = colnames(pvals.data)))
  pos <- which(kp)
  pvals <- pvals.data[, c(1, 2, 5, 6, 8, 9, pos)]
  means <- means.data[, c(1, 2, 5, 6, 8, 9, pos)]
  pvals <- pvals[rowSums(pvals[, 7:ncol(pvals)] < p.cutoff) > 0, ]
  means <- means[means$id_cp_interaction %in% pvals$id_cp_interaction, ]
  df <- merge(reshape2::melt(pvals, id.vars = "interacting_pair"),
              reshape2::melt(means, id.vars = "interacting_pair"),
              by = c("interacting_pair", "variable"))
  ggplot(df, aes(variable, interacting_pair)) +
    geom_point(aes(size = -log10(value.x + 1e-4), colour = log2(value.y + 1))) +
    scale_colour_gradientn(colours = c("#3A5978", "#F6B31D", "#DA2328")) +
    theme_bw() + labs(x = "", y = "", title = paste("CellPhoneDB:", key))
}

for (s in c("ST_WT", "ST_Kras")) {
  d <- file.path(cpdb_dir, s)
  df.net <- read.table(file.path(d, "count_network.txt"), header = TRUE, sep = "\t")
  df.net <- spread(df.net, TARGET, count)
  rownames(df.net) <- df.net$SOURCE
  df.net <- as.matrix(df.net[, -1])
  pvals_stat <- read.delim(file.path(d, "pvalues.txt"), check.names = FALSE)
  means_stat <- read.delim(file.path(d, "means.txt"), check.names = FALSE)

  pdf(file.path(fig_dir, paste0("CellPhoneDB_net_circle_", s, ".pdf")))
  netVisual_circle(df.net, weight.scale = TRUE, label.edge = FALSE)
  dev.off()

  pdf(file.path(fig_dir, paste0("CellPhoneDB_DotPlot_HSC_", s, ".pdf")),
      width = 8, height = 10)
  print(cellphoneDB_Dotplot(pvals_stat, means_stat, key = "HSC", target.cells_1 = c("HSC")))
  dev.off()
}

# --- 2. CellChat inference on deconvoluted spatial data ----------------------
# CellChat was applied to the RCTD-annotated spatial data (CellChatDB.mouse).
# Build a meta object per sample from the spatial annotations and run the
# standard CellChat workflow (createCellChat -> identifyOverExpressedGenes ->
# computeCommunProb -> filterCommunication -> netVisual_circle).
# NOTE: exact CellChat parameters (nboot, type) used in the original analysis
#       are to be confirmed; defaults of CellChat v1.6.1 were used.
