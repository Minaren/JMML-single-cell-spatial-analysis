# ============================================================================
# Script: 08_cell_cell_communication.R
# Purpose: Visualisation of CellPhoneDB communication results inferred from
#          the spatial transcriptomic data: interaction networks and
#          HSC-centred ligand-receptor dot plots.
# Inputs:  output/cellphonedb/ST_WT/ and output/cellphonedb/ST_Kras/
#          (CellPhoneDB output files: count_network.txt, pvalues.txt,
#          means.txt, significant_means.txt)
#          output/spatial/ST_*_RCTD-annotated spot labels (Spatial_CellType.tsv)
# Outputs: figures/CellPhoneDB_net_circle_*.pdf,
#          figures/CellPhoneDB_DotPlot_HSC_*.pdf
# Run order: after 07_spatial_transcriptomics_analysis.R, after running
#          CellPhoneDB (see README.md "Cell-cell communication" section for
#          the required command)
# NOTE 1: CellPhoneDB is run separately after mouse-to-human orthologue
#         conversion by 08_prepare_cellphonedb_inputs.py and the versioned
#         08_run_cellphonedb.py entry point, which records package/database
#         versions and checksums.
# NOTE 2: HSC-centred dot plots are filtered at p.cutoff = 0.05.
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(igraph))

cpdb_dir <- file.path(output_dir, "cellphonedb")

# --- 1. CellPhoneDB network circle plots and HSC dot plots -------------------
cellphoneDB_Dotplot <- function(pvals.data, means.data, key,
                                target.cells_1, p.cutoff = 0.05) {
  colnames(pvals.data) <- str_replace_all(colnames(pvals.data), "\\.", "_")
  colnames(means.data) <- str_replace_all(colnames(means.data), "\\.", "_")
  pair_cols <- intersect(grep("\\|", colnames(pvals.data), value = TRUE),
                         grep("\\|", colnames(means.data), value = TRUE))
  keep_pair <- Reduce(`|`, lapply(target.cells_1, function(cell) {
    grepl(paste0("(^|\\|)", cell, "(\\||$)"), pair_cols)
  }))
  pair_cols <- pair_cols[keep_pair]
  if (!length(pair_cols)) stop("No HSC-centred CellPhoneDB cell-pair columns were found.")

  id_cols <- intersect(c("id_cp_interaction", "interacting_pair"),
                       intersect(colnames(pvals.data), colnames(means.data)))
  if (!"interacting_pair" %in% id_cols) {
    stop("CellPhoneDB outputs lack the interacting_pair column.")
  }
  p_long <- pvals.data %>%
    select(all_of(c(id_cols, pair_cols))) %>%
    pivot_longer(all_of(pair_cols), names_to = "cell_pair", values_to = "p_value")
  m_long <- means.data %>%
    select(all_of(c(id_cols, pair_cols))) %>%
    pivot_longer(all_of(pair_cols), names_to = "cell_pair", values_to = "mean")
  df <- inner_join(p_long, m_long, by = c(id_cols, "cell_pair")) %>%
    filter(is.finite(p_value), is.finite(mean), p_value < p.cutoff)
  if (!nrow(df)) stop("No HSC-centred interactions passed P < ", p.cutoff, ".")

  ggplot(df, aes(x = cell_pair, y = interacting_pair)) +
    geom_point(aes(size = -log10(pmax(p_value, 1e-300)),
                   colour = log2(mean + 1))) +
    scale_colour_gradientn(colours = c("#3A5978", "#F6B31D", "#DA2328")) +
    theme_bw() + labs(x = "", y = "", title = paste("CellPhoneDB:", key))
}

for (s in c("ST_WT", "ST_Kras")) {
  d <- file.path(cpdb_dir, s)
  require_files(file.path(d, c("count_network.txt", "pvalues.txt", "means.txt")))
  df.net <- read.table(file.path(d, "count_network.txt"), header = TRUE, sep = "\t")
  required_network_cols <- c("SOURCE", "TARGET", "count")
  if (!all(required_network_cols %in% colnames(df.net))) {
    stop("count_network.txt for ", s, " lacks SOURCE, TARGET or count.")
  }
  cell_types <- sort(unique(c(df.net$SOURCE, df.net$TARGET)))
  df.net <- df.net %>%
    complete(SOURCE = cell_types, TARGET = cell_types, fill = list(count = 0)) %>%
    arrange(match(SOURCE, cell_types), match(TARGET, cell_types)) %>%
    pivot_wider(names_from = TARGET, values_from = count, values_fill = 0)
  rownames(df.net) <- df.net$SOURCE
  df.net <- as.matrix(df.net[, cell_types, drop = FALSE])
  pvals_stat <- read.delim(file.path(d, "pvalues.txt"), check.names = FALSE)
  means_stat <- read.delim(file.path(d, "means.txt"), check.names = FALSE)

  pdf(file.path(fig_dir, paste0("CellPhoneDB_net_circle_", s, ".pdf")))
  graph <- igraph::graph_from_adjacency_matrix(df.net, mode = "directed",
                                               weighted = TRUE, diag = FALSE)
  edge_weights <- igraph::E(graph)$weight
  edge_width <- if (length(edge_weights) && max(edge_weights) > 0) {
    0.5 + 4.5 * edge_weights / max(edge_weights)
  } else {
    0.5
  }
  plot(graph, layout = igraph::layout_in_circle(graph), edge.width = edge_width,
       edge.arrow.size = 0.25, vertex.size = 18, vertex.label.cex = 0.75,
       main = paste("CellPhoneDB interaction counts:", s))
  dev.off()

  pdf(file.path(fig_dir, paste0("CellPhoneDB_DotPlot_HSC_", s, ".pdf")),
      width = 8, height = 10)
  print(cellphoneDB_Dotplot(pvals_stat, means_stat, key = "HSC", target.cells_1 = c("HSC")))
  dev.off()
}
