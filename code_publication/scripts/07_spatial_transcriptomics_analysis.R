# ============================================================================
# Script: 07_spatial_transcriptomics_analysis.R
# Purpose: Spatial transcriptomic data processing: spot-level cell-type
#          deconvolution with RCTD (spacexr) using the merged mouse bone-marrow
#          single-cell atlas as reference, and HSC-neighborhood composition
#          analysis.
# Inputs:  data/spatial/ST_WT/ and data/spatial/ST_Kras/
#          (BMKMANU S1000 output: matrix.mtx.gz, barcodes.tsv.gz,
#          features.tsv.gz, barcodes_pos.tsv.gz, barcodes_read.tsv.gz;
#          see data/README.md for the format note)
#          data/spatial/sc_meta.txt      (single-cell reference counts,
#                                         genes x cells)
#          data/spatial/ref_cell_anno    (single-cell reference annotation,
#                                         barcode <TAB> cell_type)
# Outputs: output/spatial/ST_WT/ and output/spatial/ST_Kras/
#          (RCTD.rds, Spatial_CellType.tsv, HSC_neighbor_cell_proportion.tsv)
# Run order: after 00_setup.R; single-cell reference prepared by scripts 01-03
# NOTE 1: Spatial experiment was performed by Biomarker Technologies (BMKMANU
#         S1000); raw data are not public (contact the authors / company).
# NOTE 2: The reference cell-type annotation used for RCTD must contain the
#         cell types of interest (e.g., HSC); cell types with <=25 cells in
#         the reference are excluded (as in the original analysis).
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(spacexr))   # RCTD deconvolution

spatial_out <- file.path(output_dir, "spatial")
dir.create(spatial_out, showWarnings = FALSE, recursive = TRUE)

scmeta      <- file.path(spatial_dir, "sc_meta.txt")
anno        <- file.path(spatial_dir, "ref_cell_anno")
min_cells   <- 5
min_features<- 100
HSC_name    <- "HSC"
radius      <- 100    # HSC neighborhood radius (coordinate units)

# --- 1. RCTD reference -------------------------------------------------------
sc_counts <- read.table(scmeta, header = TRUE, row.names = 1, check.names = FALSE)
sc_nUMI   <- colSums(sc_counts)
cellType  <- read.table(anno, header = FALSE, sep = "\t", check.names = FALSE)
colnames(cellType) <- c("barcode", "cell_type")
cellType$cell_type <- as.factor(cellType$cell_type)
cellType <- cellType[cellType$cell_type %in%
                       names(table(cellType$cell_type)[table(cellType$cell_type) > 25]), ]
cell_types <- setNames(as.character(cellType$cell_type), cellType$barcode)
cell_types <- as.factor(cell_types)
sc_counts <- sc_counts[, colnames(sc_counts) %in% names(cell_types)]
sc_nUMI   <- sc_nUMI[names(cell_types)]
reference <- Reference(sc_counts, cell_types, sc_nUMI)

# --- 2. per-sample RCTD deconvolution ----------------------------------------
for (s in c("ST_WT", "ST_Kras")) {
  cat("Processing", s, "\n")
  stmat <- file.path(spatial_dir, s)
  coords <- read.table(gzfile(file.path(stmat, "barcodes_pos.tsv.gz")),
                       sep = "\t", header = FALSE)
  colnames(coords) <- c("barcode", "x", "y")
  coords$y <- -coords$y              # y-flip applied in the original analysis
  rownames(coords) <- coords$barcode

  expr <- Read10X(stmat, gene.column = 2)
  sp_obj <- CreateSeuratObject(expr, assay = "Spatial",
                               min.cells = min_cells, min.features = min_features)
  sp_counts <- as.matrix(sp_obj@assays$Spatial@counts)
  sp_nUMI   <- colSums(sp_counts)
  puck <- SpatialRNA(coords, sp_counts, sp_nUMI)

  myRCTD <- create.RCTD(puck, reference, max_cores = 8, CELL_MIN_INSTANCE = 20)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  saveRDS(myRCTD, file.path(spatial_out, s, "RCTD.rds"))

  res_df <- myRCTD@results$results_df
  valid_bc <- rownames(res_df[res_df$spot_class != "reject" & puck@nUMI >= 1, ])
  anno_df <- puck@coords[valid_bc, ]
  anno_df$cell_type <- res_df[valid_bc, "first_type"]
  write.table(anno_df %>% rownames_to_column("barcode"),
              file.path(spatial_out, s, "Spatial_CellType.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --- 3. HSC neighborhood composition ---------------------------------------
  hsc_coords <- anno_df[anno_df$cell_type == HSC_name, ]
  dist_mat <- as.matrix(dist(anno_df[, c("x", "y")]))
  neighbor_idx <- which(apply(dist_mat[, rownames(hsc_coords)], 1,
                              function(x) any(x <= radius)))
  neighbor_cells <- anno_df[neighbor_idx, ]
  neighbor_cells <- neighbor_cells[neighbor_cells$cell_type != HSC_name, ]
  prop_df <- neighbor_cells %>% count(cell_type) %>% mutate(prop = n / sum(n))
  write.table(prop_df, file.path(spatial_out, s, "HSC_neighbor_cell_proportion.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
