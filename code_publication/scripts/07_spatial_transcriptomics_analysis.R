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
#          output/mouse_HSPC_annotated.rds (script 01 output; RCTD reference)
# Outputs: output/spatial/ST_WT/ and output/spatial/ST_Kras/
#          (RCTD.rds, Spatial_CellType.tsv, HSC_neighbor_cell_proportion.tsv)
# Run order: after script 01 (which prepares the annotated single-cell reference)
# NOTE 1: Spatial experiment was performed by Biomarker Technologies (BMKMANU
#         S1000); the deposited study data are available under GSE313878.
# NOTE 2: The reference cell-type annotation used for RCTD must contain the
#         cell types of interest (e.g., HSC); cell types with <=25 cells in
#         the reference are excluded (as in the original analysis).
# ============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(spacexr))   # RCTD deconvolution
suppressPackageStartupMessages(library(RANN))      # memory-safe nearest-neighbour search

spatial_out <- file.path(output_dir, "spatial")
dir.create(spatial_out, showWarnings = FALSE, recursive = TRUE)

reference_rds <- file.path(output_dir, "mouse_HSPC_annotated.rds")
min_cells   <- 5
min_features<- 100
HSC_name    <- "HSC"
radius      <- 100    # HSC neighborhood radius (coordinate units)

# --- 1. RCTD reference -------------------------------------------------------
# Rebuild the reference directly from the versioned script-01 output. This
# avoids the historical unversioned dense exports sc_meta.txt/ref_cell_anno and
# guarantees that the RCTD labels match the released annotation workflow.
require_files(reference_rds)
reference_object <- readRDS(reference_rds)
if (!inherits(reference_object, "Seurat")) {
  stop("mouse_HSPC_annotated.rds is not a Seurat object.")
}
if (!"celltype" %in% colnames(reference_object@meta.data)) {
  stop("mouse_HSPC_annotated.rds does not contain the required celltype column.")
}
cell_type_by_cell <- setNames(as.character(reference_object$celltype),
                              colnames(reference_object))
if (anyNA(cell_type_by_cell) || any(!nzchar(cell_type_by_cell))) {
  stop("The single-cell reference contains missing or empty cell-type labels.")
}
type_counts <- table(cell_type_by_cell)
retained_types <- names(type_counts[type_counts > 25])
keep_cells <- names(cell_type_by_cell)[cell_type_by_cell %in% retained_types]
if (!length(keep_cells) || !HSC_name %in% retained_types) {
  stop("RCTD reference filtering removed every cell or the HSC reference type.")
}
sc_counts <- get_assay_matrix(reference_object, assay = "RNA", slot = "counts")
sc_counts <- sc_counts[, keep_cells, drop = FALSE]
if (anyNA(sc_counts@x) || any(sc_counts@x < 0)) {
  stop("Single-cell RCTD reference contains missing or negative counts.")
}
cell_types <- factor(cell_type_by_cell[keep_cells])
names(cell_types) <- keep_cells
sc_nUMI <- Matrix::colSums(sc_counts)
if (!identical(colnames(sc_counts), names(cell_types)) ||
    !identical(colnames(sc_counts), names(sc_nUMI))) {
  stop("Single-cell reference counts, annotations and nUMI are not aligned.")
}
reference <- Reference(sc_counts, cell_types, sc_nUMI)

# --- 2. per-sample RCTD deconvolution ----------------------------------------
for (s in c("ST_WT", "ST_Kras")) {
  cat("Processing", s, "\n")
  stmat <- file.path(spatial_dir, s)
  sample_out <- file.path(spatial_out, s)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)
  require_files(c(file.path(stmat, "barcodes_pos.tsv.gz"),
                  file.path(stmat, "matrix.mtx.gz"),
                  file.path(stmat, "barcodes.tsv.gz"),
                  file.path(stmat, "features.tsv.gz")))
  coords <- read.table(gzfile(file.path(stmat, "barcodes_pos.tsv.gz")),
                       sep = "\t", header = FALSE)
  colnames(coords) <- c("barcode", "x", "y")
  coords$y <- -coords$y              # y-flip applied in the original analysis
  rownames(coords) <- coords$barcode

  expr <- Read10X(stmat, gene.column = 2)
  sp_obj <- CreateSeuratObject(expr, assay = "Spatial",
                               min.cells = min_cells, min.features = min_features)
  sp_counts <- get_assay_matrix(sp_obj, assay = "Spatial", slot = "counts")
  common_barcodes <- intersect(colnames(sp_counts), rownames(coords))
  if (!length(common_barcodes)) {
    stop("No shared barcodes between the count matrix and coordinates for ", s, ".")
  }
  sp_counts <- sp_counts[, common_barcodes, drop = FALSE]
  coords <- coords[common_barcodes, c("x", "y"), drop = FALSE]
  if (!identical(colnames(sp_counts), rownames(coords))) {
    stop("Count-matrix columns and coordinate rows are not aligned for ", s, ".")
  }
  sp_nUMI   <- colSums(sp_counts)
  puck <- SpatialRNA(coords, sp_counts, sp_nUMI)

  myRCTD <- create.RCTD(puck, reference, max_cores = 8, CELL_MIN_INSTANCE = 20)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  saveRDS(myRCTD, file.path(sample_out, "RCTD.rds"))

  res_df <- myRCTD@results$results_df
  result_barcodes <- intersect(rownames(res_df), names(puck@nUMI))
  if (!length(result_barcodes)) {
    stop("RCTD results and spatial nUMI have no shared barcodes for ", s, ".")
  }
  valid_bc <- result_barcodes[
    res_df[result_barcodes, "spot_class"] != "reject" & puck@nUMI[result_barcodes] >= 1
  ]
  if (!length(valid_bc)) stop("RCTD retained no non-rejected spots for ", s, ".")
  anno_df <- puck@coords[valid_bc, ]
  anno_df$cell_type <- as.character(res_df[valid_bc, "first_type"])
  write.table(anno_df %>% rownames_to_column("barcode"),
              file.path(sample_out, "Spatial_CellType.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # --- 3. HSC neighborhood composition ---------------------------------------
  hsc_coords <- anno_df[anno_df$cell_type == HSC_name, ]
  if (!nrow(hsc_coords)) {
    stop("No spots annotated as HSC for ", s,
         "; HSC-neighbourhood composition cannot be calculated.")
  }
  # Nearest-HSC lookup avoids constructing an all-spots-by-all-spots distance
  # matrix, which is not feasible for the high-density BMKMANU S1000 output.
  nn <- RANN::nn2(data = as.matrix(hsc_coords[, c("x", "y")]),
                  query = as.matrix(anno_df[, c("x", "y")]), k = 1)
  neighbor_idx <- which(nn$nn.dists[, 1] <= radius)
  neighbor_cells <- anno_df[neighbor_idx, ]
  neighbor_cells <- neighbor_cells[neighbor_cells$cell_type != HSC_name, ]
  if (!nrow(neighbor_cells)) {
    stop("No non-HSC neighbours found within radius ", radius, " for ", s, ".")
  }
  prop_df <- neighbor_cells %>% count(cell_type) %>% mutate(prop = n / sum(n))
  write.table(prop_df, file.path(sample_out, "HSC_neighbor_cell_proportion.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
