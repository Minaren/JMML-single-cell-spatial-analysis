# ============================================================================
# Script: 02_mouse_annotation_crossvalidation.R
# Purpose: Optional cross-validation of the marker-based cell-type annotation
#          (script 01) by label transfer from the healthy mouse bone-marrow
#          reference GSE122465 (Baryawno et al., Cell 2019; scmap-based cell
#          type labels provided by the data source).
# Inputs:  output/mouse_HSPC_annotated.rds   (annotated mouse HSPC object)
#          data/reference/GSE122465/notlabel.RDS
#          data/reference/GSE122465/metaInfo.txt
# Outputs: output/mouse_label_transfer_predictions.rds (prediction scores)
#          figures/umap_label_transfer.pdf
# Run order: after 01_mouse_HSPC_processing.R (optional)
# NOTE: The final annotation used in the manuscript is the marker-based one
#       from script 01. This label-transfer result was used only as a
#       cross-validation; the two approaches are compared in the README.
# ============================================================================

source("scripts/00_setup.R")

# --- 1. build the reference -------------------------------------------------
sc <- readRDS(file.path(ref_dir, "GSE122465", "notlabel.RDS"))
sc <- UpdateSeuratObject(sc)

meta_info <- read.delim(file.path(ref_dir, "GSE122465", "metaInfo.txt"),
                        header = TRUE, row.names = 1, sep = "\t")
sc@meta.data$CellType <- meta_info[rownames(sc@meta.data), "CellType..based.on.scmap."]

# --- 2. label transfer ------------------------------------------------------
query <- readRDS(file.path(output_dir, "mouse_HSPC_annotated.rds"))
anchors <- FindTransferAnchors(reference = sc, query = query, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = sc$CellType, dims = 1:30)
query$predicted.id <- predictions$predicted.id

saveRDS(predictions, file.path(output_dir, "mouse_label_transfer_predictions.rds"))

p <- DimPlot(query, reduction = "umap", group.by = "predicted.id", label = TRUE) +
  ggtitle("Label transfer from GSE122465")
ggsave(file.path(fig_dir, "umap_label_transfer.pdf"), p, width = 12, height = 10)