# ============================================================================
# Script: 05_human_data_integration.R
# Purpose: Human JMML vs normal bone-marrow scRNA-seq integration and
#          analysis: QC, Harmony integration, clustering (resolution 1.2),
#          UMAP, marker-based cell-type annotation, GSVA validation, and
#          cell-type proportions.
# Inputs:  data/raw/human_JMML/ and data/raw/human_PB/
#          (10x Cell Ranger output directories; JMML patient and normal
#          pediatric bone-marrow samples; see data/README.md for accessions)
#          data/gene_sets/gsva_human_cluster.csv
# Outputs: output/human_bone_marrow.rds (integrated, annotated object)
#          output/gsva_res_human.csv
# Run order: after 00_setup.R
# NOTE: the original code read two local 10x directories (JMMLID5, PBM2).
#       Confirm the correspondence between these local samples and the public
#       accessions GSE111895 (JMML) and GSE155259 (normal pediatric BM).
# ============================================================================

source("scripts/00_setup.R")

# --- 1. read and merge samples ---------------------------------------------
dirs <- c(file.path(raw_dir, "human_JMML"),
          file.path(raw_dir, "human_PB"))
samples_name <- c("JMML", "Healthy")

scRNAlist <- list()
for (i in seq_along(dirs)) {
  counts <- Read10X(data.dir = dirs[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                       min.cells = 3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
  scRNAlist[[i]][["percent.MT"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
}
sc_merge <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

# --- 2. QC, integration, clustering -----------------------------------------
sc_filt <- subset(sc_merge, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.MT < 10)
sc_n <- NormalizeData(sc_filt)
sc_n <- FindVariableFeatures(sc_n, selection.method = "vst", nfeatures = 2000)
sc_n <- ScaleData(sc_n)
sc_n <- RunPCA(sc_n, npcs = 50)
sc_n <- RunUMAP(sc_n, dims = 1:30, reduction = "pca",
                reduction.name = "umap_pca", reduction.key = "UMAPpca_")
p_no_harmony <- DimPlot(sc_n, reduction = "umap_pca", group.by = "orig.ident") +
  ggtitle("Human bone marrow before Harmony (condition-confounded samples)")
ggsave(file.path(fig_dir, "sensitivity_human_BM_preHarmony.pdf"),
       p_no_harmony, width = 7, height = 6)
sce_har <- RunHarmony(sc_n, group.by.vars = "orig.ident")
sce.har <- FindNeighbors(sce_har, dims = 1:30, reduction = "harmony")
sce.har <- FindClusters(sce.har, resolution = 1.2, algorithm = 1)
sce.harm <- RunUMAP(sce.har, dims = 1:30, reduction = "harmony")

# --- 3. marker-based cell-type annotation -----------------------------------
# Manual cluster->celltype mapping retained from the original analysis.
observed_cluster_ids <- sort(unique(as.integer(as.character(sce.harm$seurat_clusters))))
celltype <- data.frame(ClusterID = observed_cluster_ids,
                       celltype = paste0("Unassigned_", observed_cluster_ids))
celltype[celltype$ClusterID == 9, "celltype"] <- "HSC"
celltype[celltype$ClusterID == 0, "celltype"] <- "MPP"
celltype[celltype$ClusterID == 1, "celltype"] <- "CMP"
celltype[celltype$ClusterID == 2, "celltype"] <- "GMP"
celltype[celltype$ClusterID == 3, "celltype"] <- "MEP"
celltype[celltype$ClusterID == 4, "celltype"] <- "CLP"
celltype[celltype$ClusterID == 5, "celltype"] <- "B"
celltype[celltype$ClusterID == 6, "celltype"] <- "T"
celltype[celltype$ClusterID == 7, "celltype"] <- "NK"
celltype[celltype$ClusterID == 8, "celltype"] <- "Mono"
# Clusters 10-19 were not assigned in the retained original code. They are
# kept explicitly as Unassigned_<cluster> rather than converted to NA or
# silently discarded. Their markers are exported for author review.

sce.harm$celltype <- NA_character_
for (i in seq_len(nrow(celltype))) {
  sce.harm@meta.data[which(sce.harm@meta.data$seurat_clusters == celltype$ClusterID[i]),
                     "celltype"] <- celltype$celltype[i]
}
Idents(sce.harm) <- "celltype"
saveRDS(sce.harm, file.path(output_dir, "human_bone_marrow.rds"))

unassigned_cells <- rownames(sce.harm@meta.data)[grepl("^Unassigned_", sce.harm$celltype)]
if (length(unassigned_cells)) {
  unassigned <- subset(sce.harm, cells = unassigned_cells)
  Idents(unassigned) <- "seurat_clusters"
  unassigned_markers <- FindAllMarkers(unassigned, only.pos = TRUE,
                                       logfc.threshold = 0.25, min.pct = 0.1)
  write.csv(unassigned_markers,
            file.path(output_dir, "human_unassigned_cluster_markers.csv"),
            row.names = FALSE)
}

# --- 4. GSVA validation of annotations --------------------------------------
genesets <- read_gene_sets(file.path(gene_dir, "gsva_human_cluster.csv"))
expr <- AverageExpression(sce.harm, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr) > 0, ]
expr <- as.matrix(expr)
gsva.res <- run_gsva(expr, genesets, kcdf = "Gaussian")
write.csv(data.frame(Genesets = rownames(gsva.res), gsva.res, check.names = FALSE),
          file.path(output_dir, "gsva_res_human.csv"), row.names = FALSE)

# --- 5. cell-type proportions ------------------------------------------------
Cellratio <- prop.table(table(Idents(sce.harm), sce.harm$orig.ident), margin = 2)
write.csv(as.data.frame(Cellratio), file.path(output_dir, "human_celltype_proportions.csv"))
