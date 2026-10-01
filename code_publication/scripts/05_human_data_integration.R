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
# Full resolution-1.2 map restored from the retained manuscript-generation
# script 投稿2.0.R (initial repository commit 01da390, lines 756-768).
celltype[celltype$ClusterID == 0, "celltype"] <- "MPP"
celltype[celltype$ClusterID == 1, "celltype"] <- "CMP"
celltype[celltype$ClusterID %in% c(2, 6, 16), "celltype"] <- "Monocyte"
celltype[celltype$ClusterID %in% c(3, 5, 10, 15), "celltype"] <- "Granulocyte"
celltype[celltype$ClusterID %in% c(4, 17), "celltype"] <- "CLP"
celltype[celltype$ClusterID == 8, "celltype"] <- "LMPP"
celltype[celltype$ClusterID == 9, "celltype"] <- "HSC"
celltype[celltype$ClusterID == 11, "celltype"] <- "proB"
celltype[celltype$ClusterID %in% c(7, 12, 13, 19), "celltype"] <- "Erythroid_Progenitor"
celltype[celltype$ClusterID == 14, "celltype"] <- "Macrophage"
celltype[celltype$ClusterID == 18, "celltype"] <- "MK"
if (any(grepl("^Unassigned_", celltype$celltype))) {
  stop("The restored human cluster map does not cover every observed cluster: ",
       paste(celltype$ClusterID[grepl("^Unassigned_", celltype$celltype)], collapse = ", "))
}

sce.harm$celltype <- NA_character_
for (i in seq_len(nrow(celltype))) {
  sce.harm@meta.data[which(as.character(sce.harm@meta.data$seurat_clusters) ==
                             as.character(celltype$ClusterID[i])),
                     "celltype"] <- celltype$celltype[i]
}
if (anyNA(sce.harm$celltype)) stop("At least one human cell did not receive a cell-type label.")
Idents(sce.harm) <- "celltype"
saveRDS(sce.harm, file.path(output_dir, "human_bone_marrow.rds"))

write.csv(celltype, file.path(output_dir, "human_cluster_annotation_map.csv"),
          row.names = FALSE)

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
