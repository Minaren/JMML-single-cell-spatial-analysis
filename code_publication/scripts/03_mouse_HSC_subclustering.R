# ============================================================================
# Script: 03_mouse_HSC_subclustering.R
# Purpose: Re-analysis of the HSC compartment: sub-clustering (resolution 0.8),
#          functional module and transcription-factor activity scoring
#          (AddModuleScore), cell-surface marker module mapping, Slingshot
#          pseudotime with GAM-smoothed trends, and module correlations.
# Inputs:  output/mouse_HSPC_annotated.rds (annotated mouse HSPC object)
# Outputs: output/mouse_HSC.rds (HSC subset with scores and trajectory)
#          figures/Fig_HSC_*.pdf
# Run order: after 01_mouse_HSPC_processing.R
# NOTE: CD69-high vs CD69-low HSC classification is performed downstream by
#       flow cytometry (Methods); here the three transcriptomic HSC clusters
#       are characterised.
# ============================================================================

source("scripts/00_setup.R")

sce.harm <- readRDS(file.path(output_dir, "mouse_HSPC_annotated.rds"))

# --- 1. HSC sub-clustering --------------------------------------------------
HSC <- subset(sce.harm, celltype == "HSC")
HSC <- ScaleData(HSC, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
HSC <- FindVariableFeatures(HSC, nfeatures = 4000)
HSC <- RunPCA(HSC, npcs = 50, verbose = FALSE)
HSC <- FindNeighbors(HSC, reduction = "pca", dims = 1:50)
HSC <- FindClusters(HSC, resolution = 0.8)
HSC <- RunUMAP(HSC, reduction = "pca", dims = 1:50)
Idents(HSC) <- "seurat_clusters"

# --- 2. functional module scoring (Hallmark-style curated gene sets) --------
inflam_genes <- c("CCL20","CCL4","CCL5","CCL7","CCL8","CD69","CD74","CSF1","CSF2","CSF3",
                  "CXCL1","CXCL10","CXCL2","CXCL3","CXCL8","ICAM1","IL12B","IL15","IL18",
                  "IL1A","IL1B","IL6","IRF7","LTA","PTGS2","SELE","SELP","TNF","VCAM1")
ifna_genes <- c("BST2","EIF2AK2","GBP1","GBP2","IFI27","IFI35","IFI6","IFIT1","IFIT2",
                "IFIT3","IFITM1","IFITM2","IFITM3","IRF1","IRF7","IRF9","ISG15","MX1",
                "MX2","OAS1","OAS2","OAS3","RSAD2","STAT1","STAT2","USP18")
ifng_genes <- c("CASP1","CCL2","CCL5","CCL7","CD274","CIITA","CXCL10","CXCL9","GBP1",
                "GBP2","GBP4","GBP5","ICAM1","IFNG","IFNGR1","IFNGR2","IRF1","IRF7",
                "IRF9","ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","RSAD2","SOCS1",
                "STAT1","STAT2","TAP1","TAP2","TNF","USP18","VCAM1")
kras_genes <- c("ADAMTS1","AKT1","APAF1","AREG","ATF3","BCL2","BMP2","BTG2","CCL2",
                "CCL20","CCL4","CCL5","CCND1","CD44","CDKN1A","CFLAR","CLU","CSF1",
                "CSF2","CSF3","CTGF","CXCL1","CXCL10","CXCL2","CXCL3","CXCL8","CYR61",
                "DUSP1","DUSP4","DUSP5","DUSP6","EDN1","EFNB1","EGR1","EGR2","EGR3",
                "EIF4EBP1","EPHA2","EREG","ETV5","FGF2","FOS","FOSB","FOSL1","FOSL2",
                "FST","GADD45A","GADD45B","GAPDH","GEM","HBEGF","HMGA2","HMOX1",
                "ICAM1","IER2","IGF1R","IL13RA1","IL1A","IL1B","IL6","IL7R","IRF1",
                "JUN","JUNB","KDM5B","KLF2","KLF4","KLF6","LIF","MAP2K3","MAPK8",
                "MCL1","MEST","MUC1","MYC","NFKB1","NFKB2","NR4A1","NR4A2","NR4A3",
                "PDGFB","PLAUR","PLK2","PLK3","PMEPA1","PPP1R15A","PTGER4","PTGS2",
                "PTPRE","RAC1","RAP1A","REL","RELB","RHOB","SAT1","SERPINB2","SERPINE1",
                "SKP2","SLC2A1","SLC2A3","SLC4A7","SNHG12","SPRY1","SPRY2","SPRY4",
                "SYK","TBX2","TFPI2","TGFBR2","THBS1","TIMP1","TNF","TRIB3","UPP1",
                "VEGFA","VIM","WT1")
myeloid_genes <- c("AIF1","ASXL1","BCL11A","BMP4","CCR1","CCR2","CD34","CEBPA",
                   "CEBPB","CEBPE","CLEC4D","CSF1","CSF1R","CSF2","CSF2RA","CSF2RB",
                   "CSF3","CSF3R","CTSG","ELANE","FCER2","FLT3","GATA1","GATA2",
                   "GFI1","HOXB4","ID2","IFNG","IKZF1","IL3","IL5","IRF8","ITGAM",
                   "JAG1","KIT","KITLG","LY6G","LYL1","MAFB","MEIS1","MPO","MYB",
                   "NFE2","NOTCH1","PRTN3","PU1","RUNX1","S100A8","S100A9","SPI1",
                   "STAT3","STAT5A","STAT5B","TAL1","TLR2","TNF","ZFP36")
stemness_genes <- c("Procr", "Ly6a", "Mllt3", "Mecom", "Hlf", "Gfi1",
                    "Tek", "Fgd5", "Hoxb5", "Mpl", "Cd34", "Kit",
                    "Hoxa9", "Meis1", "Fli1", "Gata2", "Tcf15", "Cdkn1c")
module_list <- list(
  Stemness              = stemness_genes,
  Myeloid_Priming       = myeloid_genes,
  Inflammatory_Response = inflam_genes,
  IFN_Response          = unique(c(ifna_genes, ifng_genes)),
  KRAS_Signaling        = kras_genes
)

old_cols <- grep("_Score$", colnames(HSC@meta.data), value = TRUE)
if (length(old_cols)) HSC@meta.data <- HSC@meta.data[, !colnames(HSC@meta.data) %in% old_cols]

for (nm in names(module_list)) {
  HSC <- AddModuleScore(HSC, features = module_list[nm],
                        name = paste0(nm, "_Score"), ctrl = min(100, nrow(HSC)))
  colnames(HSC@meta.data)[ncol(HSC@meta.data)] <- paste0(nm, "_Score")
}

# --- 3. transcription-factor target scoring ---------------------------------
tf_targets <- list(
  Nfkb1 = c("Ccl2","Ccl5","Cxcl1","Cxcl2","Il6","Tnf","Icam1","Nfkbia","Ccl20","Cxcl10"),
  Fos   = c("Fos","Fosb","Egr1","Junb","Dusp1","Dusp5","Dusp6","Nr4a1"),
  Jun   = c("Jun","Junb","Jund","Fos","Fosb","Ccl2","Cd44","Vegfa","Timp1"),
  Stat3 = c("Socs3","Il6","Ccl2","Ccl5","Mcl1","Bcl2l1","Cebpb","Il1b")
)
for (tf in names(tf_targets)) {
  HSC <- AddModuleScore(HSC, features = tf_targets[tf],
                        name = paste0(tf, "_Score"), ctrl = min(100, nrow(HSC)))
  colnames(HSC@meta.data)[ncol(HSC@meta.data)] <- paste0(tf, "_Score")
}
# NOTE: TF target gene lists above are examples; the original script read
#       tf_targets from a workspace object. Confirm the exact target lists if
#       exact reproduction of Fig. 2F is required.

# --- 4. Slingshot trajectory and GAM-smoothed trends ------------------------
# Embedding from the HSC UMAP; start cluster set to "1" as in the original.
rd <- Embeddings(HSC, "umap")
cl <- as.character(Idents(HSC))
sds <- slingshot(rd, cl, start.clus = "1")
HSC$pseudotime <- sds@curves[[1]]$lambda[rownames(Embeddings(HSC, "umap"))]

# GAM fit of module scores and Cd69 against pseudotime (per-cluster trends)
pseudo_df <- data.frame(
  pseudotime = HSC$pseudotime,
  cluster    = Idents(HSC),
  Cd69       = GetAssayData(HSC, assay = "RNA", layer = "data")["Cd69", ],
  FetchData(HSC, vars = paste0(names(module_list), "_Score"))
)
p_dynamics <- pseudo_df %>%
  pivot_longer(cols = c("Cd69", paste0(names(module_list), "_Score")),
               names_to = "module", values_to = "value") %>%
  ggplot(aes(x = pseudotime, y = value)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.1, alpha = 0.15) +
  facet_wrap(~ module, scales = "free_y", ncol = 2) +
  theme_bw()
ggsave(file.path(fig_dir, "Fig_Dynamics_Modules_Cd69.pdf"), p_dynamics,
       width = 7, height = 9)

# --- 5. module correlation (Spearman) ---------------------------------------
cor_data <- FetchData(HSC, vars = c(paste0(names(module_list), "_Score"),
                                    paste0(names(tf_targets), "_Score")))
cor_mat <- cor(cor_data, method = "spearman")
write.csv(cor_mat, file.path(output_dir, "HSC_module_correlation.csv"))

saveRDS(HSC, file.path(output_dir, "mouse_HSC.rds"))