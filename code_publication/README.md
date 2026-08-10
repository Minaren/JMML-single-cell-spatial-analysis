# Treg–IL-10–CD69 axis in Kras-driven JMML: analysis code

## Project overview

Single-cell and spatial transcriptomic analyses of a KrasG12D/+ mouse model of
juvenile myelomonocytic leukemia (JMML) and of human JMML samples. The code
covers: (1) mouse bone-marrow HSPC single-cell analysis, (2) mouse bone-marrow
T-cell (Treg) analysis, (3) human JMML single-cell analysis and CD69high HSC
survival analysis, (4) BMKMANU S1000 spatial transcriptomic analysis
(deconvolution, HSC niche, cell-cell communication).

## Repository structure

```
code_publication/
├── scripts/
│   ├── 00_setup.R                              # paths + shared packages
│   ├── 01_mouse_HSPC_processing.R              # mouse HSPC QC/integration/annotation
│   ├── 02_mouse_annotation_crossvalidation.R   # label-transfer cross-validation (optional)
│   ├── 03_mouse_HSC_subclustering.R            # HSC sub-clusters, modules, trajectory
│   ├── 04_mouse_Tcell_analysis.R               # Treg/Tconv, DE, GO, GSVA
│   ├── 05_human_data_integration.R             # human JMML + normal BM integration
│   ├── 06_human_HSC_survival_analysis.R        # CD69high HSC signature + survival
│   ├── 07_spatial_transcriptomics_analysis.R   # RCTD deconvolution + HSC niche
│   └── 08_cell_cell_communication.R            # CellPhoneDB/CellChat visualisation
├── data/          # input data (see data/README.md; raw data not included)
├── output/        # analysis outputs (written by scripts)
├── figures/       # figures (written by scripts)
└── environment/   # packages.R (package list and version reporting)
```

## Environment

- R >= 4.2 (analysis was performed with R 4.x; R 4.6.1 used during preparation)
- Key R packages: Seurat (v4.1.0 used in the study; the code is compatible with
  Seurat v5 when `JoinLayers()` is called before layer access),
  harmony, spacexr (RCTD), CellChat (v1.6.1), GSVA (v2.2.0), slingshot, mgcv,
  clusterProfiler, limma, org.Mm.eg.db, survival, survminer
- External software:
  - Cell Ranger v6.1.1 (10x Genomics) for scRNA-seq preprocessing
  - BSTMatrix v1.0 (Biomarker Technologies) for spatial transcriptomic upstream
    processing (reads mapped to the mouse reference genome mm10; exact build to
    be confirmed)
  - CellPhoneDB (v4+; version used by the original analysis to be confirmed)
    for ligand-receptor analysis
- See environment/packages.R for the full package list. To lock a reproducible
  environment with renv: `renv::init()`, install the packages, then
  `renv::snapshot()`.

## Run guide

Run the scripts in numerical order from the package root (the `.here` file
marks the project root, so `here::here()` resolves correctly regardless of the
working directory). Each script sources `scripts/00_setup.R`.

```
# from the package root:
Rscript scripts/01_mouse_HSPC_processing.R
Rscript scripts/02_mouse_annotation_crossvalidation.R   # optional
Rscript scripts/03_mouse_HSC_subclustering.R
Rscript scripts/04_mouse_Tcell_analysis.R
Rscript scripts/05_human_data_integration.R
Rscript scripts/06_human_HSC_survival_analysis.R
Rscript scripts/07_spatial_transcriptomics_analysis.R
# 08 requires CellPhoneDB outputs (see below), then:
Rscript scripts/08_cell_cell_communication.R
```

### Required input data

1. Mouse HSPC 10x output (this study; Kras and WT) -> data/raw/mouse_HSPC/
2. Mouse T-cell 10x output (this study; CD45+CD3+ sorted, 3 mice pooled per
   genotype) -> data/raw/mouse_Tcell/
3. Human JMML and normal pediatric BM 10x output (GSE155295, GSE111895;
   local samples were named JMMLID5 and PBM2 in the original analysis) ->
   data/raw/human_JMML/ and data/raw/human_PB/
4. Healthy mouse BM reference GSE122465 (notlabel.RDS + metaInfo.txt) ->
   data/reference/GSE122465/
5. Spatial transcriptomic data (BMKMANU S1000, Biomarker Technologies; not
   public) -> data/spatial/ST_WT/, data/spatial/ST_Kras/
6. Bulk RNA-seq + clinical data GSE71449 (ids_exprs.csv, Table_S1.xlsx) ->
   data/bulk/GSE71449/

See data/README.md for details, accessions and the BMK output format.

### Cell-cell communication (CellPhoneDB)

CellPhoneDB is a Python tool and is run separately on the RCTD-annotated
spatial spots. After script 07, generate the CellPhoneDB inputs from
output/spatial/ST_<sample>/Spatial_CellType.tsv and the spot count matrix, then
run:

```
cellphonedb method statistical_analysis \
  counts_<sample>.txt meta_<sample>.txt \
  --counts-data gene_name --threshold 0.1 --iterations 1000 \
  --output-path output/cellphonedb/ST_<sample>
```

Place the resulting count_network.txt, pvalues.txt, means.txt and
significant_means.txt files into output/cellphonedb/ST_WT/ and
output/cellphonedb/ST_Kras/, then run script 08 to reproduce the network
circle plots and the HSC-centred dot plots (significance threshold P < 0.05).

## Key parameters

| Step | Parameter | Value |
|---|---|---|
| Object creation | min.cells / min.features | 3 / 200 |
| QC | nFeature_RNA range / percent.mt | 500-6,000 / < 10% |
| Variable features | method / n | vst / 2,000 (HSC: 4,000) |
| PCA | n components | 50 |
| Integration | method / dims | Harmony / 1:30 |
| Clustering (all cells) | algorithm / resolution | Louvain (1) / 1.2 |
| Clustering (HSC) | resolution / dims | 0.8 / 1:50 (mouse PCA), 1:30 (human) |
| Marker detection | log2FC / min.pct | 0.5 / 0.1 (mouse); 0.25 / 0.1 (human HSC) |
| Treg DE | log2FC / FDR | 0.25 / < 0.05 |
| GSVA | kcdf | Poisson (counts), Gaussian (ssGSEA of bulk) |
| Trajectory | tool / start cluster | Slingshot / "1" |
| RCTD | CELL_MIN_INSTANCE / cores / doublet | 20 / 8 / "doublet" |
| HSC niche | neighbourhood radius | 100 coordinate units |
| CellPhoneDB | significance threshold | P < 0.05 |
| Survival | cutoffs | median + survminer::surv_cutpoint |

## Cell-type annotation

Two annotation approaches were used in the original analysis:
1. marker-based manual cluster assignment (final annotation used in the
   manuscript; scripts 01 and 05);
2. label transfer from the GSE122465 healthy bone-marrow reference
   (cross-validation; script 02).

The manuscript reports the marker-based annotation. Note that the mouse
annotation code referenced the meta column RNA_snn_res.0.6 while clustering was
performed at resolution 1.2; confirm the correspondence if re-running.

## Notes and uncertainties

- The mouse T-cell data were generated from three mice pooled per genotype
  (one 10x library per genotype); the Kras-vs-WT comparison is therefore based
  on one library per genotype and should be treated as exploratory.
- The original human HSC GSVA block referenced workspace variables
  (HSC_counts_filtered, mouse GO:BP sets, CD69_status) that were not defined in
  the shared scripts; script 06 rebuilds them deterministically and uses human
  MSigDB GO:BP sets. Please confirm this matches the analysis as performed.
- Spatial data are BMKMANU S1000 output; the experimental protocol is described
  in the Biomarker Technologies methods document (BMKMANU S1000 Spatial
  transcriptomics Materials and method) and summarised in data/README.md.
- Uncertain parameters (CellPhoneDB version/database, CellChat exact settings,
  BSTMatrix reference build, Monocle 2 usage, pI-pC administration in the
  transplantation model) are marked with [UNCERTAIN] in the scripts and listed
  in CHANGES.md.

## Outputs

See output/README.md for the expected outputs per script.

## License / data availability

Analysis code is shared for reproducibility. Raw sequencing data of this study
are available from the corresponding authors on request; public datasets are
listed in data/README.md. Spatial transcriptomic data and images are property
of the study and Biomarker Technologies; access on request.