# CHANGES.md

## v1.0.0 release preparation (2026-08-31)

- Added version, citation, Zenodo, release-note and validation metadata.
- Corrected the public human-data mapping to GSE111895 (JMML) and GSE155259
  (normal pediatric bone marrow).
- Added generated-data mappings for GSE313553, GSE313878 and GSE313879.
- Recorded the author-confirmed 6-week age for HSPC scRNA-seq mice.
- Aligned documented versions with the manuscript: Cell Ranger v6.0.1,
  Seurat v4.3.0, GSVA v2.2.0, CellChat v1.6.1 and BSTMatrix v1.0.
- Did not alter scientific analysis logic or resolve other author-dependent items.

Summary of changes applied when converting the original analysis code
(D:\aJMML\代码\code) into this publication package. The analysis logic and
parameters were preserved; only structure, paths, comments and dead code were
modified.

## File reorganisation

| Original | Publication |
|---|---|
| 合并样本8.R (mouse HSPC + HSC + spatial sections) | 01_mouse_HSPC_processing.R, 02_mouse_annotation_crossvalidation.R, 03_mouse_HSC_subclustering.R, 07_spatial_transcriptomics_analysis.R, 08_cell_cell_communication.R |
| Treg2.R | 04_mouse_Tcell_analysis.R |
| 合并样本5_human.R | 05_human_data_integration.R, 06_human_HSC_survival_analysis.R |
| (new) | 00_setup.R, environment/packages.R, data/README.md, output/README.md, figures/README.md, .here |

## Paths

- All absolute paths (D:/aJMML/...) replaced with project-relative paths built
  from `here::here()` (see scripts/00_setup.R). The `.here` marker file at the
  package root makes path resolution independent of the working directory.
- Intermediate objects (sce_harm_tu.RData, sce_har_anno.RData, mouse_CD4Tcells
  .RData, human_HSC.RData, etc.) are now saved to and loaded from output/ in
  .rds format.

## Code cleanup

- Removed `rm(list = ls())` statements and debug/inspection prints.
- Removed commented-out dead code and duplicate plotting variants; the main
  figure-generation calls were retained.
- Kept informative comments; added a header to each script (purpose, inputs,
  outputs, run order, prerequisites).
- Fixed obvious typos that did not change logic (e.g., variable name
  `plot.featrures` -> `plot.features`).

## Logic preserved / notes

- Mouse annotation: marker-based manual cluster assignment (final); the
  GSE122465 label-transfer was moved to a separate optional script (02) for
  cross-validation. See README "Cell-type annotation".
- T-cell data: three mice pooled per genotype, one library per genotype.
  Script 04 rebuilds the CD4+ T-cell object from raw data with the HSPC QC
  parameters when available; otherwise it loads a pre-processed object from
  output/. [UNCERTAIN] confirm the original T-cell QC parameters.
- Human HSC GSVA: the original block referenced undefined workspace variables
  (HSC_counts_filtered, mouse GO:BP sets, CD69_status). Script 06 builds these
  deterministically and uses human MSigDB C5 GO:BP sets. [UNCERTAIN] confirm
  that this matches the analysis as performed.
- CD69high/CD69low classification in human HSC: cluster 0 = CD69high,
  cluster 1 = CD69low (as in the original code comment). Confirm if the
  manuscript uses a different definition (e.g., expression percentile).
- Human data integration: original code read local 10x directories JMMLID5 and
  PBM2. README/data maps these to GSE111895 (JMML) and GSE155259 (normal
  pediatric BM), consistent with the manuscript and public GEO records.
- Spatial: y-coordinate flip (`coords$y <- -coords$y`) preserved from the
  original code. BMKMANU S1000 output format documented in data/README.md.
- CellPhoneDB: the Python analysis is not part of the R package; the
  visualisation script (08) reads CellPhoneDB outputs. Version and database of
  the original run are to be confirmed. CellChat settings: defaults of v1.6.1
  assumed. [UNCERTAIN]
- Monocle 2: the original manuscript methods mention Monocle 2 v2.36.0 for
  pseudotime, but the shared code uses Slingshot for the HSC trajectory;
  script 03 uses Slingshot. [UNCERTAIN] confirm whether Monocle 2 should be
  retained.
- Mouse annotation code referenced meta column RNA_snn_res.0.6 while clustering
  was run at resolution 1.2; the discrepancy is flagged in script 01.
  [UNCERTAIN]

## Items awaiting author confirmation

1. CellPhoneDB version and database used for the spatial analysis.
2. CellChat parameters (nboot, type) used in the original run.
3. BSTMatrix reference genome build (mouse mm10 exact version).
4. Whether Monocle 2 pseudotime should be retained.
5. Human HSC GSVA variable construction (see above).
6. CD69high definition in human HSC (cluster-based vs percentile).
7. T-cell QC parameters and whether each genotype was one pooled library.
8. pI-pC administration in the bone-marrow transplantation model.
