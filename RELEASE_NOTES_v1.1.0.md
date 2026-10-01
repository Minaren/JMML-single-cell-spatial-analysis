# Release notes: v1.1.0

Released 2026-10-01. This version archives the audited corrective code set.
The validation scope is stated below so that the release does not imply an
unperformed end-to-end numerical rerun.

## Corrected analysis code

- Replaced obsolete GSVA calls with the parameter-object API required by
  GSVA v2.2.0.
- Implemented the manuscript-specified survival signature with
  `ssgseaParam()` followed by `gsva()`.
- Removed Seurat v5 layer calls and `JoinLayers()` from the Seurat v4.3.0
  workflow.
- Applied the mouse cell-type map to the active resolution-1.2 clustering.
- Compared human HSC cluster 0 directly with cluster 1 when constructing the
  CD69-high and CD69-low signatures.
- Added separate univariable Cox output for Supplementary Table S6 and
  multivariable Cox plus proportional-hazards diagnostics for Figure 8H.
- Replaced cell-as-replicate pathway P values in pooled-library comparisons
  with clearly labelled descriptive scores.
- Kept BMKMANU S1000 matrices sparse, aligned barcodes explicitly and replaced
  the all-by-all spatial distance matrix with nearest-neighbour search.
- Added versionable mouse-to-human orthologue conversion and sparse h5ad input
  preparation for CellPhoneDB; corrected CellPhoneDB argument order.
- Rebuilt the CellPhoneDB dot-plot parser so metadata columns are not treated
  as interaction measurements.
- Removed the non-executable CellChat placeholder instead of implying that an
  undocumented default-parameter analysis was reproducible.
- Added early checks for retained cell counts, missing inputs, gene-set overlap,
  missing annotations and spreadsheet-corrupted gene symbols.
- Restored the complete human cluster 0-19 annotation map and the original
  Nfkb1/Fos/Stat3/Jun target modules from commit `01da390`.
- Repaired date-converted and I/l-corrupted mouse gene symbols, with an
  auditable repair record.
- Fixed the GSE71449 event definition to the manuscript cohort rule based on
  the Table S1 `Cause of death` field.
- Added the GEO-verified six-sample GSE313879 metadata and a versioned MGI
  mouse-human orthologue map with checksums.
- Added a fully pinned CellPhoneDB v5.0.1 Python lock file and an automated
  sparse orthologue-conversion/h5ad regression test.
- Replaced the unversioned dense RCTD reference exports with the annotated
  script-01 Seurat object and added validation/checksums for author-supplied
  BMKMANU coordinates and H&E images.

## Documentation corrections

- Spatial sequencing depth now follows the BMKMANU S1000 manufacturer wording:
  approximately 100,000 paired-end reads per 100-um spatial analysis unit and
  approximately 60-150 Gb per sample.
- GSE313878 is documented as public.
- The supported environment is Seurat 4.3.0 and GSVA 2.2.0; Seurat v5 is not
  silently accepted.

## Validation scope

GitHub Actions passed the publication-package checks, Python regression test
for CellPhoneDB input preparation, and R syntax parsing under R 4.5.1. The
author-supplied spatial matrices, coordinates, read tables, and H&E images
also passed the repository input-integrity validator.

A complete end-to-end numerical rerun against every manuscript figure and
supplementary table was not performed in the release CI because the full
deposited single-cell inputs and specialist R environment are not bundled with
the repository.
