# Draft release notes: v1.1.0

Status: development; do not cite as a released version until numerical
validation and the GitHub/Zenodo release are complete.

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

## Documentation corrections

- Spatial sequencing depth now follows the BMKMANU S1000 manufacturer wording:
  approximately 100,000 paired-end reads per 100-um spatial analysis unit and
  approximately 60-150 Gb per sample.
- GSE313878 is documented as public.
- The supported environment is Seurat 4.3.0 and GSVA 2.2.0; Seurat v5 is not
  silently accepted.

## Validation still required before release

- Restore date-converted gene symbols from the original gene-set source.
- Supply the exact Figure 2F transcription-factor target list.
- Complete the human cluster 10-19 annotation map or justify their exclusion.
- Supply and archive the exact mouse-human orthologue mapping resource.
- Record the CellPhoneDB software/database versions and BSTMatrix reference
  annotation release.
- Validate the GSE71449 event-status definition against the source table.
- Run all scripts in the locked R/Python environment and compare numerical
  outputs with the manuscript figures and supplementary tables.
