# Code validation report for v1.0.0

## Completed

- GEO mappings and the author-confirmed 6-week HSPC scRNA-seq mouse age were
  aligned with the manuscript and public GSE313553 record.
- Documented Cell Ranger, Seurat, GSVA, CellChat and BSTMatrix versions were
  aligned with the manuscript.
- Scripts use project-relative paths and document inputs, outputs and run order.
- Release metadata include VERSION, CITATION.cff, .zenodo.json, MIT LICENSE and
  release notes.

## Validation boundary

The nine R scripts passed static delimiter and UTF-8 checks. End-to-end execution
and numerical figure comparison were not possible without the complete inputs
and matching R environment.

## Author confirmation still required

1. The exact clustering column used for final mouse annotation
   (`RNA_snn_res.0.6` versus the resolution-1.2 clustering).
2. Exact transcription-factor target lists used for Figure 2F.
3. Whether Slingshot or Monocle 2 generated the final pseudotime figure.
4. Original CellPhoneDB software/database versions.
5. Non-default CellChat parameters, if any.
6. Exact BSTMatrix mm10 reference build and annotation release.
7. Human HSC GSVA construction and CD69high/CD69low definition.
8. T-cell QC parameters and pooled-library design.

Until these items are resolved, describe this as the documented analysis code,
not an independently verified exact reproduction package.
