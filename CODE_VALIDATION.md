# Code validation report for v1.1.0-dev

## Corrections checked statically

- Seurat v5-only `layer=` and `JoinLayers()` calls were removed; the workflow
  now requires Seurat >=4.3.0 and <5.0.0 and uses `slot=` consistently.
- All GSVA calls use the parameter-object interface. GSE71449 survival scoring
  specifically calls `ssgseaParam()` followed by `gsva()`.
- Mouse HSPC annotation uses the active resolution-1.2 `seurat_clusters`
  column, and expected retained-cell counts are checked explicitly.
- Figure 8H multivariable Cox output and Supplementary Table S6 univariable
  Cox output are generated separately. Proportional-hazards diagnostics are
  exported.
- Spatial counts remain sparse; count columns and coordinate rows are
  intersected and reordered explicitly; sample output directories are created.
- CellPhoneDB inputs are prepared only after an explicit mouse-human mapping;
  metadata precede counts in the documented command, and human-symbol h5ad is
  used for the high-density spatial data.
- Python syntax, JSON/YAML metadata, Git integrity and whitespace checks are
  included in the validation record produced for this working tree.

## Validation boundary

The current execution environment does not contain R or the complete deposited
input datasets. R parsing, package installation, end-to-end execution and
numerical comparison with manuscript panels have therefore not yet been
performed. This development version must not be tagged or archived as v1.1.0
until those checks pass.

## Inputs/provenance still required

1. Restore spreadsheet-corrupted gene symbols in the mouse gene-set files from
   the original source.
2. Supply the exact transcription-factor target lists used for Figure 2F as
   `data/gene_sets/tf_targets_mouse.csv`.
3. Complete or formally exclude human clusters 10-19 after reviewing the
   exported marker table.
4. Supply the exact mouse-human orthologue table and archive its database name,
   version, download date and URL.
5. Record the CellPhoneDB software/database versions and the exact BSTMatrix
   mm10 annotation release.
6. Validate the GSE71449 event-status rule against the source clinical table.
7. Provide `sample_metadata.tsv` for the three NC and three OE GSE313879
   libraries and compare DESeq2 results with the manuscript bulk-RNA outputs.
8. Confirm whether the final trajectory figure was generated with Slingshot or
   Monocle 2; the current script implements Slingshot.
