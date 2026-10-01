# Code validation report for v1.1.0

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
- Author-supplied WT/Kras spatial matrices, coordinates, read tables and H&E
  TIFFs pass the dedicated input validator; dimensions and SHA-256 checksums
  are recorded in `code_publication/data/spatial/INPUT_PROVENANCE.md`.
- Script 07 now builds the RCTD reference from script 01's annotated Seurat
  object rather than unversioned dense `sc_meta.txt/ref_cell_anno` exports.
- CellPhoneDB inputs are prepared only after an explicit mouse-human mapping;
  metadata precede counts in the documented command, and human-symbol h5ad is
  used for the high-density spatial data.
- Python syntax, JSON/YAML metadata, Git integrity and whitespace checks are
  included in the validation record produced for this working tree.
- A synthetic regression test verifies sparse mouse-to-human count aggregation,
  exclusion of ambiguous one-to-many orthologues, barcode sanitisation and
  h5ad/metadata output with the locked CellPhoneDB Python environment.
- The full human cluster map and compact TF modules were recovered from initial
  commit `01da390`; repaired gene-set values are documented in
  `code_publication/data/gene_sets/GENE_SET_REPAIR_REQUIRED.md`.
- GSE71449 overall-survival events are defined by non-empty Table S1
  `Cause of death` values, matching the retained cohort analysis record.
- The corrected CellPhoneDB rerun uses the included MGI mouse-human homology
  report mapping downloaded on 2026-09-30; source/output checksums are recorded.
- GSE313879 sample IDs and NC/OE assignments were verified against the public
  GEO family record and are included in the repository.

## Validation boundary

GitHub Actions parses all maintained R scripts under R 4.5.1 and runs the
dependency-free publication-package checks plus the locked CellPhoneDB input
regression test. The CI environment does not contain the complete deposited
single-cell inputs or the full specialist R package stack, so end-to-end
execution and numerical comparison with every manuscript panel have not been
performed. Version 1.1.0 is therefore an audited corrective code archive, not
a claim of complete numerical reproduction.
