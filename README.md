# JMML-single-cell-spatial-analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22207699.svg)](https://doi.org/10.5281/zenodo.22207699)

Analysis code for the Treg-IL-10-CD69 axis study in Kras-driven juvenile
myelomonocytic leukemia (JMML): mouse bone-marrow HSPC and T-cell single-cell
RNA-seq, human JMML single-cell analysis and CD69high HSC survival analysis,
and BMKMANU S1000 spatial transcriptomic analysis (deconvolution, HSC niche,
cell-cell communication).

The maintained and documented analysis package is in `code_publication/`.
See `code_publication/README.md` for the script list (00-09), inputs, outputs,
and run instructions. Raw sequencing data are not included in this repository;
public data accessions are listed in `code_publication/data/README.md`.

## Release

The DOI-linked software release `v1.0.0` remains available from Zenodo at
<https://doi.org/10.5281/zenodo.22207699>.

The audited `v1.1.0-dev` correction set is being prepared on a development
branch. It must be numerically rerun against the deposited inputs before being
merged and tagged as `v1.1.0`; it does not replace or rewrite the archived
v1.0.0 DOI record. See `RELEASE_NOTES_v1.1.0-draft.md` and
`CODE_VALIDATION.md`.

Suggested citation:

> Ren X, Li Q, Yue J, Zhang L, He A, Kong G. JMML single-cell and spatial
> transcriptomic analysis code. Version 1.0.0. Zenodo; 2026.
> https://doi.org/10.5281/zenodo.22207699

## Environment

- R >= 4.5 with Bioconductor 3.21 (the exact patch version was not retained)
- Key R packages: Seurat (v4.3.0), harmony, spacexr (RCTD), GSVA
  (v2.2.0), slingshot, mgcv, clusterProfiler, DESeq2, org.Mm.eg.db, survival,
  survminer
- External: Cell Ranger v7.0.0 as confirmed by the author (the GEO record
  currently says v2.1.1 and must be reconciled); BSTMatrix v1.0
  (Biomarker Technologies); CellPhoneDB v5.0.1 with cellphonedb-data v5.0.0

## License

MIT.
