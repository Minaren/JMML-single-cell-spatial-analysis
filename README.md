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

The audited `v1.1.0` corrective release updates the analysis code, provenance
records, and automated validation package. It does not replace or rewrite the
archived v1.0.0 DOI record. The precise validation boundary is disclosed in
`RELEASE_NOTES_v1.1.0.md` and `CODE_VALIDATION.md`.

Suggested citation:

> Ren X, Li Q, Yue J, Zhang L, He A, Kong G. JMML single-cell and spatial
> transcriptomic analysis code. Version 1.1.0. Zenodo; 2026.

The version-specific v1.1.0 DOI is added after the GitHub-Zenodo archive has
finished processing the release.

## Environment

- R >= 4.5 with Bioconductor 3.21 (the exact patch version was not retained)
- Key R packages: Seurat (v4.3.0), harmony, spacexr (RCTD), GSVA
  (v2.2.0), slingshot, mgcv, clusterProfiler, DESeq2, org.Mm.eg.db, survival,
  survminer
- External: bcl2fastq2 v2.20.0; Cell Ranger v7.0.0; BSTMatrix v1.0
  (Biomarker Technologies); CellPhoneDB v5.0.1 with cellphonedb-data v5.0.0

## License

MIT.
