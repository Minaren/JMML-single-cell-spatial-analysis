# JMML-single-cell-spatial-analysis

Analysis code for the Treg-IL-10-CD69 axis study in Kras-driven juvenile
myelomonocytic leukemia (JMML): mouse bone-marrow HSPC and T-cell single-cell
RNA-seq, human JMML single-cell analysis and CD69high HSC survival analysis,
and BMKMANU S1000 spatial transcriptomic analysis (deconvolution, HSC niche,
cell-cell communication).

The maintained and documented analysis package is in `code_publication/`.
See `code_publication/README.md` for the script list (00-08), inputs, outputs,
and run instructions. Raw sequencing data are not included in this repository;
public data accessions are listed in `code_publication/data/README.md`.

## Release

This repository is prepared as software release `v1.0.0`. See
`RELEASE_NOTES.md`, `CODE_VALIDATION.md`, and `CITATION.cff`.

## Environment

- R >= 4.2 (the exact analysis-time R version was not retained)
- Key R packages: Seurat (v4.3.0), harmony, spacexr (RCTD), CellChat (v1.6.1), GSVA
  (v2.2.0), slingshot, mgcv, clusterProfiler, limma, org.Mm.eg.db, survival,
  survminer
- External: Cell Ranger v6.0.1; BSTMatrix v1.0 (Biomarker Technologies);
  CellPhoneDB (exact analysis-time version/database not retained)

## License

MIT.
