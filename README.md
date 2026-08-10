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

## Environment

- R >= 4.2 (prepared with R 4.6.1)
- Key R packages: Seurat, harmony, spacexr (RCTD), CellChat (v1.6.1), GSVA
  (v2.2.0), slingshot, mgcv, clusterProfiler, limma, org.Mm.eg.db, survival,
  survminer
- External: Cell Ranger v6.1.1; CellPhoneDB (v5.0.1; cellphonedb-data v5.0.0);
  BSTMatrix v1.0 (Biomarker Technologies)

## License

MIT.