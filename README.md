# JMML single-cell & spatial transcriptomics analysis code

Analysis code for the study of a Treg–IL-10–CD69 axis in Kras-driven juvenile
myelomonocytic leukemia (JMML), covering:

1. Mouse bone-marrow HSPC single-cell RNA-seq analysis (QC, Harmony
   integration, clustering, annotation, HSC sub-clustering, pseudotime)
2. Mouse bone-marrow T-cell (Treg/Tconv) analysis (DE, GO, GSVA)
3. Human JMML single-cell analysis and CD69high HSC survival analysis
4. BMKMANU S1000 spatial transcriptomic analysis (RCTD deconvolution,
   HSC niche, CellPhoneDB/CellChat cell-cell communication)

## Repository layout

| Path | Content |
|---|---|
| `code_publication/` | **Maintained release package** (numbered, documented scripts with relative paths). Start here. See `code_publication/README.md`. |
| `code/` | Original working scripts used for the study (`Treg2.R`, `合并样本5_human.R`, `合并样本8.R`). |
| `mouse/` `human/` `human_RNAseq/` `空间转录组/` `备用/` | Original working scripts and notebooks (absolute paths; kept for reference). |
| Root `*.R` | Miscellaneous working scripts. |

The `code_publication/` package is the recommended entry point: it contains a
setup script, per-analysis scripts (01–08), an environment package list, and
documentation of inputs/outputs. The scripts in the other directories are the
original working versions and contain absolute paths (e.g., `D:/aJMML/...`);
adjust these paths before running.

## Environment

- R >= 4.2 (prepared with R 4.6.1)
- Key R packages: Seurat (v4.1.0 in the study; compatible with v5 when using
  `JoinLayers()`), harmony, spacexr (RCTD), CellChat (v1.6.1), GSVA (v2.2.0),
  slingshot, mgcv, clusterProfiler, limma, org.Mm.eg.db, survival, survminer
- External: Cell Ranger v6.1.1; CellPhoneDB (v5.0.1; database cellphonedb-data
  v5.0.0) for ligand–receptor analysis; BSTMatrix v1.0 (Biomarker
  Technologies) for spatial transcriptomic upstream processing

## Data availability

Raw data are not included in this repository. Public data:

| Dataset | Accession |
|---|---|
| JMML bone marrow scRNA-seq | GSE111895 |
| Normal pediatric bone marrow scRNA-seq | GSE155259 |
| Healthy mouse bone-marrow reference | GSE122465 |
| Bulk RNA-seq (survival cohort) | GSE71449 |

Mouse 10x scRNA-seq, mouse T-cell scRNA-seq, and BMKMANU S1000 spatial
transcriptomic data generated for this study are available from the Gene
Expression Omnibus (accessions listed in the manuscript). See
`code_publication/data/README.md` for the expected input layout.

## License

MIT. See `LICENSE`.

## Citation

If you use this code, please cite the associated manuscript (details to be
added upon publication).