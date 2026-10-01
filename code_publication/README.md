# Treg–IL-10–CD69 axis in Kras-driven JMML: analysis code

## Project overview

Single-cell and spatial transcriptomic analyses of a KrasG12D/+ mouse model of
juvenile myelomonocytic leukemia (JMML) and of human JMML samples. The code
covers: (1) mouse bone-marrow HSPC single-cell analysis, (2) mouse bone-marrow
T-cell (Treg) analysis, (3) human JMML single-cell analysis and CD69high HSC
survival analysis, (4) BMKMANU S1000 spatial transcriptomic analysis
(deconvolution, HSC niche, cell-cell communication).

## Repository structure

```
code_publication/
├── scripts/
│   ├── 00_setup.R                              # paths + shared packages
│   ├── 01_mouse_HSPC_processing.R              # mouse HSPC QC/integration/annotation
│   ├── 02_mouse_annotation_crossvalidation.R   # label-transfer cross-validation (optional)
│   ├── 03_mouse_HSC_subclustering.R            # HSC sub-clusters, modules, trajectory
│   ├── 04_mouse_Tcell_analysis.R               # Treg/Tconv, DE, GO, GSVA
│   ├── 05_human_data_integration.R             # human JMML + normal BM integration
│   ├── 06_human_HSC_survival_analysis.R        # CD69high HSC signature + survival
│   ├── 07_spatial_transcriptomics_analysis.R   # RCTD deconvolution + HSC niche
│   ├── 08_prepare_cellphonedb_inputs.py         # orthologue conversion + sparse h5ad
│   ├── 08_run_cellphonedb.py                    # versioned CellPhoneDB v5 run
│   ├── 08_cell_cell_communication.R             # CellPhoneDB visualisation
│   └── 09_NB4_bulk_RNAseq.R                     # GSE313879 NB4 DESeq2 analysis
├── data/          # input data (see data/README.md; raw data not included)
├── output/        # analysis outputs (written by scripts)
├── figures/       # figures (written by scripts)
└── environment/   # packages.R (package list and version reporting)
```

## Environment

- R >= 4.5 (Bioconductor 3.21, containing GSVA v2.2.0;
  the exact analysis-time patch version was not retained)
- Key R packages: Seurat v4.3.0 (this release deliberately stops on Seurat v5
  rather than mixing v4 slots with v5 assay layers),
  harmony, spacexr (RCTD), GSVA (v2.2.0), slingshot, mgcv, igraph,
  clusterProfiler, DESeq2, org.Mm.eg.db, RANN, survival, survminer
- External software:
  - Cell Ranger v7.0.0 (10x Genomics) for scRNA-seq preprocessing, as confirmed
    by the author. GSE313553 currently reports v2.1.1; this
    external metadata conflict is a release blocker documented below.
  - BSTMatrix v1.0 (Biomarker Technologies) for spatial transcriptomic upstream
    processing (reads mapped to the mouse reference genome mm10; exact build to
    be confirmed)
  - CellPhoneDB v5.0.1 with cellphonedb-data v5.0.0; the runner records the
    installed package version and database checksum.
- See environment/packages.R for the full package list. To lock a reproducible
  environment with renv: `renv::init()`, install the packages, then
  `renv::snapshot()`.
- `environment/requirements-cellphonedb-lock.txt` records the complete tested
  Python environment for CellPhoneDB v5.0.1. Install it in a dedicated virtual
  environment with
  `python -m pip install -r environment/requirements-cellphonedb-lock.txt`.

## Run guide

Run the scripts in numerical order from the package root (the `.here` file
marks the project root, so `here::here()` resolves correctly regardless of the
working directory). Each script sources `scripts/00_setup.R`.

```
# from the package root:
Rscript environment/validate_versions.R
Rscript scripts/01_mouse_HSPC_processing.R
Rscript scripts/02_mouse_annotation_crossvalidation.R   # optional
Rscript scripts/03_mouse_HSC_subclustering.R
Rscript scripts/04_mouse_Tcell_analysis.R
Rscript scripts/05_human_data_integration.R
Rscript scripts/06_human_HSC_survival_analysis.R
Rscript scripts/07_spatial_transcriptomics_analysis.R
python scripts/08_prepare_cellphonedb_inputs.py
python scripts/08_run_cellphonedb.py --database /path/to/versioned/cellphonedb.zip
Rscript scripts/08_cell_cell_communication.R
Rscript scripts/09_NB4_bulk_RNAseq.R              # independent bulk analysis
```

Before an analysis run, execute the release-gate checks from the repository
root:

```
python code_publication/tests/static_validation.py
python code_publication/tests/test_prepare_cellphonedb_inputs.py
```

### Required input data

1. Mouse HSPC 10x output (this study; Kras and WT) -> data/raw/mouse_HSPC/
2. Mouse T-cell 10x output (this study; CD45+CD3+ sorted, 3 mice pooled per
   genotype) -> data/raw/mouse_Tcell/
3. Human JMML and normal pediatric BM 10x output (JMML: GSE111895;
   normal pediatric BM: GSE155259; local samples were named JMMLID5 and PBM2
   in the original analysis) ->
   data/raw/human_JMML/ and data/raw/human_PB/
4. Healthy mouse BM reference GSE122465 (notlabel.RDS + metaInfo.txt) ->
   data/reference/GSE122465/
5. Spatial transcriptomic data (BMKMANU S1000, Biomarker Technologies;
   GSE313878) -> data/spatial/ST_WT/, data/spatial/ST_Kras/. Script 07 builds
   its RCTD reference directly from script 01's
   `output/mouse_HSPC_annotated.rds`.
6. Bulk RNA-seq + clinical data GSE71449 (ids_exprs.csv, Table_S1.xlsx) ->
   data/bulk/GSE71449/
7. NB4 processed count matrix GSE313879 plus the included GEO-verified
   `sample_metadata.tsv` -> data/bulk/GSE313879/

See data/README.md for details, accessions and the BMK output format.

### Cell-cell communication (CellPhoneDB)

CellPhoneDB uses human ligand-receptor identifiers. The release includes
`data/reference/mouse_to_human_orthologues.tsv`, generated from the official
MGI mouse-human homology report downloaded on 2026-09-30, together with source
and output checksums. Script 08_prepare_cellphonedb_inputs.py excludes ambiguous
one-mouse-to-many-human mappings, sums many-mouse-to-one-human counts, and
writes sparse h5ad plus metadata. The versioned runner then supplies metadata
first and counts second through the official CellPhoneDB v5 Python API:

```
python scripts/08_run_cellphonedb.py \
  --database /path/to/versioned/cellphonedb.zip \
  --iterations 1000 --threshold 0.1 --p-cutoff 0.05 --seed 220625
```

Place the resulting count_network.txt, pvalues.txt, means.txt and
significant_means.txt files into output/cellphonedb/ST_WT/ and
output/cellphonedb/ST_Kras/, then run script 08 to reproduce the network
circle plots and the HSC-centred dot plots (significance threshold P < 0.05).
The runner also writes `cellphonedb_run_metadata.json`, containing package
versions and SHA-256 checksums of the database and analysis inputs.

## Key parameters

| Step | Parameter | Value |
|---|---|---|
| Object creation | min.cells / min.features | 3 / 200 |
| QC | nFeature_RNA range / percent.mt | 500-6,000 / < 10% |
| Variable features | method / n | vst / 2,000 (HSC: 4,000) |
| PCA | n components | 50 |
| Integration | method / dims | Harmony / 1:30 |
| Clustering (all cells) | algorithm / resolution | Louvain (1) / 1.2 |
| Clustering (HSC) | resolution / dims | 0.8 / 1:50 (mouse PCA), 1:30 (human) |
| Marker detection | log2FC / min.pct | 0.5 / 0.1 (mouse); 0.25 / 0.1 (human HSC) |
| Treg DE | log2FC / FDR | 0.25 / < 0.05 |
| GSVA | algorithm / input | `gsvaParam()` / average log-normalized expression |
| Survival signature | algorithm | `ssgseaParam()` followed by `gsva()` |
| Trajectory | tool / start cluster | Slingshot / "1" |
| RCTD | CELL_MIN_INSTANCE / cores / doublet | 20 / 8 / "doublet" |
| HSC niche | neighbourhood radius | 100 coordinate units |
| CellPhoneDB | significance threshold | P < 0.05 |
| Survival | cutoffs | median + survminer::surv_cutpoint |

## Cell-type annotation

Two annotation approaches were used in the original analysis:
1. marker-based manual cluster assignment (final annotation used in the
   manuscript; scripts 01 and 05);
2. label transfer from the GSE122465 healthy bone-marrow reference
   (cross-validation; script 02).

The manuscript reports the marker-based annotation. Mouse annotation is
applied to the active `seurat_clusters` generated at resolution 1.2. The full
human 0-19 cluster map was restored from the retained manuscript-generation
script in initial Git commit `01da390`; script 05 stops if an observed cluster
falls outside that map.

## Notes and uncertainties

- The mouse T-cell data were generated from three mice pooled per genotype
  (one 10x library per genotype); the Kras-vs-WT comparison is therefore based
  on one library per genotype and should be treated as exploratory.
- Human HSC and mouse Treg pathway comparisons are descriptive because the
  retained datasets do not provide independent library-level replication for
  these cell-level contrasts.
- Spatial data are BMKMANU S1000 output and are publicly deposited under
  GSE313878; the experimental protocol is described
  in the Biomarker Technologies methods document (BMKMANU S1000 Spatial
  transcriptomics Materials and method) and summarised in data/README.md.
- The corrected spatial rerun uses the included versioned MGI mouse-human map.
  Script 08 records its SHA-256 checksum; no undocumented online lookup is
  performed. The public GSE313878 record confirms mm10 but omits the exact
  BSTMatrix annotation release and the coordinate/image inputs required by
  script 07; see `data/GEO_PROVENANCE_AUDIT.md`.
- The maintained HSC pseudotime implementation is Slingshot. Historical
  Monocle 2/3 exploratory scripts are not part of this release workflow.

## Reproducibility record

Run `Rscript environment/capture_session_info.R` after installing the required
packages. The DOI-linked `v1.0.0` remains immutable; these corrections are
prepared for a later `v1.1.0` release after numerical validation.

For locally supplied BMKMANU inputs, run
`python tests/validate_spatial_inputs.py`; the author-supplied file checksums
used for the v1.1.0 rerun are recorded in
`data/spatial/INPUT_PROVENANCE.md`.

## Outputs

See output/README.md for the expected outputs per script.

## License / data availability

Analysis code is shared for reproducibility. Raw sequencing data of this study
are deposited under the GEO accessions listed in data/README.md. Supporting
materials not included in GEO are available from the corresponding authors on
reasonable request.
