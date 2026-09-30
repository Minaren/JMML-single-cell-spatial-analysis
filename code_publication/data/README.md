# data/

Contains all auxiliary input data required by the analysis scripts. Raw
sequencing data (fastq, 10x matrices, spatial transcriptomic images) are NOT
included in this repository; see README.md and the notes below for how to
obtain and where to place them.

## Directory layout (as expected by the scripts)

```
data/
├── raw/
│   ├── mouse_HSPC/Kras/          # 10x Cell Ranger output (barcodes/features/matrix)
│   ├── mouse_HSPC/WT/            # 10x Cell Ranger output
│   ├── mouse_Tcell/Kras/         # CD45+CD3+ sorted T cells (3 mice pooled per genotype)
│   ├── mouse_Tcell/WT/
│   ├── human_JMML/               # JMML bone marrow 10x output (local sample "JMMLID5")
│   └── human_PB/                 # normal pediatric bone-marrow 10x output (local sample "PBM2")
├── reference/
│   ├── GSE122465/
│       ├── notlabel.RDS          # healthy mouse BM reference (Seurat object)
│       └── metaInfo.txt          # cell type labels (scmap-based and paper-based)
│   └── mouse_to_human_orthologues.tsv  # explicit mapping used for CellPhoneDB
├── gene_sets/
│   ├── gsva_mouse_cluster.csv    # curated module gene sets (mouse)
│   ├── gsva_human_cluster.csv    # curated module gene sets (human)
│   ├── function_mouse_HSC.csv    # HSC functional gene sets (mouse)
│   └── tf_targets_mouse.csv      # exact Figure 2F targets (set,gene; author supplied)
├── spatial/
│   ├── ST_WT/  ST_Kras/          # BMKMANU S1000 output (see format note below)
│   ├── sc_meta.txt               # single-cell reference counts matrix (genes x cells)
│   └── ref_cell_anno             # barcode <TAB> cell_type reference annotation
└── bulk/
    ├── GSE71449/
    │   ├── ids_exprs.csv         # bulk expression matrix
    │   └── Table_S1.xlsx         # clinical data
    └── GSE313879/
        ├── gene_count_matrix.tsv.gz
        └── sample_metadata.tsv   # sample_id <TAB> condition (NC or OE)
```

## Public data accessions

| Dataset | Accession | Contents |
|---|---|---|
| Human JMML scRNA-seq | GSE111895 | JMML bone-marrow samples, including JMML_ID5 |
| Normal pediatric BM scRNA-seq | GSE155259 | developmental hematopoiesis dataset, including pediatric BM |
| JMML bulk RNA-seq + clinical | GSE71449 | bulk expression + Table_S1 clinical data |
| Healthy mouse BM reference | GSE122465 | bone marrow cell types incl. stroma (Baryawno et al., Cell 2019) |
| Mouse c-Kit+ HSPC scRNA-seq generated in this study | GSE313553 | WT and KrasG12D/+ HSPCs from 6-week-old mice |
| Mouse bone-marrow spatial transcriptomics generated in this study | GSE313878 | WT and KrasG12D/+ spatial matrices |
| NB4 CD69-overexpression bulk RNA-seq generated in this study | GSE313879 | three control and three CD69-OE samples |

The JMML/normal mapping and 6-week mouse age were checked against public GEO
records and confirmed by the author. The original analysis code read local 10x
directories JMMLID5 and PBM2.

For GSE313879, `sample_metadata.tsv` must contain exactly three `NC` and three
`OE` rows, and `sample_id` must exactly match the processed count-matrix column
names. The script never infers biological groups from filenames.

The retained original repository did not contain the exact transcription-factor
target lists used for Figure 2F. Add them as a two-column, header-free
`tf_targets_mouse.csv` (`set,gene`). Script 03 skips that panel with a warning
when the file is absent rather than using illustrative targets.

## Spatial transcriptomic data

Spatial transcriptomics was performed by Biomarker Technologies Corporation
(Beijing, China) with the BMKMANU S1000 platform (fresh-frozen sections,
10 um, H&E staining, 20x bright-field imaging, Illumina NovaSeq 6000 PE150,
approximately 100,000 paired-end reads per 100 um spatial analysis unit and
approximately 60-150 Gb per sample, according to the manufacturer's protocol).
The study's spatial data are publicly available under
GSE313878. Place the downloaded files into the directory layout below.

BMKMANU S1000 output format (as received, one folder per sample):
- matrix.mtx.gz        # sparse count matrix (genes x spots)
- features.tsv.gz      # 3 columns: Ensembl ID, gene symbol, feature type
- barcodes.tsv.gz      # spot barcodes (one per line)
- barcodes_pos.tsv.gz  # barcode, x, y spot coordinates
- barcodes_read.tsv.gz # barcode, number of reads
- barcodes_cluster.tsv.gz  # bin-to-spot cluster mapping (informational)
- he-<sample>.tif      # H&E bright-field image

The scripts expect the folders to be named `ST_WT` and `ST_Kras` under
data/spatial/ (if your folders differ, adjust the paths in
07_spatial_transcriptomics_analysis.R).

## Reference files used by RCTD (data/spatial/)

- sc_meta.txt  : gene x cell counts matrix of the merged mouse bone-marrow
                 single-cell reference (used to build the RCTD Reference object).
- ref_cell_anno: two-column file (barcode <TAB> cell_type) matching sc_meta.txt.
                 Cell types represented by <=25 cells are excluded automatically
                 in the script (as in the original analysis).

## Mouse-to-human orthologue mapping for CellPhoneDB

CellPhoneDB requires human identifiers. Provide
`data/reference/mouse_to_human_orthologues.tsv` with two tab-separated columns:
`mouse_symbol` and `human_symbol`. Archive a companion README containing the
mapping database name, version/release, download date and source URL. The
preparation script excludes ambiguous one-mouse-to-many-human mappings and
sums many-mouse-to-one-human counts. Do not substitute an undocumented online
lookup, because that would make the interaction analysis non-versioned.
