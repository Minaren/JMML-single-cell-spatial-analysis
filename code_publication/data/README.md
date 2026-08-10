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
│   └── GSE122465/
│       ├── notlabel.RDS          # healthy mouse BM reference (Seurat object)
│       └── metaInfo.txt          # cell type labels (scmap-based and paper-based)
├── gene_sets/
│   ├── gsva_mouse_cluster.csv    # curated module gene sets (mouse)
│   ├── gsva_human_cluster.csv    # curated module gene sets (human)
│   └── function_mouse_HSC.csv    # HSC functional gene sets (mouse)
├── spatial/
│   ├── ST_WT/  ST_Kras/          # BMKMANU S1000 output (see format note below)
│   ├── sc_meta.txt               # single-cell reference counts matrix (genes x cells)
│   └── ref_cell_anno             # barcode <TAB> cell_type reference annotation
└── bulk/
    └── GSE71449/
        ├── ids_exprs.csv         # bulk expression matrix
        └── Table_S1.xlsx         # clinical data
```

## Public data accessions

| Dataset | Accession | Contents |
|---|---|---|
| Human JMML scRNA-seq | GSE155295 | JMML bone-marrow mononuclear cells |
| Normal pediatric BM scRNA-seq | GSE111895 | normal pediatric bone marrow |
| JMML bulk RNA-seq + clinical | GSE71449 | bulk expression + Table_S1 clinical data |
| Healthy mouse BM reference | GSE122465 | bone marrow cell types incl. stroma (Baryawno et al., Cell 2019) |

NOTE: the manuscript Figure 8 legend lists the JMML scRNA-seq accession as
GSE155295; the earlier draft methods text contained the typo "GSE155259".
Please confirm the correct accessions against GEO before submission. The
original analysis code read local 10x directories (JMMLID5, PBM2); confirm
the correspondence between these local samples and the accessions above.

## Spatial transcriptomic data (not public)

Spatial transcriptomics was performed by Biomarker Technologies Corporation
(Beijing, China) with the BMKMANU S1000 platform (fresh-frozen sections,
10 um, H&E staining, 20x bright-field imaging, Illumina NovaSeq 6000 PE150,
>=50,000 reads per spot). The raw data (fastq, HE images, count matrices) are
not publicly deposited; access is available from the corresponding authors or
Biomarker Technologies.

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