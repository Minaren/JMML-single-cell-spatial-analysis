# GEO provenance audit (2026-09-30)

This record compares the public GEO metadata with the maintained analysis
package. It is a release gate because repository documentation must not silently
contradict the deposited datasets.

## GSE313553: mouse c-Kit+ HSPC scRNA-seq

Confirmed from the public family SOFT record:

- two libraries: KrasG12D/+ (`GSM9370090`) and WT (`GSM9370091`);
- HSPCs from 6-week-old mice;
- 10x Genomics Chromium Single Cell 3' v3;
- mouse mm10/GRCm38 reference.

Unresolved conflict:

- GEO states Cell Ranger v2.1.1, whereas the author confirmed Cell Ranger
  v7.0.0 on 2026-10-01. The repository records v7.0.0, but the GEO record must
  be corrected or the actual pipeline log must establish which version was
  used.

Additional unresolved demultiplexing record:

- the author supplied `bcl2fastq 5.0.1`; no such release appears in Illumina's
  official bcl2fastq2 version series. The original run log is required before
  this version is reported in the manuscript or release metadata.

## GSE313878: BMKMANU spatial transcriptomics

Confirmed from the public family SOFT record:

- one KrasG12D/+ mouse (`GSM9377043`) and one WT mouse (`GSM9377044`);
- OCT-embedded frozen 10-um sections;
- BMKMANU S1000 is stated for the Kras sample;
- Illumina NovaSeq 6000, PE150;
- assembly recorded as mouse mm10;
- public status since 2026-03-31.

Issues requiring correction or supplementation:

- the WT sample description says `BMKMANU S1001`, inconsistent with the series
  description, the Kras sample and the manufacturer protocol; this appears to
  be a GEO metadata typo and should be corrected to `BMKMANU S1000`;
- the GEO supplement contains only barcodes, features and count matrices. The
  author supplied the coordinate files and H&E images on 2026-10-01; their
  dimensions and checksums are recorded in `spatial/INPUT_PROVENANCE.md`.
  These supporting files should still be added to GEO. Script 07 now derives
  its RCTD reference from the versioned script-01 output rather than separate
  unversioned reference exports;
- GEO records BSTMatrix and mm10 but not a BSTMatrix version or annotation
  release. The manuscript records BSTMatrix v1.0; retain the vendor run report
  or pipeline log as supporting provenance.

## GSE313879: NB4 bulk RNA-seq

Confirmed from the public family SOFT record and processed count matrix:

| Matrix column | Condition | GEO sample |
|---|---|---|
| NB4_NC1 | NC | GSM9377045 |
| NB4_NC2 | NC | GSM9377046 |
| NB4_NC3 | NC | GSM9377047 |
| NB4_OE1 | OE | GSM9377048 |
| NB4_OE2 | OE | GSM9377049 |
| NB4_OE3 | OE | GSM9377050 |

The verified mapping is included as
`data/bulk/GSE313879/sample_metadata.tsv`. The downloaded processed matrix had
SHA-256 `8da6edc29b1c579b3306fc581883ad10b88924f2b51e12b1e801da33e89f2a65`.

## Source records

- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313553>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313878>
- <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313879>
- <https://www.10xgenomics.com/support/software/cell-ranger/latest/release-notes/cr-release-notes>
- <https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html>
