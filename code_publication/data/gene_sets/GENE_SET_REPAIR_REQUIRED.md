# Gene-set repair record

The v1.1.0 audit repaired spreadsheet date conversion and capital-I/lowercase-l
transcription errors in the retained mouse gene-set files. The repaired symbols
use the historical mouse symbols expected by the mm10-era expression matrices;
some now have newer MGI symbols (for example, `Sept1` is now `Septin1`, and
`Sep15` is now `Selenof`).

| Corrupted value | Restored value |
|---|---|
| `1-Sep`, `2-Sep`, `6-Sep`, `7-Sep`, `11-Sep` | `Sept1`, `Sept2`, `Sept6`, `Sept7`, `Sept11` |
| `15-Sep` | `Sep15` |
| `6-Mar`, `7-Mar` | `March6`, `March7` |
| leading lowercase `l` in known symbols | leading uppercase `I` (for example `Ifitm1`, `Irf8`, `Itga2b`) |
| embedded capital `I` in known lowercase-l positions | lowercase `l` (for example `Epb41l4b`, `Ifi27l2a`, `Pcp4l1`, `Med13l`, `mt-Nd4l`) |

The original unmodified values remain recoverable from initial Git commit
`01da390`. NCBI Gene/MGI records were used to check symbol identity and aliases:

- <https://www.ncbi.nlm.nih.gov/gene/54204> (`Septin1`; alias `Sept1`)
- <https://www.ncbi.nlm.nih.gov/gene/93684> (`Selenof`; alias `Sep15`)
- <https://www.ncbi.nlm.nih.gov/gene/223455> (`Marchf6`; alias `March6`)
- <https://www.ncbi.nlm.nih.gov/gene/57438> (`Marchf7`; alias `March7`)

`read_gene_sets()` still stops if a date-like value reappears. Duplicate
set-gene pairs are de-duplicated in memory by `read_gene_sets()` and are left in
the source file so the retained list remains auditable against its origin.
