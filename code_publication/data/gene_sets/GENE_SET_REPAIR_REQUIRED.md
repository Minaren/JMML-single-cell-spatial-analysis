# Gene-set source repair required

The following entries have the characteristic form of spreadsheet date
conversion and must be restored from the original gene-set source before the
v1.1.0 analysis is run. They are intentionally not guessed in this repository.

| File | Suspicious entries |
|---|---|
| `gsva_mouse_cluster.csv` | `1-Sep` (two rows) |
| `function_mouse_HSC.csv` | `1-Sep`, `7-Sep`, `2-Sep`, `7-Mar`, `15-Sep`, `6-Mar`, `6-Sep`, `11-Sep` |

After restoration, record the source publication/database and perform a
case-sensitive match against the mouse expression matrix. `read_gene_sets()`
stops when any date-like entry remains, preventing silent use of corrupted gene
symbols.

`function_mouse_HSC.csv` also contains duplicated set-gene pairs. These should
be deduplicated only after the original source has been recovered and checked.
