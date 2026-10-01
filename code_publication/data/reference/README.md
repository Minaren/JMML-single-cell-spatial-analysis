# Versioned reference data

## Mouse-human orthologues

`mouse_to_human_orthologues.tsv` was generated from the official Mouse Genome
Informatics `HOM_MouseHumanSequence.rpt` report downloaded on 2026-09-30:

<https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt>

The source report described mouse coordinates on GRCm39 and human coordinates
on GRCh38. The orthologue file contains symbol pairs sharing an MGI homology
class; script `08_prepare_cellphonedb_inputs.py` excludes mouse symbols with
multiple human counterparts before preparing CellPhoneDB inputs. Source and
output SHA-256 checksums are stored in
`mouse_to_human_orthologues.metadata.json`.

Regenerate the file with:

```bash
python scripts/00_build_mgi_orthologues.py /path/to/HOM_MouseHumanSequence.rpt \
  --download-date 2026-09-30
```

This mapping defines the corrected v1.1.0 rerun; it is not represented as the
unrecorded mapping resource used by the historical analysis.
