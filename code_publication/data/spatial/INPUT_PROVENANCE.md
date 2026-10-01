# Local spatial-input provenance (2026-10-01)

The raw BMKMANU S1000 files and H&E images are not committed to Git because of
their size. They were supplied by the author for the v1.1.0 numerical rerun and
were validated with `tests/validate_spatial_inputs.py`.

| Sample | Matrix dimensions | Nonzero entries | Matrix SHA-256 | Coordinate SHA-256 | H&E image | Image dimensions | Image SHA-256 |
|---|---:|---:|---|---|---|---:|---|
| ST_WT | 54,838 × 64,240 | 9,359,852 | `30367794928bcd371a72c5f68815866626e93383cdfd6162361346b9160e577d` | `4bc9a4b3e4adf0afa3c451d09ff4105afcb83bf10ea7796f05ea45752afa66c9` | `he-WT.tif` | 5,287 × 5,337 | `126d0f48147df9a006cdf39fac35696b6632d716f46a34f3191114cfbd767614` |
| ST_Kras | 54,838 × 64,240 | 22,834,202 | `1ddb4cb3ee36b69495d9f8dcd09353d8a93dc039cf536ea45b5a74dd68c575d4` | `199692331dfccbacf40df27e74876ecce48a60d34e9c4f0857a370e329cec57b` | `he-Kras.tif` | 10,544 × 10,665 | `d81d7d57db6be96f22116238d6d8936584469ac5a6af0ed494605e8a112178ee` |

For both samples, the matrix dimensions agree with 54,838 feature rows and
64,240 barcode rows. Barcode order is identical across `barcodes.tsv.gz`,
`barcodes_pos.tsv.gz` and `barcodes_read.tsv.gz`. WT and Kras use the same
64,240-position coordinate grid.

The uploaded archive containing the count/coordinate files had SHA-256
`7c27eb8c9af43de4c5ef5115394adadc23986f9477d4275c54e36510e56b8fa0`.
