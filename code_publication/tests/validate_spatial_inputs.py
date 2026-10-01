#!/usr/bin/env python3
"""Validate local BMKMANU S1000 matrices, coordinates and H&E inputs.

The raw files are intentionally excluded from Git. Run this validator after
placing them under ``data/spatial/ST_WT`` and ``data/spatial/ST_Kras``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gzip_rows(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            yield line_no, line.rstrip("\r\n").split("\t")


def matrix_market_shape(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"{path}: not a coordinate Matrix Market file")
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"{path}: invalid dimensions line")
            return tuple(map(int, fields))
    raise ValueError(f"{path}: dimensions line is absent")


def validate_sample(sample_dir: Path) -> dict[str, object]:
    required = {
        "matrix": sample_dir / "matrix.mtx.gz",
        "features": sample_dir / "features.tsv.gz",
        "barcodes": sample_dir / "barcodes.tsv.gz",
        "coordinates": sample_dir / "barcodes_pos.tsv.gz",
        "reads": sample_dir / "barcodes_read.tsv.gz",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing spatial input(s):\n" + "\n".join(missing))

    barcodes = []
    for line_no, row in read_gzip_rows(required["barcodes"]):
        if len(row) != 1 or not row[0]:
            raise ValueError(f"{required['barcodes']}:{line_no}: invalid barcode row")
        barcodes.append(row[0])
    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"{required['barcodes']}: duplicate barcodes")

    feature_count = 0
    for line_no, row in read_gzip_rows(required["features"]):
        if len(row) < 2 or not row[1]:
            raise ValueError(f"{required['features']}:{line_no}: invalid feature row")
        feature_count += 1

    coordinate_barcodes = []
    for line_no, row in read_gzip_rows(required["coordinates"]):
        if len(row) != 3:
            raise ValueError(f"{required['coordinates']}:{line_no}: expected 3 columns")
        float(row[1]); float(row[2])
        coordinate_barcodes.append(row[0])

    read_barcodes = []
    total_reads = 0.0
    for line_no, row in read_gzip_rows(required["reads"]):
        if len(row) != 2:
            raise ValueError(f"{required['reads']}:{line_no}: expected 2 columns")
        value = float(row[1])
        if value < 0:
            raise ValueError(f"{required['reads']}:{line_no}: negative read count")
        read_barcodes.append(row[0])
        total_reads += value

    if barcodes != coordinate_barcodes or barcodes != read_barcodes:
        raise ValueError(f"{sample_dir}: barcode, coordinate and read rows are not aligned")
    n_features, n_barcodes, nnz = matrix_market_shape(required["matrix"])
    if (n_features, n_barcodes) != (feature_count, len(barcodes)):
        raise ValueError(
            f"{sample_dir}: matrix shape {(n_features, n_barcodes)} does not match "
            f"{feature_count} features x {len(barcodes)} barcodes"
        )

    image_files = sorted(sample_dir.glob("*.tif"))
    if len(image_files) != 1:
        raise ValueError(f"{sample_dir}: expected exactly one H&E TIFF image")
    with image_files[0].open("rb") as handle:
        if handle.read(4) not in (b"II*\x00", b"MM\x00*"):
            raise ValueError(f"{image_files[0]}: invalid TIFF signature")

    return {
        "sample": sample_dir.name,
        "features": n_features,
        "barcodes": n_barcodes,
        "nonzero_entries": nnz,
        "total_spot_reads": total_reads,
        "matrix_sha256": sha256(required["matrix"]),
        "coordinates_sha256": sha256(required["coordinates"]),
        "image": image_files[0].name,
        "image_sha256": sha256(image_files[0]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="code_publication project root",
    )
    args = parser.parse_args()
    summaries = [
        validate_sample(args.project / "data" / "spatial" / sample)
        for sample in ("ST_WT", "ST_Kras")
    ]
    print(json.dumps(summaries, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
