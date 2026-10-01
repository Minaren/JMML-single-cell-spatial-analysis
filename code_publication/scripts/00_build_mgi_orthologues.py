#!/usr/bin/env python3
"""Build the versioned mouse-human symbol map from the official MGI report.

Download ``HOM_MouseHumanSequence.rpt`` from the MGI reports site, then pass
the local file to this script. Every mouse-human pair in the same MGI homology
class is emitted. The CellPhoneDB input-preparation script subsequently excludes
mouse symbols with more than one human counterpart and records the mapping
file checksum used for each sample.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import date
from pathlib import Path


SOURCE_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path, help="Downloaded HOM_MouseHumanSequence.rpt")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "data" / "reference" / "mouse_to_human_orthologues.tsv",
    )
    parser.add_argument("--download-date", default=date.today().isoformat())
    args = parser.parse_args()

    groups: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"mouse": set(), "human": set()}
    )
    with args.report.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"DB Class Key", "Common Organism Name", "Symbol"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unexpected MGI report columns")
        for row in reader:
            organism = row["Common Organism Name"]
            symbol = row["Symbol"].strip()
            if not symbol:
                continue
            if organism == "mouse, laboratory":
                groups[row["DB Class Key"]]["mouse"].add(symbol)
            elif organism == "human":
                groups[row["DB Class Key"]]["human"].add(symbol)

    pairs = sorted(
        {
            (mouse, human)
            for group in groups.values()
            for mouse in group["mouse"]
            for human in group["human"]
        }
    )
    if len(pairs) < 18_000:
        raise RuntimeError(f"Unexpectedly small MGI mouse-human map: {len(pairs)} pairs")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["mouse_symbol", "human_symbol"])
        writer.writerows(pairs)

    metadata = {
        "source": "Mouse Genome Informatics HOM_MouseHumanSequence report",
        "source_url": SOURCE_URL,
        "download_date": args.download_date,
        "source_sha256": sha256(args.report),
        "output_file": args.output.name,
        "output_sha256": sha256(args.output),
        "pair_count": len(pairs),
        "filtering": "all mouse-human pairs sharing an MGI DB Class Key",
    }
    metadata_path = args.output.with_suffix(".metadata.json")
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
