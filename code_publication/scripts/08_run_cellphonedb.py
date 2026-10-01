#!/usr/bin/env python3
"""Run CellPhoneDB v5 statistical analysis with versioned inputs.

The official Python API is used so that the metadata/counts order and database
path are explicit. Standardised output filenames are written for the R plotting
script, together with software versions, parameters and SHA-256 checksums.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def write_standard_outputs(results: dict, output: Path, p_cutoff: float) -> None:
    key_map = {
        "pvalues": "pvalues.txt",
        "means": "means.txt",
        "significant_means": "significant_means.txt",
        "deconvoluted": "deconvoluted.txt",
    }
    for key, filename in key_map.items():
        value = results.get(key)
        if isinstance(value, pd.DataFrame):
            value.to_csv(output / filename, sep="\t", index=False)

    pvalues = results.get("pvalues")
    if not isinstance(pvalues, pd.DataFrame):
        raise RuntimeError("CellPhoneDB did not return a pvalues DataFrame")
    pair_columns = [column for column in pvalues.columns if "|" in str(column)]
    if not pair_columns:
        raise RuntimeError("No cell-pair columns were found in CellPhoneDB pvalues")

    network_rows = []
    for pair in pair_columns:
        source, target = str(pair).split("|", maxsplit=1)
        values = pd.to_numeric(pvalues[pair], errors="coerce")
        network_rows.append(
            {"SOURCE": source, "TARGET": target, "count": int((values < p_cutoff).sum())}
        )
    pd.DataFrame(network_rows).to_csv(output / "count_network.txt", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--database", type=Path, required=True,
                        help="Versioned CellPhoneDB database ZIP")
    parser.add_argument("--samples", nargs="+", default=["ST_WT", "ST_Kras"])
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--seed", type=int, default=220625)
    args = parser.parse_args()

    project = args.project.resolve()
    database = args.database.resolve()
    if not database.exists():
        raise FileNotFoundError(database)

    versions = {
        name: package_version(name)
        for name in ("cellphonedb", "anndata", "numpy", "pandas", "scipy")
    }
    if versions["cellphonedb"] != "5.0.1":
        raise RuntimeError(
            "This release requires CellPhoneDB 5.0.1; installed version is "
            + versions["cellphonedb"]
        )

    for sample in args.samples:
        output = project / "output" / "cellphonedb" / sample
        output.mkdir(parents=True, exist_ok=True)
        meta = output / f"meta_{sample}.tsv"
        counts = output / f"counts_{sample}.h5ad"
        missing = [str(path) for path in (meta, counts) if not path.exists()]
        if missing:
            raise FileNotFoundError("Missing prepared input(s):\n" + "\n".join(missing))

        results = cpdb_statistical_analysis_method.call(
            cpdb_file_path=str(database),
            meta_file_path=str(meta),
            counts_file_path=str(counts),
            counts_data="hgnc_symbol",
            output_path=str(output),
            iterations=args.iterations,
            threshold=args.threshold,
            pvalue=args.p_cutoff,
            threads=args.threads,
            debug_seed=args.seed,
        )
        write_standard_outputs(results, output, args.p_cutoff)

        record = {
            "sample": sample,
            "run_utc": datetime.now(timezone.utc).isoformat(),
            "database_file": database.name,
            "database_sha256": sha256(database),
            "meta_sha256": sha256(meta),
            "counts_sha256": sha256(counts),
            "parameters": {
                "iterations": args.iterations,
                "threshold": args.threshold,
                "p_cutoff": args.p_cutoff,
                "threads": args.threads,
                "debug_seed": args.seed,
                "counts_data": "hgnc_symbol",
            },
            "python_packages": versions,
        }
        (output / "cellphonedb_run_metadata.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
