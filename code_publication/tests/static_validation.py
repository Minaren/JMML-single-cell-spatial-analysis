#!/usr/bin/env python3
"""Dependency-free release-gate checks for the publication package."""

from __future__ import annotations

import csv
import hashlib
import json
import py_compile
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "code_publication"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def validate_gene_sets() -> None:
    date_pattern = re.compile(
        r"^[0-9]{1,2}-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$",
        re.IGNORECASE,
    )
    known_corruptions = {
        "l7Rn6", "lgf2bp2", "lrf8", "lfitm6", "lgf1", "lrak2",
        "lfitm7", "lfitm3", "lfitm1", "Epb41I4b", "mt-Nd4I",
        "lrf2bp2", "ltgb5", "lghm", "lkzf2", "lfi203", "lkzf1",
        "lfi27I2a", "lqgap1", "ll17ra", "ll12a", "Pcp4I1", "lfi27",
        "ltga2b", "ll7r", "Med13I",
    }
    for path in sorted((PACKAGE / "data" / "gene_sets").glob("*.csv")):
        with path.open(newline="", encoding="utf-8") as handle:
            for line_no, row in enumerate(csv.reader(handle), 1):
                require(len(row) == 2, f"{path}:{line_no}: expected exactly two columns")
                gene = row[1].strip()
                require(not date_pattern.fullmatch(gene),
                        f"{path}:{line_no}: date-converted gene symbol {gene}")
                require(gene not in known_corruptions,
                        f"{path}:{line_no}: known I/l-corrupted gene symbol {gene}")


def validate_r_sources() -> None:
    setup = read(PACKAGE / "scripts" / "00_setup.R")
    survival = read(PACKAGE / "scripts" / "06_human_HSC_survival_analysis.R")
    human = read(PACKAGE / "scripts" / "05_human_data_integration.R")
    mouse_hsc = read(PACKAGE / "scripts" / "03_mouse_HSC_subclustering.R")
    spatial = read(PACKAGE / "scripts" / "07_spatial_transcriptomics_analysis.R")

    require("GSVA::ssgseaParam" in setup and "GSVA::gsva(param" in setup,
            "The GSVA 2.2.0 ssGSEA parameter-object API is not present")
    require("Table_S1_nonempty_Cause_of_death_is_event" in survival,
            "The prespecified GSE71449 event rule is missing")
    require("status_candidates" not in survival,
            "Speculative survival-status-column selection remains")

    for cluster_id in range(20):
        require(re.search(rf"(?:==|%in% c\([^)]*)\s*{cluster_id}(?:\s|,|\))", human),
                f"Human cluster {cluster_id} is absent from the restored map")
    require("Unassigned_" in human and "stop(" in human,
            "Human mapping does not retain its coverage assertion")

    for tf in ("Nfkb1", "Fos", "Stat3", "Jun"):
        require(re.search(rf"\b{tf}\s*=\s*c\(", mouse_hsc),
                f"Restored {tf} target module is missing")
    require("tf_targets_mouse.csv" not in mouse_hsc,
            "Script 03 still depends on the incorrectly reported missing TF file")
    require("mouse_HSPC_annotated.rds" in spatial,
            "Script 07 does not use the versioned script-01 RCTD reference")
    require('file.path(spatial_dir, "sc_meta.txt")' not in spatial and
            'file.path(spatial_dir, "ref_cell_anno")' not in spatial,
            "Script 07 still depends on unversioned dense RCTD exports")

    for path in sorted((PACKAGE / "scripts").glob("*.R")):
        source = read(path)
        executable_source = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("#")
        )
        require(not re.search(r"(?<![A-Za-z])[A-Za-z]:[/\\\\]", executable_source),
                f"Absolute Windows path remains in {path}")
        require("JoinLayers(" not in executable_source and
                "layer=" not in executable_source.replace(" ", ""),
                f"Seurat v5-only API remains in {path}")


def validate_python_sources() -> None:
    python_paths = list((PACKAGE / "scripts").glob("*.py"))
    python_paths += list((PACKAGE / "tests").glob("*.py"))
    for path in sorted(python_paths):
        py_compile.compile(str(path), doraise=True)
    requirements = read(PACKAGE / "environment" / "requirements-cellphonedb.txt")
    lock = read(PACKAGE / "environment" / "requirements-cellphonedb-lock.txt")
    runner = read(PACKAGE / "scripts" / "08_run_cellphonedb.py")
    for pin in ("cellphonedb==5.0.1", "anndata==0.10.9", "numpy==1.26.4",
                "pandas==2.2.3", "scipy==1.13.1"):
        require(pin in requirements and pin in lock,
                f"Required Python dependency is not pinned consistently: {pin}")
    require('versions["cellphonedb"] != "5.0.1"' in runner,
            "The CellPhoneDB runner does not enforce 5.0.1")
    versions = read(PACKAGE / "environment" / "software_versions.tsv")
    require("Cell Ranger\t7.0.0" in versions,
            "Cell Ranger is not recorded as the author-confirmed v7.0.0")
    require("bcl2fastq2\t2.20.0" in versions,
            "bcl2fastq2 is not recorded as the author-confirmed v2.20.0")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_versioned_inputs() -> None:
    ref_dir = PACKAGE / "data" / "reference"
    mapping = ref_dir / "mouse_to_human_orthologues.tsv"
    metadata_path = ref_dir / "mouse_to_human_orthologues.metadata.json"
    metadata = json.loads(read(metadata_path))
    require(file_sha256(mapping) == metadata["output_sha256"],
            "Mouse-human mapping checksum does not match its metadata")
    with mapping.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    require(rows[0] == ["mouse_symbol", "human_symbol"],
            "Mouse-human mapping header is invalid")
    require(len(rows) - 1 == metadata["pair_count"] == 24_584,
            "Mouse-human mapping pair count is invalid")

    sample_metadata = PACKAGE / "data" / "bulk" / "GSE313879" / "sample_metadata.tsv"
    with sample_metadata.open(newline="", encoding="utf-8") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    require([row["sample_id"] for row in samples] ==
            ["NB4_NC1", "NB4_NC2", "NB4_NC3", "NB4_OE1", "NB4_OE2", "NB4_OE3"],
            "GSE313879 sample IDs are not in the GEO-verified order")
    require([row["condition"] for row in samples] == ["NC"] * 3 + ["OE"] * 3,
            "GSE313879 conditions are invalid")


def main() -> None:
    validate_gene_sets()
    validate_r_sources()
    validate_python_sources()
    validate_versioned_inputs()
    print("Static release-gate checks passed.")


if __name__ == "__main__":
    main()
