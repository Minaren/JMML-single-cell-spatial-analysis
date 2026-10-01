#!/usr/bin/env python3
"""Prepare sparse human-orthologue CellPhoneDB inputs from BMKMANU S1000 data.

The mapping table is deliberately an explicit input rather than an implicit
online lookup so the orthologue resource, version and one-to-many decisions can
be archived with the release.  Required columns are ``mouse_symbol`` and
``human_symbol``. Ambiguous one-mouse-to-many-human mappings are excluded;
many-mouse-to-one-human counts are summed.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


def read_mtx(path: Path) -> sparse.csr_matrix:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        return sparse.csr_matrix(mmread(handle))


def prepare_sample(project: Path, sample: str, mapping_path: Path) -> None:
    source = project / "data" / "spatial" / sample
    annotation_path = project / "output" / "spatial" / sample / "Spatial_CellType.tsv"
    output = project / "output" / "cellphonedb" / sample
    output.mkdir(parents=True, exist_ok=True)

    required = [
        source / "matrix.mtx.gz",
        source / "features.tsv.gz",
        source / "barcodes.tsv.gz",
        annotation_path,
        mapping_path,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(missing))

    matrix = read_mtx(source / "matrix.mtx.gz")
    features = pd.read_csv(source / "features.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(source / "barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    annotations = pd.read_csv(annotation_path, sep="\t", dtype=str)
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)

    if features.shape[1] < 2:
        raise ValueError("features.tsv.gz must contain a gene-symbol column in position 2")
    if not {"barcode", "cell_type"}.issubset(annotations.columns):
        raise ValueError("Spatial_CellType.tsv must contain barcode and cell_type columns")
    if not {"mouse_symbol", "human_symbol"}.issubset(mapping.columns):
        raise ValueError("Orthologue table must contain mouse_symbol and human_symbol columns")

    genes = features.iloc[:, 1].astype(str)
    if matrix.shape == (len(barcodes), len(genes)):
        matrix = matrix.T
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match {len(genes)} genes x "
            f"{len(barcodes)} barcodes"
        )

    annotations = annotations.dropna(subset=["barcode", "cell_type"]).drop_duplicates("barcode")
    barcode_to_col = pd.Series(np.arange(len(barcodes)), index=barcodes)
    common = annotations.loc[annotations["barcode"].isin(barcode_to_col.index)].copy()
    if common.empty:
        raise ValueError(f"No annotated barcodes overlap the count matrix for {sample}")
    matrix = matrix[:, barcode_to_col.loc[common["barcode"]].to_numpy()]

    # Aggregate duplicated feature symbols before applying the orthologue map.
    # This prevents ambiguous pandas indexing and preserves the total counts.
    unique_mouse = pd.Index(pd.unique(genes))
    mouse_to_row = pd.Series(np.arange(len(unique_mouse)), index=unique_mouse)
    aggregate_mouse = sparse.csr_matrix(
        (
            np.ones(len(genes), dtype=np.int8),
            (mouse_to_row.loc[genes].to_numpy(), np.arange(len(genes))),
        ),
        shape=(len(unique_mouse), len(genes)),
    )
    matrix = aggregate_mouse @ matrix
    genes = unique_mouse

    mapping = mapping[["mouse_symbol", "human_symbol"]].dropna().copy()
    mapping["mouse_symbol"] = mapping["mouse_symbol"].str.strip()
    mapping["human_symbol"] = mapping["human_symbol"].str.strip()
    mapping = mapping[(mapping["mouse_symbol"] != "") & (mapping["human_symbol"] != "")]
    mapping = mapping.drop_duplicates()

    # Exclude ambiguous one-to-many mappings. Multiple mouse symbols mapping
    # to one human symbol are combined by sparse matrix multiplication.
    unambiguous_mouse = mapping.groupby("mouse_symbol")["human_symbol"].nunique()
    unambiguous_mouse = unambiguous_mouse[unambiguous_mouse == 1].index
    mapping = mapping[mapping["mouse_symbol"].isin(unambiguous_mouse)]
    mapping = mapping[mapping["mouse_symbol"].isin(set(genes))]
    if mapping.empty:
        raise ValueError("No count-matrix genes matched the supplied orthologue table")

    gene_to_row = pd.Series(np.arange(len(genes)), index=genes)
    human_symbols = pd.Index(sorted(mapping["human_symbol"].unique()))
    human_to_row = pd.Series(np.arange(len(human_symbols)), index=human_symbols)
    map_matrix = sparse.csr_matrix(
        (
            np.ones(len(mapping), dtype=np.int8),
            (
                human_to_row.loc[mapping["human_symbol"]].to_numpy(),
                gene_to_row.loc[mapping["mouse_symbol"]].to_numpy(),
            ),
        ),
        shape=(len(human_symbols), len(genes)),
    )
    human_counts = map_matrix @ matrix

    obs = common[["barcode", "cell_type"]].copy().set_index("barcode")
    # Force Python-object strings rather than pandas' nullable StringArray.
    # anndata 0.10.x cannot serialise StringArray-backed indices/columns when
    # used with newer pandas releases.
    obs.index = pd.Index(
        [str(value).replace("-", "_") for value in obs.index],
        dtype=object,
        name="Cell",
    )
    obs["cell_type"] = obs["cell_type"].map(str).astype(object)
    if obs.index.has_duplicates:
        raise ValueError("Barcode sanitisation produced duplicate CellPhoneDB cell names")
    var = pd.DataFrame(index=pd.Index([str(value) for value in human_symbols],
                                     dtype=object, name="gene"))
    adata = ad.AnnData(X=human_counts.T.tocsr(), obs=obs, var=var)
    adata.write_h5ad(output / f"counts_{sample}.h5ad", compression="gzip")

    meta = obs.reset_index()
    meta.to_csv(output / f"meta_{sample}.tsv", sep="\t", index=False)
    mapping_sha256 = hashlib.sha256(mapping_path.read_bytes()).hexdigest()
    pd.DataFrame(
        {
            "metric": [
                "mapping_file",
                "mapping_sha256",
                "input_mouse_genes",
                "mapped_mouse_genes",
                "output_human_genes",
                "annotated_spots",
            ],
            "value": [
                str(mapping_path),
                mapping_sha256,
                str(len(genes)),
                str(mapping["mouse_symbol"].nunique()),
                str(len(human_symbols)),
                str(len(obs)),
            ],
        }
    ).to_csv(output / f"orthologue_mapping_report_{sample}.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--mapping",
        type=Path,
        default=None,
        help="TSV with mouse_symbol and human_symbol columns",
    )
    parser.add_argument("--samples", nargs="+", default=["ST_WT", "ST_Kras"])
    args = parser.parse_args()
    mapping = args.mapping or args.project / "data" / "reference" / "mouse_to_human_orthologues.tsv"
    for sample in args.samples:
        prepare_sample(args.project.resolve(), sample, mapping.resolve())


if __name__ == "__main__":
    main()
