#!/usr/bin/env python3
"""Synthetic regression test for CellPhoneDB input preparation."""

from __future__ import annotations

import gzip
import importlib.util
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "code_publication" / "scripts" / "08_prepare_cellphonedb_inputs.py"


def load_prepare_sample():
    spec = importlib.util.spec_from_file_location("prepare_cpdb", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.prepare_sample


def write_gzip_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)


def main() -> None:
    prepare_sample = load_prepare_sample()
    with tempfile.TemporaryDirectory() as tmp:
        project = Path(tmp)
        sample = "ST_test"
        source = project / "data" / "spatial" / sample
        source.mkdir(parents=True)

        # M1 and M2 map to the same human symbol and must be summed. M3 is
        # ambiguous (H2/H3) and must be excluded. M4 maps one-to-one to H4.
        counts = sparse.coo_matrix(
            np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype=np.int64)
        )
        matrix_path = source / "matrix.mtx.gz"
        with gzip.open(matrix_path, "wb") as handle:
            mmwrite(handle, counts)
        write_gzip_text(source / "features.tsv.gz", "g1\tM1\ng2\tM2\ng3\tM3\ng4\tM4\n")
        write_gzip_text(source / "barcodes.tsv.gz", "b-1\nc-1\n")

        annotation = project / "output" / "spatial" / sample / "Spatial_CellType.tsv"
        annotation.parent.mkdir(parents=True)
        annotation.write_text("barcode\tcell_type\nb-1\tHSC\nc-1\tTreg\n", encoding="utf-8")
        mapping = project / "mapping.tsv"
        mapping.write_text(
            "mouse_symbol\thuman_symbol\n"
            "M1\tH1\nM2\tH1\nM3\tH2\nM3\tH3\nM4\tH4\n",
            encoding="utf-8",
        )

        prepare_sample(project, sample, mapping)
        output = project / "output" / "cellphonedb" / sample
        adata = ad.read_h5ad(output / f"counts_{sample}.h5ad")

        assert adata.shape == (2, 2)
        assert list(adata.obs_names) == ["b_1", "c_1"]
        assert list(adata.var_names) == ["H1", "H4"]
        np.testing.assert_array_equal(adata.X.toarray(), np.array([[4, 7], [6, 8]]))
        meta = pd.read_csv(output / f"meta_{sample}.tsv", sep="\t")
        assert meta.to_dict("records") == [
            {"Cell": "b_1", "cell_type": "HSC"},
            {"Cell": "c_1", "cell_type": "Treg"},
        ]

    print("CellPhoneDB input-preparation regression test passed.")


if __name__ == "__main__":
    main()
