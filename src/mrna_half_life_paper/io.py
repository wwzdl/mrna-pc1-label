from __future__ import annotations

from pathlib import Path

import pandas as pd


REQUIRED_METADATA_COLUMNS = ["sample_id", "study", "method", "cell_type", "species"]


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"Unsupported table format: {path}")


def load_matrix(path: str | Path) -> pd.DataFrame:
    df = read_table(path)
    if df.shape[1] < 2:
        raise ValueError("Matrix file must contain at least one gene column and one sample column.")
    gene_col = df.columns[0]
    df = df.set_index(gene_col)
    return df.apply(pd.to_numeric, errors="coerce")


def load_metadata(path: str | Path) -> pd.DataFrame:
    meta = read_table(path)
    missing = [col for col in REQUIRED_METADATA_COLUMNS if col not in meta.columns]
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")
    meta = meta.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    return meta


def align_matrix_and_metadata(matrix: pd.DataFrame, metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_ids = metadata["sample_id"].tolist()
    missing_in_matrix = [sample for sample in sample_ids if sample not in matrix.columns]
    if missing_in_matrix:
        raise ValueError(f"Samples declared in metadata are missing from matrix: {missing_in_matrix[:10]}")

    common_samples = [sample for sample in matrix.columns if sample in set(sample_ids)]
    if not common_samples:
        raise ValueError("No shared samples found between matrix columns and metadata.")

    aligned_matrix = matrix.loc[:, common_samples].copy()
    aligned_metadata = metadata.set_index("sample_id").loc[common_samples].reset_index()
    return aligned_matrix, aligned_metadata


def save_table(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    df.to_csv(path, sep=sep, index=True)
