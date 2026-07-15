from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


@dataclass
class PreprocessConfig:
    log_offset: float = 0.1
    apply_log: bool = True
    min_observed_per_gene: int = 3
    n_components: int | None = None
    max_iter: int = 30
    tol: float = 1e-6


def safe_log_transform(matrix: pd.DataFrame, offset: float = 0.1) -> pd.DataFrame:
    values = matrix.astype(float)
    min_value = np.nanmin(values.to_numpy())
    if min_value <= -offset:
        raise ValueError(
            f"Matrix contains values <= {-offset}; cannot safely apply log with offset {offset}."
        )
    return np.log(values + offset)


def filter_genes_by_coverage(matrix: pd.DataFrame, min_observed_per_gene: int) -> pd.DataFrame:
    keep = matrix.notna().sum(axis=1) >= min_observed_per_gene
    return matrix.loc[keep].copy()


def zscore_by_sample(matrix: pd.DataFrame) -> pd.DataFrame:
    centered = matrix.subtract(matrix.mean(axis=0, skipna=True), axis=1)
    scale = matrix.std(axis=0, skipna=True, ddof=0).replace(0, 1.0)
    return centered.divide(scale, axis=1)


def iterative_pca_impute(
    matrix: pd.DataFrame,
    n_components: int | None = None,
    max_iter: int = 30,
    tol: float = 1e-6,
) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    mask = np.isnan(values)
    if not mask.any():
        return matrix.copy()

    filled = values.copy()
    col_means = np.nanmean(filled, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means)
    filled[mask] = np.take(col_means, np.where(mask)[1])

    n_rows, n_cols = filled.shape
    max_rank = max(1, min(n_rows, n_cols) - 1)
    k = n_components if n_components is not None else min(5, max_rank)
    k = max(1, min(k, max_rank))

    prev = filled.copy()
    for _ in range(max_iter):
        pca = PCA(n_components=k, svd_solver="full", random_state=0)
        scores = pca.fit_transform(filled)
        reconstructed = pca.inverse_transform(scores)
        filled[mask] = reconstructed[mask]
        delta = np.nanmax(np.abs(filled - prev))
        prev = filled.copy()
        if float(delta) < tol:
            break

    return pd.DataFrame(filled, index=matrix.index, columns=matrix.columns)


def quantile_normalize(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    sorter = np.argsort(values, axis=0)
    sorted_vals = np.take_along_axis(values, sorter, axis=0)
    mean_sorted = sorted_vals.mean(axis=1)

    ranks = np.argsort(sorter, axis=0)
    normalized = np.take_along_axis(mean_sorted[:, None], ranks, axis=0)
    return pd.DataFrame(normalized, index=matrix.index, columns=matrix.columns)


def preprocess_matrix(matrix: pd.DataFrame, config: PreprocessConfig) -> pd.DataFrame:
    out = matrix.copy()
    out = filter_genes_by_coverage(out, config.min_observed_per_gene)
    if config.apply_log:
        out = safe_log_transform(out, offset=config.log_offset)
    out = zscore_by_sample(out)
    out = iterative_pca_impute(
        out,
        n_components=config.n_components,
        max_iter=config.max_iter,
        tol=config.tol,
    )
    out = quantile_normalize(out)
    return out
