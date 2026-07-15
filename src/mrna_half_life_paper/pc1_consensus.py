from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def compute_mean_consensus(matrix: pd.DataFrame) -> pd.Series:
    return matrix.mean(axis=1).rename("mean_consensus")


def _align_pc1_sign(pc1: pd.Series, reference: pd.Series) -> pd.Series:
    corr = np.corrcoef(pc1.to_numpy(), reference.to_numpy())[0, 1]
    if np.isnan(corr):
        corr = 1.0
    if corr < 0:
        pc1 = -pc1
    return pc1


def compute_pc1_consensus(matrix: pd.DataFrame) -> pd.Series:
    pca = PCA(n_components=1, svd_solver="full", random_state=0)
    scores = pca.fit_transform(matrix.to_numpy(dtype=float)).ravel()
    pc1 = pd.Series(scores, index=matrix.index, name="pc1_consensus")
    return _align_pc1_sign(pc1, matrix.mean(axis=1))


def compute_study_weighted_pc1(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "study",
) -> pd.Series:
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include 'sample_id'.")
    if group_col not in metadata.columns:
        raise ValueError(f"Metadata must include '{group_col}'.")

    sample_order = list(matrix.columns)
    aligned_meta = metadata.set_index("sample_id").loc[sample_order]
    group_sizes = aligned_meta[group_col].value_counts()
    weights = aligned_meta[group_col].map(lambda value: 1.0 / np.sqrt(group_sizes[value])).to_numpy(dtype=float)
    weights = weights / weights.mean()

    weighted_matrix = matrix.multiply(weights, axis=1)
    pca = PCA(n_components=1, svd_solver="full", random_state=0)
    scores = pca.fit_transform(weighted_matrix.to_numpy(dtype=float)).ravel()
    pc1 = pd.Series(scores, index=matrix.index, name="study_weighted_pc1")
    return _align_pc1_sign(pc1, weighted_matrix.mean(axis=1))


def compute_group_collapsed_pc1(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "study",
) -> pd.Series:
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include 'sample_id'.")
    if group_col not in metadata.columns:
        raise ValueError(f"Metadata must include '{group_col}'.")

    aligned_meta = metadata.set_index("sample_id").loc[list(matrix.columns)]
    grouped = aligned_meta[group_col]
    collapsed_cols: dict[str, pd.Series] = {}
    for group_value in grouped.unique().tolist():
        sample_ids = grouped.index[grouped == group_value].tolist()
        collapsed_cols[str(group_value)] = matrix.loc[:, sample_ids].mean(axis=1)
    collapsed = pd.DataFrame(collapsed_cols, index=matrix.index)
    pc1 = compute_pc1_consensus(collapsed)
    pc1.name = "group_collapsed_pc1"
    return pc1


def compute_jackknife_pc1(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "study",
) -> pd.Series:
    if "sample_id" not in metadata.columns:
        raise ValueError("Metadata must include 'sample_id'.")
    if group_col not in metadata.columns:
        raise ValueError(f"Metadata must include '{group_col}'.")

    base_pc1 = compute_pc1_consensus(matrix).rename("base_pc1")
    labels: list[pd.Series] = []
    groups = sorted(metadata[group_col].unique().tolist())
    if len(groups) < 2:
        raise ValueError("Jackknife PC1 requires at least two groups.")

    for group_value in groups:
        keep_meta = metadata.loc[metadata[group_col] != group_value].copy()
        keep_samples = keep_meta["sample_id"].tolist()
        subset_matrix = matrix.loc[:, keep_samples]
        subset_pc1 = compute_pc1_consensus(subset_matrix)
        subset_pc1 = _align_pc1_sign(subset_pc1, base_pc1.loc[subset_pc1.index])
        std = float(subset_pc1.std(ddof=0))
        if std == 0 or np.isnan(std):
            standardized = subset_pc1 - float(subset_pc1.mean())
        else:
            standardized = (subset_pc1 - float(subset_pc1.mean())) / std
        labels.append(standardized.rename(str(group_value)))

    jackknife = pd.concat(labels, axis=1).mean(axis=1)
    jackknife.name = "jackknife_pc1"
    return _align_pc1_sign(jackknife, base_pc1.loc[jackknife.index])
