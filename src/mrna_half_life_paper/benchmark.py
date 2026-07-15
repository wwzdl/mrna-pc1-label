from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def label_correlation(a: pd.Series, b: pd.Series) -> dict[str, float]:
    common = a.index.intersection(b.index)
    if len(common) < 2:
        return {"n": float(len(common)), "pearson": np.nan, "spearman": np.nan}
    a_common = a.loc[common].astype(float)
    b_common = b.loc[common].astype(float)
    pearson = float(a_common.corr(b_common, method="pearson"))
    spearman = float(spearmanr(a_common, b_common, nan_policy="omit").statistic)
    return {"n": float(len(common)), "pearson": pearson, "spearman": spearman}


def leave_one_group_out(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    estimator_name: str,
    estimator_fn: Callable[[pd.DataFrame, pd.DataFrame], pd.Series],
) -> pd.DataFrame:
    full = estimator_fn(matrix, metadata)
    rows: list[dict[str, float | str]] = []

    for group_value in sorted(metadata[group_col].unique().tolist()):
        keep_meta = metadata.loc[metadata[group_col] != group_value].copy()
        keep_samples = keep_meta["sample_id"].tolist()
        subset_matrix = matrix.loc[:, keep_samples]
        subset_labels = estimator_fn(subset_matrix, keep_meta)
        corr = label_correlation(full, subset_labels)
        rows.append(
            {
                "estimator": estimator_name,
                "group_type": group_col,
                "left_out": group_value,
                "n_genes": corr["n"],
                "pearson": corr["pearson"],
                "spearman": corr["spearman"],
            }
        )

    return pd.DataFrame(rows)
