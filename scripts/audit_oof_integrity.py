#!/usr/bin/env python3
"""Verify coverage, uniqueness, and completeness of the key OOF predictions."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def _audit_group(label: str, relative_root: str, pattern: str, expected_n: int) -> None:
    paths = sorted((ROOT / relative_root).glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No OOF files found for {label}: {relative_root}/{pattern}")

    reference_ids: tuple[str, ...] | None = None
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        if len(frame) != expected_n:
            raise ValueError(f"{path}: expected {expected_n} rows, found {len(frame)}")
        if frame["gene_id"].nunique() != expected_n:
            raise ValueError(f"{path}: gene_id values are not unique")

        prediction_cols = [
            column
            for column in frame.columns
            if column == "prediction" or column.startswith("prediction_trained_on_")
        ]
        if not prediction_cols:
            raise ValueError(f"{path}: no prediction column found")
        if frame[prediction_cols].isna().any().any():
            raise ValueError(f"{path}: prediction columns contain missing values")

        gene_ids = tuple(frame["gene_id"].astype(str))
        if reference_ids is None:
            reference_ids = gene_ids
        elif gene_ids != reference_ids:
            raise ValueError(f"{path}: gene order or evaluation universe differs within {label}")

    print(f"[PASS] {label}: {len(paths)} files x {expected_n:,} unique genes; predictions complete")


def main() -> None:
    _audit_group(
        "global fixed-target",
        "results/tenfold_global_prior/oof_predictions",
        "*.tsv",
        12_916,
    )
    _audit_group(
        "target-shrinkage cross-target",
        "results/ortholog_regularized_cross_target/oof_predictions",
        "*.tsv.gz",
        12_307,
    )
    _audit_group(
        "target-shrinkage benchmark",
        "results/ortholog_regularized_label_10fold_main/oof_predictions",
        "*.tsv",
        12_307,
    )

    residual_path = ROOT / "results/prior_residual_analysis/oof_predictions.tsv"
    residual = pd.read_csv(residual_path, sep="\t")
    if len(residual) != 11_107 or residual["gene_id"].nunique() != 11_107:
        raise ValueError(f"{residual_path}: expected 11,107 unique genes")
    prediction_cols = [
        "prior_only_prediction",
        "compact_all_only_prediction",
        "compact_residual_prediction",
        "prior_plus_compact_residual_prediction",
    ]
    if residual[prediction_cols].isna().any().any():
        raise ValueError(f"{residual_path}: prediction components contain missing values")
    if sorted(residual["fold"].unique().tolist()) != [1, 2, 3, 4, 5]:
        raise ValueError(f"{residual_path}: expected outer folds 1-5")
    print("[PASS] residual decomposition: 11,107 unique genes; folds 1-5 and all components complete")


if __name__ == "__main__":
    main()
