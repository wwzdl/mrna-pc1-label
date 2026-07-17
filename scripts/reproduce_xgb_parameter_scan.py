#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from xgboost import XGBRegressor

from mrna_half_life_paper.config import RESULTS_DIR
from mrna_half_life_paper.saluki_naming import SALUKI_HUMAN_PC1
from mrna_half_life_paper.sequence_baseline import _load_target_tables


SCAN_CONFIGS: tuple[dict[str, float | int], ...] = (
    {
        "config_id": 1,
        "max_depth": 6,
        "learning_rate": 0.02,
        "min_child_weight": 4,
        "subsample": 0.90,
        "colsample_bytree": 0.75,
        "reg_alpha": 0.0,
        "reg_lambda": 1.0,
    },
    {
        "config_id": 2,
        "max_depth": 6,
        "learning_rate": 0.03,
        "min_child_weight": 3,
        "subsample": 0.90,
        "colsample_bytree": 0.70,
        "reg_alpha": 0.0,
        "reg_lambda": 1.0,
    },
    {
        "config_id": 3,
        "max_depth": 6,
        "learning_rate": 0.03,
        "min_child_weight": 5,
        "subsample": 0.90,
        "colsample_bytree": 0.80,
        "reg_alpha": 0.0,
        "reg_lambda": 1.5,
    },
    {
        "config_id": 4,
        "max_depth": 6,
        "learning_rate": 0.03,
        "min_child_weight": 4,
        "subsample": 0.85,
        "colsample_bytree": 0.70,
        "reg_alpha": 0.0,
        "reg_lambda": 1.0,
    },
    {
        "config_id": 5,
        "max_depth": 7,
        "learning_rate": 0.03,
        "min_child_weight": 4,
        "subsample": 0.85,
        "colsample_bytree": 0.70,
        "reg_alpha": 0.0,
        "reg_lambda": 1.0,
    },
)


def _load_scan_table(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    features = pd.read_csv(feature_path, sep="\t")
    targets = _load_target_tables("human")
    common_gene_ids = set(features["gene_id"].astype(str))
    for target in targets.values():
        common_gene_ids &= set(target["gene_id"].astype(str))

    target = targets[SALUKI_HUMAN_PC1]
    merged = (
        features.loc[features["gene_id"].isin(common_gene_ids)]
        .merge(target, on="gene_id", how="inner")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    feature_cols = [
        column
        for column in merged.columns
        if column not in {"gene_id", "transcript_id", SALUKI_HUMAN_PC1}
    ]
    return merged, feature_cols


def _run_config(
    frame: pd.DataFrame,
    feature_cols: list[str],
    params: dict[str, float | int],
    *,
    device: str,
    cv_splits: int,
    random_state: int,
) -> dict[str, float | int]:
    x = frame[feature_cols].to_numpy(dtype=np.float32)
    y = frame[SALUKI_HUMAN_PC1].to_numpy(dtype=np.float32)
    prediction = np.full(len(y), np.nan, dtype=np.float32)
    best_iterations: list[int] = []
    kfold = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)
    started = time.perf_counter()

    for fold_idx, (outer_train, outer_test) in enumerate(kfold.split(x), start=1):
        rng = np.random.RandomState(random_state + fold_idx)
        permutation = rng.permutation(len(outer_train))
        validation_size = max(512, int(0.1 * len(outer_train)))
        validation = outer_train[permutation[:validation_size]]
        fit = outer_train[permutation[validation_size:]]

        model = XGBRegressor(
            n_estimators=5000,
            early_stopping_rounds=120,
            objective="reg:squarederror",
            tree_method="hist",
            device=device,
            max_bin=256,
            max_depth=int(params["max_depth"]),
            learning_rate=float(params["learning_rate"]),
            min_child_weight=float(params["min_child_weight"]),
            subsample=float(params["subsample"]),
            colsample_bytree=float(params["colsample_bytree"]),
            reg_alpha=float(params["reg_alpha"]),
            reg_lambda=float(params["reg_lambda"]),
            random_state=random_state,
        )
        model.fit(
            x[fit],
            y[fit],
            eval_set=[(x[validation], y[validation])],
            verbose=False,
        )
        best_iteration = int(getattr(model, "best_iteration", 4999))
        best_iterations.append(best_iteration)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                category=UserWarning,
            )
            prediction[outer_test] = model.predict(
                x[outer_test],
                iteration_range=(0, best_iteration + 1),
            )

    return {
        **params,
        "pearson": float(np.corrcoef(y, prediction)[0, 1]),
        "spearman": float(spearmanr(y, prediction).statistic),
        "r2": float(r2_score(y, prediction)),
        "mean_best_iter": float(np.mean(best_iterations)),
        "sec": float(time.perf_counter() - started),
    }


def run_scan(
    *,
    feature_path: Path,
    output_path: Path,
    metadata_path: Path,
    device: str,
    cv_splits: int,
    random_state: int,
) -> tuple[Path, Path]:
    frame, feature_cols = _load_scan_table(feature_path)
    rows = [
        _run_config(
            frame,
            feature_cols,
            params,
            device=device,
            cv_splits=cv_splits,
            random_state=random_state,
        )
        for params in SCAN_CONFIGS
    ]
    result = pd.DataFrame(rows).sort_values(["pearson", "r2"], ascending=False).reset_index(drop=True)
    result = result.drop(columns=["config_id"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    selected = result.iloc[0]
    metadata = {
        "scope": "Preliminary base-sequence XGBoost scan on the P2 common-gene universe.",
        "target": SALUKI_HUMAN_PC1,
        "genes": int(len(frame)),
        "features": int(len(feature_cols)),
        "cv_splits": int(cv_splits),
        "random_state": int(random_state),
        "inner_validation": "max(512, 10% of outer-training genes)",
        "selection_metric": "highest OOF Pearson, with OOF R2 as a deterministic tie-breaker",
        "selected_parameters": {
            key: int(selected[key]) if key in {"max_depth", "min_child_weight"} else float(selected[key])
            for key in (
                "max_depth",
                "learning_rate",
                "min_child_weight",
                "subsample",
                "colsample_bytree",
                "reg_alpha",
                "reg_lambda",
            )
        },
        "downstream_use": (
            "The selected setting was frozen for compact, prior-enhanced, shuffled-prior, "
            "cross-target, and residual analyses."
        ),
        "inference_boundary": (
            "Selection and OOF evaluation use the same compendium, so absolute performance is "
            "descriptive rather than a nested-CV estimate. Matched real-versus-control comparisons "
            "use identical frozen parameters."
        ),
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return output_path, metadata_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reproduce the preliminary XGBoost parameter scan.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS_DIR / "human_xgboost_gpu_saluki_refine.tsv",
    )
    parser.add_argument(
        "--metadata-output",
        type=Path,
        default=RESULTS_DIR / "human_xgboost_parameter_provenance.json",
    )
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--cv-splits", type=int, default=5)
    parser.add_argument("--random-state", type=int, default=42)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output, metadata = run_scan(
        feature_path=args.feature_path,
        output_path=args.output,
        metadata_path=args.metadata_output,
        device=args.device,
        cv_splits=args.cv_splits,
        random_state=args.random_state,
    )
    print(output)
    print(metadata)
    print(pd.read_csv(output, sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
