from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from mrna_half_life_paper.config import RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _load_compact_all_features
from mrna_half_life_paper.label_competition import _xgb_oof_predictions
from mrna_half_life_paper.ortholog_regularized_label import (
    _load_human_base_label,
    _load_reconstructed_mouse_pc1,
    _load_saluki_label,
    _zscore,
    make_ortholog_regularized_label,
)


HUMAN_TARGET = "human_no_gejman_pc1"
ORTHOREG_TARGET = "orthoreg_reconstructed_mouse_pc1_0.10"
SALUKI_REFERENCE = "saluki_human_pc1"


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    mse = mean_squared_error(y_true, y_pred)
    return {
        "pearson": float(np.corrcoef(y_true, y_pred)[0, 1]),
        "spearman": float(spearmanr(y_true, y_pred).statistic),
        "r2": float(r2_score(y_true, y_pred)),
        "rmse": float(np.sqrt(mse)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
    }


def _bootstrap_correlation_difference(
    y_true: np.ndarray,
    prediction_a: np.ndarray,
    prediction_b: np.ndarray,
    *,
    n_bootstrap: int,
    random_state: int,
) -> dict[str, float | int]:
    """Paired gene bootstrap for Pearson(A, y) - Pearson(B, y)."""

    rng = np.random.default_rng(random_state)
    n = y_true.shape[0]
    observed_a = float(np.corrcoef(y_true, prediction_a)[0, 1])
    observed_b = float(np.corrcoef(y_true, prediction_b)[0, 1])
    differences = np.empty(n_bootstrap, dtype=float)
    for idx in range(n_bootstrap):
        sample = rng.integers(0, n, size=n)
        corr_a = np.corrcoef(y_true[sample], prediction_a[sample])[0, 1]
        corr_b = np.corrcoef(y_true[sample], prediction_b[sample])[0, 1]
        differences[idx] = corr_a - corr_b
    lower, upper = np.quantile(differences, [0.025, 0.975])
    return {
        "n": int(n),
        "n_bootstrap": int(n_bootstrap),
        "pearson_a": observed_a,
        "pearson_b": observed_b,
        "difference_a_minus_b": observed_a - observed_b,
        "difference_ci95_low": float(lower),
        "difference_ci95_high": float(upper),
        "bootstrap_fraction_le_zero": float(np.mean(differences <= 0.0)),
    }


def _bootstrap_label_correlation(
    label_a: np.ndarray,
    label_b: np.ndarray,
    *,
    n_bootstrap: int,
    random_state: int,
) -> dict[str, float | int]:
    rng = np.random.default_rng(random_state)
    n = label_a.shape[0]
    observed = float(np.corrcoef(label_a, label_b)[0, 1])
    values = np.empty(n_bootstrap, dtype=float)
    for idx in range(n_bootstrap):
        sample = rng.integers(0, n, size=n)
        values[idx] = np.corrcoef(label_a[sample], label_b[sample])[0, 1]
    lower, upper = np.quantile(values, [0.025, 0.975])
    return {
        "n": int(n),
        "n_bootstrap": int(n_bootstrap),
        "pearson": observed,
        "pearson_ci95_low": float(lower),
        "pearson_ci95_high": float(upper),
    }


def run_cross_target_analysis(
    *,
    feature_path: Path,
    results_root: Path,
    random_states: list[int],
    cv_splits: int,
    device: str | None,
    n_bootstrap: int,
) -> dict[str, Path]:
    results_root.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    features = _load_compact_all_features(feature_path).sort_values("gene_id").reset_index(drop=True)
    feature_cols = [column for column in features.columns if column != "gene_id"]
    human = _zscore(_load_human_base_label()).rename(HUMAN_TARGET)
    mouse = _load_reconstructed_mouse_pc1()
    saluki = _zscore(_load_saluki_label("human")).rename(SALUKI_REFERENCE)
    orthoreg = make_ortholog_regularized_label(
        human_label=human,
        mouse_label=mouse,
        lambda_=0.10,
        saluki_reference=saluki,
    ).rename(ORTHOREG_TARGET)

    target_frame = pd.concat([human, orthoreg, saluki], axis=1).dropna().reset_index().rename(columns={"index": "gene_id"})
    merged = (
        features.merge(target_frame, on="gene_id", how="inner")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )

    summary_rows: list[dict[str, float | int | str]] = []
    seed_predictions: list[pd.DataFrame] = []
    for random_state in random_states:
        prediction_frame = merged[["gene_id", HUMAN_TARGET, ORTHOREG_TARGET, SALUKI_REFERENCE]].copy()
        for training_target in (HUMAN_TARGET, ORTHOREG_TARGET):
            prediction, best_iterations, device_used = _xgb_oof_predictions(
                merged,
                feature_cols,
                training_target,
                random_state=random_state,
                cv_splits=cv_splits,
                device=device,
            )
            prediction_col = f"prediction_trained_on_{training_target}"
            prediction_frame[prediction_col] = prediction
            for evaluation_target in (HUMAN_TARGET, ORTHOREG_TARGET, SALUKI_REFERENCE):
                metrics = _metrics(
                    merged[evaluation_target].to_numpy(dtype=float),
                    prediction.astype(float),
                )
                summary_rows.append(
                    {
                        "training_target": training_target,
                        "evaluation_target": evaluation_target,
                        "random_state": int(random_state),
                        "cv_splits": int(cv_splits),
                        "n": int(len(merged)),
                        "n_features": int(len(feature_cols)),
                        "mean_best_iteration": float(np.mean(best_iterations)),
                        "device_used": device_used,
                        **metrics,
                    }
                )
        prediction_frame.to_csv(
            oof_dir / f"cross_target__folds{cv_splits}__seed{random_state}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        seed_predictions.append(prediction_frame)

    summary_by_seed = pd.DataFrame(summary_rows)
    summary_by_seed_path = results_root / "summary_by_seed.tsv"
    summary_by_seed.to_csv(summary_by_seed_path, sep="\t", index=False)

    metric_cols = ["pearson", "spearman", "r2", "rmse", "mae", "mean_best_iteration"]
    aggregate = (
        summary_by_seed.groupby(["training_target", "evaluation_target"], sort=False)
        .agg(
            n_random_states=("random_state", "nunique"),
            cv_splits=("cv_splits", "first"),
            n=("n", "first"),
            n_features=("n_features", "first"),
            **{f"{metric}_mean": (metric, "mean") for metric in metric_cols},
            **{f"{metric}_sd": (metric, "std") for metric in metric_cols},
        )
        .reset_index()
    )
    aggregate_path = results_root / "summary_aggregate.tsv"
    aggregate.to_csv(aggregate_path, sep="\t", index=False)

    ensemble = seed_predictions[0][["gene_id", HUMAN_TARGET, ORTHOREG_TARGET, SALUKI_REFERENCE]].copy()
    for training_target in (HUMAN_TARGET, ORTHOREG_TARGET):
        prediction_col = f"prediction_trained_on_{training_target}"
        ensemble[prediction_col] = np.mean(
            [frame[prediction_col].to_numpy(dtype=float) for frame in seed_predictions],
            axis=0,
        )
    ensemble_path = results_root / "oof_three_seed_ensemble.tsv.gz"
    ensemble.to_csv(ensemble_path, sep="\t", index=False, compression="gzip")

    human_prediction = ensemble[f"prediction_trained_on_{HUMAN_TARGET}"].to_numpy(dtype=float)
    orthoreg_prediction = ensemble[f"prediction_trained_on_{ORTHOREG_TARGET}"].to_numpy(dtype=float)
    bootstrap_rows: list[dict[str, float | int | str]] = []
    for evaluation_target in (HUMAN_TARGET, ORTHOREG_TARGET, SALUKI_REFERENCE):
        result = _bootstrap_correlation_difference(
            ensemble[evaluation_target].to_numpy(dtype=float),
            orthoreg_prediction,
            human_prediction,
            n_bootstrap=n_bootstrap,
            random_state=20260713,
        )
        bootstrap_rows.append(
            {
                "analysis": "paired_prediction_comparison",
                "evaluation_target": evaluation_target,
                "model_a": f"trained_on_{ORTHOREG_TARGET}",
                "model_b": f"trained_on_{HUMAN_TARGET}",
                **result,
            }
        )
    label_result = _bootstrap_label_correlation(
        ensemble[ORTHOREG_TARGET].to_numpy(dtype=float),
        ensemble[HUMAN_TARGET].to_numpy(dtype=float),
        n_bootstrap=n_bootstrap,
        random_state=20260714,
    )
    bootstrap_rows.append(
        {
            "analysis": "label_geometry",
            "evaluation_target": HUMAN_TARGET,
            "model_a": ORTHOREG_TARGET,
            "model_b": HUMAN_TARGET,
            **label_result,
        }
    )
    bootstrap_path = results_root / "paired_bootstrap_ci.tsv"
    pd.DataFrame(bootstrap_rows).to_csv(bootstrap_path, sep="\t", index=False)

    metadata = {
        "analysis": "human-only cross-target evaluation of weak ortholog regularization",
        "human_target": HUMAN_TARGET,
        "ortholog_regularized_target": ORTHOREG_TARGET,
        "lambda": 0.10,
        "mouse_label_source": "locally reconstructed mouse PC1 under Saluki-like coverage",
        "features": "compact_all human sequence/regulatory features only",
        "n_features": len(feature_cols),
        "n_genes": len(merged),
        "cv_splits": cv_splits,
        "random_states": random_states,
        "bootstrap_iterations": n_bootstrap,
        "prediction_boundary": "No mouse prior values are used as model inputs in this analysis.",
    }
    metadata_path = results_root / "analysis_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return {
        "summary_by_seed": summary_by_seed_path,
        "summary_aggregate": aggregate_path,
        "bootstrap": bootstrap_path,
        "ensemble": ensemble_path,
        "metadata": metadata_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cross-target OOF evaluation for the lambda=0.10 ortholog-regularized label."
    )
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR / "ortholog_regularized_cross_target",
    )
    parser.add_argument("--random-states", nargs="+", type=int, default=[13, 42, 101])
    parser.add_argument("--cv-splits", type=int, default=10)
    parser.add_argument("--device", default=None)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_cross_target_analysis(
        feature_path=args.feature_path,
        results_root=args.results_root,
        random_states=args.random_states,
        cv_splits=args.cv_splits,
        device=args.device,
        n_bootstrap=args.bootstrap_iterations,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
