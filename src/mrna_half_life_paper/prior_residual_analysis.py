from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _load_global_table
from mrna_half_life_paper.saluki_naming import HAS_SALUKI_MOUSE_PRIOR, SALUKI_HUMAN_PC1, SALUKI_MOUSE_PRIOR

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None


PRIOR_COLS = ["mouse_pc1", SALUKI_MOUSE_PRIOR]


def _metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "pearson": float(np.corrcoef(y, pred)[0, 1]),
        "spearman": float(spearmanr(y, pred).statistic),
        "r2": float(r2_score(y, pred)),
    }


def _make_xgb(random_state: int, device: str) -> XGBRegressor:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for prior residual analysis.")
    return XGBRegressor(
        n_estimators=3000,
        early_stopping_rounds=100,
        objective="reg:squarederror",
        tree_method="hist",
        device=device,
        max_bin=256,
        max_depth=6,
        learning_rate=0.02,
        min_child_weight=4,
        subsample=0.9,
        colsample_bytree=0.75,
        reg_alpha=0.0,
        reg_lambda=1.0,
        random_state=random_state,
    )


def _prediction(model: XGBRegressor, x: np.ndarray, best_iteration: int) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
            category=UserWarning,
        )
        return model.predict(x, iteration_range=(0, best_iteration + 1))


def _load_both_prior_table(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    df, base_cols = _load_global_table(feature_path)
    both = df.loc[(df["has_mouse_pc1"] > 0.5) & (df[HAS_SALUKI_MOUSE_PRIOR] > 0.5)].copy()
    both = both.sort_values("gene_id").reset_index(drop=True)
    return both, base_cols


def run_prior_residual_analysis(
    feature_path: Path | None = None,
    results_root: Path | None = None,
    random_state: int = 42,
    cv_splits: int = 5,
    nuisance_cv_splits: int = 5,
    device: str = "cuda",
) -> dict[str, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "prior_residual_analysis"
    results_root.mkdir(parents=True, exist_ok=True)

    df, human_cols = _load_both_prior_table(feature_path)
    y = df[SALUKI_HUMAN_PC1].to_numpy(dtype=np.float32)
    x_prior = df[PRIOR_COLS].to_numpy(dtype=np.float32)
    x_human = df[human_cols].to_numpy(dtype=np.float32)

    prior_pred = np.full(len(y), np.nan, dtype=np.float32)
    compact_pred = np.full(len(y), np.nan, dtype=np.float32)
    residual_pred = np.full(len(y), np.nan, dtype=np.float32)
    residual_target_oof = np.full(len(y), np.nan, dtype=np.float32)
    fold_id = np.zeros(len(y), dtype=np.int16)
    residual_best_iterations: list[int] = []
    compact_best_iterations: list[int] = []

    kfold = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)
    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(np.arange(len(y))), start=1):
        fold_id[test_idx] = fold_idx

        train_prior_pred = np.full(len(train_idx), np.nan, dtype=np.float32)
        nuisance_kfold = KFold(
            n_splits=nuisance_cv_splits,
            shuffle=True,
            random_state=random_state + 1000 + fold_idx,
        )
        for nuisance_train_rel, nuisance_test_rel in nuisance_kfold.split(train_idx):
            nuisance_model = RidgeCV(alphas=np.logspace(-3, 3, 15))
            nuisance_model.fit(x_prior[train_idx[nuisance_train_rel]], y[train_idx[nuisance_train_rel]])
            train_prior_pred[nuisance_test_rel] = nuisance_model.predict(
                x_prior[train_idx[nuisance_test_rel]]
            ).astype(np.float32)
        if not np.isfinite(train_prior_pred).all():
            raise RuntimeError(f"Incomplete cross-fitted nuisance predictions in outer fold {fold_idx}.")

        prior_model = RidgeCV(alphas=np.logspace(-3, 3, 15))
        prior_model.fit(x_prior[train_idx], y[train_idx])
        test_prior_pred = prior_model.predict(x_prior[test_idx]).astype(np.float32)
        prior_pred[test_idx] = test_prior_pred

        train_residual = y[train_idx] - train_prior_pred
        rng = np.random.RandomState(random_state + fold_idx)
        perm = rng.permutation(len(train_idx))
        val_n = max(512, int(0.1 * len(train_idx)))
        val_rel = perm[:val_n]
        fit_rel = perm[val_n:]
        fit_idx = train_idx[fit_rel]
        val_idx = train_idx[val_rel]
        residual_by_index = pd.Series(train_residual, index=train_idx)

        residual_model = _make_xgb(random_state=random_state, device=device)
        residual_model.fit(
            x_human[fit_idx],
            residual_by_index.loc[fit_idx].to_numpy(dtype=np.float32),
            eval_set=[(x_human[val_idx], residual_by_index.loc[val_idx].to_numpy(dtype=np.float32))],
            verbose=False,
        )
        residual_best = int(getattr(residual_model, "best_iteration", 2999))
        residual_best_iterations.append(residual_best)
        residual_pred[test_idx] = _prediction(residual_model, x_human[test_idx], residual_best)
        residual_target_oof[test_idx] = y[test_idx] - test_prior_pred

        compact_model = _make_xgb(random_state=random_state, device=device)
        compact_model.fit(
            x_human[fit_idx],
            y[fit_idx],
            eval_set=[(x_human[val_idx], y[val_idx])],
            verbose=False,
        )
        compact_best = int(getattr(compact_model, "best_iteration", 2999))
        compact_best_iterations.append(compact_best)
        compact_pred[test_idx] = _prediction(compact_model, x_human[test_idx], compact_best)

    combined_pred = prior_pred + residual_pred
    prior_metrics = _metrics(y, prior_pred)
    compact_metrics = _metrics(y, compact_pred)
    combined_metrics = _metrics(y, combined_pred)
    remaining_fraction = (combined_metrics["r2"] - prior_metrics["r2"]) / (1.0 - prior_metrics["r2"])

    summary = pd.DataFrame(
        [
            {
                "analysis": "prior_residual_analysis",
                "model": "prior_only_linear",
                "n": int(len(df)),
                "cv_splits": int(cv_splits),
                "nuisance_cv_splits": int(nuisance_cv_splits),
                "random_state": int(random_state),
                "n_prior_features": int(len(PRIOR_COLS)),
                "n_human_features": 0,
                **prior_metrics,
                "delta_r2_vs_prior": 0.0,
                "remaining_variance_explained": 0.0,
                "mean_best_iter": np.nan,
            },
            {
                "analysis": "prior_residual_analysis",
                "model": "compact_all_only",
                "n": int(len(df)),
                "cv_splits": int(cv_splits),
                "nuisance_cv_splits": int(nuisance_cv_splits),
                "random_state": int(random_state),
                "n_prior_features": 0,
                "n_human_features": int(len(human_cols)),
                **compact_metrics,
                "delta_r2_vs_prior": compact_metrics["r2"] - prior_metrics["r2"],
                "remaining_variance_explained": np.nan,
                "mean_best_iter": float(np.mean(compact_best_iterations)),
            },
            {
                "analysis": "prior_residual_analysis",
                "model": "prior_plus_compact_residual",
                "n": int(len(df)),
                "cv_splits": int(cv_splits),
                "nuisance_cv_splits": int(nuisance_cv_splits),
                "random_state": int(random_state),
                "n_prior_features": int(len(PRIOR_COLS)),
                "n_human_features": int(len(human_cols)),
                **combined_metrics,
                "delta_r2_vs_prior": combined_metrics["r2"] - prior_metrics["r2"],
                "remaining_variance_explained": float(remaining_fraction),
                "mean_best_iter": float(np.mean(residual_best_iterations)),
            },
        ]
    )

    oof = df[["gene_id", SALUKI_HUMAN_PC1, "mouse_pc1", SALUKI_MOUSE_PRIOR]].copy()
    oof["fold"] = fold_id
    oof["prior_only_prediction"] = prior_pred
    oof["prior_residual_target"] = residual_target_oof
    oof["compact_all_only_prediction"] = compact_pred
    oof["compact_residual_prediction"] = residual_pred
    oof["prior_plus_compact_residual_prediction"] = combined_pred

    summary_path = results_root / "summary.tsv"
    oof_path = results_root / "oof_predictions.tsv"
    json_path = results_root / "summary.json"
    note_path = MANUSCRIPT_DIR.parent / "docs" / "prior_residual_analysis.md"
    result_table_path = results_root / "residual_decomposition_table.md"

    summary.to_csv(summary_path, sep="\t", index=False, na_rep="NA")
    oof.to_csv(oof_path, sep="\t", index=False, na_rep="NA")
    json_records = summary.astype(object).where(pd.notna(summary), None).to_dict(orient="records")
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(json_records, handle, indent=2, allow_nan=False)

    _write_note(summary, note_path)
    _write_result_table(summary, result_table_path)

    return {
        "summary": summary_path,
        "oof_predictions": oof_path,
        "summary_json": json_path,
        "note": note_path,
        "result_table": result_table_path,
    }


def _write_note(summary: pd.DataFrame, note_path: Path) -> Path:
    prior = summary.loc[summary["model"] == "prior_only_linear"].iloc[0]
    compact = summary.loc[summary["model"] == "compact_all_only"].iloc[0]
    combined = summary.loc[summary["model"] == "prior_plus_compact_residual"].iloc[0]
    lines = [
        "# Prior Residual Analysis",
        "",
        "## Purpose",
        "",
        "该分析把 mouse-prior 可解释部分与 human compact_all 特征的额外修正分开估计，",
        "用于回答 dual-prior transfer model 是否只是读取 mouse prior。",
        "",
        "## Reproducible command",
        "",
        "```bash",
        "source scripts/activate_env.sh",
        "python -m mrna_half_life_paper.prior_residual_analysis",
        "```",
        "",
        "## Algorithm",
        "",
        "1. 仅保留同时具有重建 mouse prior 和 Saluki mouse PC1 的 `both_priors_available` genes。",
        "2. 在每个 outer fold 中，用两列 mouse priors 训练 RidgeCV，并对 held-out fold 得到 `prior_only_prediction`。",
        "3. 在 outer-training genes 内再做 5-fold cross-fitting，使用 nuisance OOF prior prediction 构造 residual target。",
        "4. 用 human `compact_all` features 训练 XGBoost/CUDA 去预测 residual target。",
        "5. held-out fold 的最终预测为 `prior_only_prediction + compact_residual_prediction`。",
        "",
        "## Key results",
        "",
        f"- N={int(prior['n'])}, {int(prior['cv_splits'])}-fold OOF, random_state={int(prior['random_state'])}.",
        f"- prior-only linear: Pearson={prior['pearson']:.4f}, R2={prior['r2']:.4f}.",
        f"- compact_all only: Pearson={compact['pearson']:.4f}, R2={compact['r2']:.4f}.",
        f"- prior + compact residual: Pearson={combined['pearson']:.4f}, R2={combined['r2']:.4f}.",
        f"- Remaining variance explained after prior: {100 * combined['remaining_variance_explained']:.2f}%.",
        "",
        "## Interpretation",
        "",
        "该结果支持的边界是：mouse prior 本身很强，但 human sequence/regulatory features 仍能在 prior 之外提供额外修正；",
        "它不能被写成 pure sequence-only result，也不应被解释为 prior 完全没有贡献。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def _write_result_table(summary: pd.DataFrame, table_path: Path) -> Path:
    rows = [
        "# Prior residual decomposition result table",
        "",
        "| model | genes | prior_features | human_features | pearson | spearman | r2 | delta_r2_vs_prior | remaining_variance_explained |",
        "|:--|--:|--:|--:|--:|--:|--:|--:|--:|",
    ]
    for _, row in summary.iterrows():
        remaining = "" if pd.isna(row["remaining_variance_explained"]) else f"{100 * row['remaining_variance_explained']:.2f}%"
        rows.append(
            "| "
            + " | ".join(
                [
                    str(row["model"]),
                    str(int(row["n"])),
                    str(int(row["n_prior_features"])),
                    str(int(row["n_human_features"])),
                    f"{row['pearson']:.3f}",
                    f"{row['spearman']:.3f}",
                    f"{row['r2']:.3f}",
                    f"{row['delta_r2_vs_prior']:.3f}",
                    remaining,
                ]
            )
            + " |"
        )
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return table_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Separate mouse-prior signal from human compact feature residual correction."
    )
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "prior_residual_analysis")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--cv-splits", type=int, default=5)
    parser.add_argument("--nuisance-cv-splits", type=int, default=5)
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_prior_residual_analysis(
        feature_path=args.feature_path,
        results_root=args.results_root,
        random_state=args.random_state,
        cv_splits=args.cv_splits,
        nuisance_cv_splits=args.nuisance_cv_splits,
        device=args.device,
    )
    for label, path in outputs.items():
        print(f"{label}: {path}")
    print(pd.read_csv(outputs["summary"], sep="\t").to_string(index=False))


if __name__ == "__main__":
    main()
