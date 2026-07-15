from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _build_global_settings, _load_compact_all_features
from mrna_half_life_paper.moesm3_saluki_label_benchmark import TARGET_SPECS, _read_target

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None


DEFAULT_TARGETS = [
    "saluki_human_pc1",
    "moesm3_all_samples_as_is_pc1",
    "moesm3_no_gejman_as_is_pc1",
    "moesm3_no_gejman_rerun_pipeline_pc1",
]


def _metrics(y_true: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "pearson": float(np.corrcoef(y_true, pred)[0, 1]),
        "spearman": float(spearmanr(y_true, pred).statistic),
        "r2": float(r2_score(y_true, pred)),
    }


def _correlation(y_left: np.ndarray, y_right: np.ndarray) -> dict[str, float]:
    return {
        "pearson": float(np.corrcoef(y_left, y_right)[0, 1]),
        "spearman": float(spearmanr(y_left, y_right).statistic),
    }


def _load_target_tables(target_names: list[str]) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for target_name in target_names:
        if target_name not in TARGET_SPECS:
            raise ValueError(f"Unknown target: {target_name}")
        tables[target_name] = _read_target(TARGET_SPECS[target_name], target_name)
    return tables


def _load_global_feature_table(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    features = _load_compact_all_features(feature_path)

    mouse_common = pd.read_csv(
        RESULTS_DIR / "ortholog_analysis" / "common_ortholog_pairs.tsv",
        sep="\t",
    )[["human_gene_id", "mouse_pc1", "is_high_confidence"]].rename(columns={"human_gene_id": "gene_id"})
    saluki_mouse = pd.read_csv(
        RESULTS_DIR / "ortholog_analysis" / "saluki_ortholog_pairs.tsv",
        sep="\t",
    )[["human_gene_id", "saluki_mouse_prior", "is_high_confidence"]].rename(
        columns={"human_gene_id": "gene_id", "is_high_confidence": "saluki_mouse_high_confidence"}
    )

    df = features.merge(mouse_common, on="gene_id", how="left").merge(saluki_mouse, on="gene_id", how="left")
    df["has_mouse_pc1"] = df["mouse_pc1"].notna().astype(np.float32)
    df["mouse_pc1_highconf"] = df["is_high_confidence"].fillna(0).astype(np.float32)
    df["mouse_pc1"] = df["mouse_pc1"].fillna(0.0)
    df["has_saluki_mouse_prior"] = df["saluki_mouse_prior"].notna().astype(np.float32)
    df["saluki_mouse_high_confidence"] = df["saluki_mouse_high_confidence"].fillna(0).astype(np.float32)
    df["saluki_mouse_prior"] = df["saluki_mouse_prior"].fillna(0.0)
    df = df.drop(columns=["is_high_confidence"]).sort_values("gene_id").reset_index(drop=True)

    feature_cols = [col for col in features.columns if col != "gene_id"]
    return df, feature_cols


def _xgb_oof_predictions(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    random_state: int = 42,
    cv_splits: int = 5,
) -> tuple[np.ndarray, list[int]]:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for the MOESM3 Saluki label prior benchmark")

    x = df[feature_cols].to_numpy(dtype=np.float32)
    y = df[target_col].to_numpy(dtype=np.float32)
    predictions = np.full(len(y), np.nan, dtype=np.float32)
    best_iterations: list[int] = []
    kfold = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)

    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(np.arange(len(y))), start=1):
        rng = np.random.RandomState(random_state + fold_idx)
        perm = rng.permutation(train_idx.size)
        val_n = max(512, int(0.1 * train_idx.size))
        val_rel = perm[:val_n]
        fit_rel = perm[val_n:]
        fit_idx = train_idx[fit_rel]
        val_idx = train_idx[val_rel]

        model = XGBRegressor(
            n_estimators=3000,
            early_stopping_rounds=100,
            objective="reg:squarederror",
            tree_method="hist",
            device="cuda",
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
        model.fit(
            x[fit_idx],
            y[fit_idx],
            eval_set=[(x[val_idx], y[val_idx])],
            verbose=False,
        )
        best_iteration = int(getattr(model, "best_iteration", 2999))
        best_iterations.append(best_iteration)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                category=UserWarning,
            )
            predictions[test_idx] = model.predict(x[test_idx], iteration_range=(0, best_iteration + 1))
    return predictions, best_iterations


def write_benchmark_note(summary: pd.DataFrame, note_path: Path) -> Path:
    by_setting = summary.loc[summary["target"] == "saluki_human_pc1"].sort_values(
        ["pearson_vs_saluki_reference"], ascending=False
    )
    baseline = by_setting.loc[by_setting["setting"] == "compact_all_global"].iloc[0]
    best_prior = by_setting.loc[by_setting["setting"] == "compact_all_plus_both_mouse_priors_global"].iloc[0]
    best_overall = summary.sort_values(["pearson_vs_saluki_reference", "r2_vs_saluki_reference"], ascending=False).iloc[0]
    best_non_saluki = summary.loc[summary["target"] != "saluki_human_pc1"].sort_values(
        ["pearson_vs_saluki_reference", "r2_vs_saluki_reference"], ascending=False
    ).iloc[0]

    lines = [
        "# MOESM3 Saluki Label Prior Benchmark",
        "",
        "更新时间：2026-05-16",
        "",
        "## 目的",
        "",
        "测试 `MOESM3` Saluki 矩阵衍生标签在带 `mouse prior` 的 cross-species transfer setting 下，",
        "是否比直接训练 Saluki 发布的 `saluki_human_pc1` 更适合作为 human 主任务标签。",
        "",
        "## 关键结果",
        "",
        f"- `saluki_human_pc1 + compact_all_global`：Pearson(Saluki ref)={baseline['pearson_vs_saluki_reference']:.4f}。",
        f"- `saluki_human_pc1 + compact_all_plus_both_mouse_priors_global`：Pearson(Saluki ref)={best_prior['pearson_vs_saluki_reference']:.4f}。",
        f"- 全部组合里，对 Saluki `saluki_human_pc1` 恢复最好的方案是 "
        f"`{best_overall['target']} + {best_overall['setting']}`，"
        f"Pearson={best_overall['pearson_vs_saluki_reference']:.4f}。",
        f"- 在其他派生标签里，最能恢复 `saluki_human_pc1` 的方案是 "
        f"`{best_non_saluki['target']} + {best_non_saluki['setting']}`，"
        f"Pearson={best_non_saluki['pearson_vs_saluki_reference']:.4f}。",
        "",
        "## 如何解读",
        "",
        "- 如果最优组合仍然是 `saluki_human_pc1`，说明 `mouse prior` 没有改变“Saluki 发布标签最适合作为主监督目标”这个结论。",
        "- 如果某个 `MOESM3` 派生标签只在 `prior-enhanced` 设定下略优，它更适合被描述为补充探索，而不是替换主线标签。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_moesm3_saluki_label_prior_benchmark(
    feature_path: Path | None = None,
    results_root: Path | None = None,
    targets: list[str] | None = None,
) -> tuple[Path, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "moesm3_saluki_label_prior_benchmark"
    results_root.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    target_names = targets if targets is not None else DEFAULT_TARGETS
    target_tables = _load_target_tables(target_names)
    saluki_reference = _read_target(TARGET_SPECS["saluki_human_pc1"], "saluki_reference_pc1")
    feature_df, base_cols = _load_global_feature_table(feature_path)
    settings = _build_global_settings(base_cols)

    common_gene_ids: set[str] = set(feature_df["gene_id"].astype(str))
    common_gene_ids &= set(saluki_reference["gene_id"].astype(str))
    for target_df in target_tables.values():
        common_gene_ids &= set(target_df["gene_id"].astype(str))
    common_gene_ids = set(sorted(common_gene_ids))
    pd.DataFrame({"gene_id": sorted(common_gene_ids)}).to_csv(
        results_root / "common_gene_universe.tsv", sep="\t", index=False
    )

    summary_rows: list[dict[str, object]] = []
    feature_df = feature_df.loc[feature_df["gene_id"].isin(common_gene_ids)].copy()
    feature_df = feature_df.sort_values("gene_id").reset_index(drop=True)
    for target_name, target_df in target_tables.items():
        merged = (
            feature_df.merge(target_df, on="gene_id", how="inner")
            .merge(saluki_reference, on="gene_id", how="inner")
            .sort_values("gene_id")
            .reset_index(drop=True)
        )
        target_values = merged[target_name].to_numpy(dtype=np.float32)
        saluki_values = merged["saluki_reference_pc1"].to_numpy(dtype=np.float32)
        label_corr = _correlation(target_values, saluki_values)
        for setting_name, cols in settings.items():
            pred, best_iterations = _xgb_oof_predictions(merged, cols, target_name)
            target_metrics = _metrics(target_values, pred)
            saluki_metrics = _metrics(saluki_values, pred)

            oof = merged[["gene_id", target_name, "saluki_reference_pc1"]].copy()
            oof["prediction"] = pred
            oof.to_csv(oof_dir / f"{target_name}__{setting_name}.tsv", sep="\t", index=False)

            summary_rows.append(
                {
                    "target": target_name,
                    "setting": setting_name,
                    "n": int(len(merged)),
                    "n_features": int(len(cols)),
                    "n_with_mouse_pc1": int(merged["has_mouse_pc1"].sum()),
                    "n_with_saluki_mouse_prior": int(merged["has_saluki_mouse_prior"].sum()),
                    "target_vs_saluki_pearson": label_corr["pearson"],
                    "target_vs_saluki_spearman": label_corr["spearman"],
                    "pearson_vs_target": target_metrics["pearson"],
                    "spearman_vs_target": target_metrics["spearman"],
                    "r2_vs_target": target_metrics["r2"],
                    "pearson_vs_saluki_reference": saluki_metrics["pearson"],
                    "spearman_vs_saluki_reference": saluki_metrics["spearman"],
                    "r2_vs_saluki_reference": saluki_metrics["r2"],
                    "mean_best_iter": float(np.mean(best_iterations)),
                }
            )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["pearson_vs_saluki_reference", "r2_vs_saluki_reference"], ascending=False
    )
    summary_path = results_root / "summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_rows, handle, indent=2)

    note_path = write_benchmark_note(summary, MANUSCRIPT_DIR.parent / "docs" / "moesm3_saluki_label_prior_benchmark.md")
    return summary_path, note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark MOESM3-derived labels with mouse priors.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR / "moesm3_saluki_label_prior_benchmark",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help="Optional subset of targets to benchmark.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary_path, note_path = run_moesm3_saluki_label_prior_benchmark(
        feature_path=args.feature_path,
        results_root=args.results_root,
        targets=args.targets,
    )
    print(summary_path)
    print(pd.read_csv(summary_path, sep="\t").to_string(index=False))
    print(note_path)


if __name__ == "__main__":
    main()
