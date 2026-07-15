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

from mrna_half_life_paper.config import MANUSCRIPT_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _load_compact_all_features
from mrna_half_life_paper.saluki_naming import saluki_label_filename

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None


TARGET_SPECS = {
    "saluki_human_pc1": PROCESSED_DIR / "real_data" / "human" / saluki_label_filename("human"),
    "moesm3_all_samples_as_is_pc1": RESULTS_DIR / "moesm3_saluki_matrix_check" / "human" / "all_samples" / "as_is_pc1.tsv",
    "moesm3_all_samples_zscore_only_pc1": RESULTS_DIR
    / "moesm3_saluki_matrix_check"
    / "human"
    / "all_samples"
    / "zscore_only_pc1.tsv",
    "moesm3_all_samples_rerun_pipeline_pc1": RESULTS_DIR
    / "moesm3_saluki_matrix_check"
    / "human"
    / "all_samples"
    / "rerun_pipeline_pc1.tsv",
    "moesm3_no_gejman_as_is_pc1": RESULTS_DIR
    / "moesm3_saluki_matrix_check"
    / "human"
    / "no_gejman"
    / "as_is_pc1.tsv",
    "moesm3_no_gejman_zscore_only_pc1": RESULTS_DIR
    / "moesm3_saluki_matrix_check"
    / "human"
    / "no_gejman"
    / "zscore_only_pc1.tsv",
    "moesm3_no_gejman_rerun_pipeline_pc1": RESULTS_DIR
    / "moesm3_saluki_matrix_check"
    / "human"
    / "no_gejman"
    / "rerun_pipeline_pc1.tsv",
}


def _read_target(path: Path, column_name: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    value_col = "saluki_human_pc1" if "saluki_human_pc1" in frame.columns else "pc1"
    return frame[["gene_id", value_col]].rename(columns={value_col: column_name})


def _load_targets(target_names: list[str] | None = None) -> dict[str, pd.DataFrame]:
    selected = target_names if target_names is not None else list(TARGET_SPECS.keys())
    targets: dict[str, pd.DataFrame] = {}
    for target_name in selected:
        if target_name not in TARGET_SPECS:
            raise ValueError(f"Unknown target: {target_name}")
        targets[target_name] = _read_target(TARGET_SPECS[target_name], target_name)
    return targets


def _load_feature_blocks(feature_path: Path) -> dict[str, pd.DataFrame]:
    base_sequence = pd.read_csv(feature_path, sep="\t").drop(columns=["transcript_id"], errors="ignore")
    compact_all = _load_compact_all_features(feature_path)
    return {
        "base_sequence": base_sequence,
        "compact_all": compact_all,
    }


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


def _xgb_oof_predictions(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    random_state: int = 42,
    cv_splits: int = 5,
) -> tuple[np.ndarray, list[int]]:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for the MOESM3 Saluki label benchmark")

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
    compact = summary.loc[summary["feature_set"] == "compact_all"].copy()
    base = summary.loc[summary["feature_set"] == "base_sequence"].copy()
    compact_saluki = compact.loc[compact["target"] == "saluki_human_pc1"].iloc[0]
    base_saluki = base.loc[base["target"] == "saluki_human_pc1"].iloc[0]
    best_saluki_recovery = summary.sort_values(
        ["pearson_vs_saluki_reference", "r2_vs_saluki_reference"], ascending=False
    ).iloc[0]
    best_self = summary.sort_values(["pearson_vs_target", "r2_vs_target"], ascending=False).iloc[0]
    best_non_saluki = summary.loc[summary["target"] != "saluki_human_pc1"].sort_values(
        ["pearson_vs_saluki_reference", "r2_vs_saluki_reference"], ascending=False
    ).iloc[0]

    lines = [
        "# MOESM3 Saluki Label Downstream Benchmark",
        "",
        "更新时间：2026-05-16",
        "",
        "## 目的",
        "",
        "在不覆盖现有 `saluki_human_pc1` 主线 benchmark 的前提下，单独测试：",
        "",
        "1. 如果把 `MOESM3` Saluki 处理矩阵重新算出来的替代标签当训练目标，下游模型会不会更容易学；",
        "2. 这些替代标签训练出来的 OOF 预测，最终能不能更好恢复 Saluki 发布的 `saluki_human_pc1`。",
        "",
        "## 评估口径",
        "",
        "- 共同基因宇宙：所有 MOESM3 衍生标签与两套 human 特征的交集。",
        "- 特征块：`base_sequence`（526 维）和 `compact_all`（1802 维）。",
        "- 训练器：与主线一致的 `XGBoost/CUDA` 五折 OOF。",
        "- 同时报告两种相关：",
        "  - `pearson_vs_target`：对训练标签本身的拟合能力。",
        "  - `pearson_vs_saluki_reference`：对 Saluki 发布的 `saluki_human_pc1` 的恢复能力。",
        "",
        "## 关键结果",
        "",
        f"- `compact_all + saluki_human_pc1` 基线：Pearson(target)={compact_saluki['pearson_vs_target']:.4f}，"
        f"Pearson(Saluki ref)={compact_saluki['pearson_vs_saluki_reference']:.4f}。",
        f"- `base_sequence + saluki_human_pc1` 基线：Pearson(target)={base_saluki['pearson_vs_target']:.4f}，"
        f"Pearson(Saluki ref)={base_saluki['pearson_vs_saluki_reference']:.4f}。",
        f"- 全部组合里，对 Saluki `saluki_human_pc1` 恢复最好的方案是 "
        f"`{best_saluki_recovery['feature_set']} + {best_saluki_recovery['target']}`，"
        f"Pearson={best_saluki_recovery['pearson_vs_saluki_reference']:.4f}。",
        f"- 全部组合里，对训练标签本身最好学的方案是 "
        f"`{best_self['feature_set']} + {best_self['target']}`，"
        f"Pearson={best_self['pearson_vs_target']:.4f}。",
        f"- 在其他派生标签里，最能恢复 `saluki_human_pc1` 的方案是 "
        f"`{best_non_saluki['feature_set']} + {best_non_saluki['target']}`，"
        f"Pearson={best_non_saluki['pearson_vs_saluki_reference']:.4f}。",
        "",
        "## 如何解读",
        "",
        "- 如果某个替代标签的 `pearson_vs_target` 更高，但 `pearson_vs_saluki_reference` 没超过 `saluki_human_pc1` 基线，说明它更容易学，但不一定更适合作为主任务标签。",
        "- 如果某个替代标签同时提升了 `pearson_vs_target` 和 `pearson_vs_saluki_reference`，仍需独立检验其与 Saluki human PC1 的 estimand 差异。",
        "- 这组结果作为补充控制，用于区分标签处理带来的结构变化与单纯改变目标定义。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_moesm3_saluki_label_benchmark(
    feature_path: Path | None = None,
    results_root: Path | None = None,
    targets: list[str] | None = None,
) -> tuple[Path, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "moesm3_saluki_label_benchmark"
    results_root.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    target_tables = _load_targets(targets)
    saluki_reference = _read_target(TARGET_SPECS["saluki_human_pc1"], "saluki_reference_pc1")
    feature_blocks = _load_feature_blocks(feature_path)

    common_gene_ids: set[str] | None = None
    for feature_df in feature_blocks.values():
        gene_ids = set(feature_df["gene_id"].astype(str))
        common_gene_ids = gene_ids if common_gene_ids is None else common_gene_ids & gene_ids
    for target_df in target_tables.values():
        common_gene_ids &= set(target_df["gene_id"].astype(str))
    common_gene_ids &= set(saluki_reference["gene_id"].astype(str))
    if common_gene_ids is None:
        raise RuntimeError("Failed to construct common gene universe")

    common_gene_ids = set(sorted(common_gene_ids))
    pd.DataFrame({"gene_id": sorted(common_gene_ids)}).to_csv(
        results_root / "common_gene_universe.tsv", sep="\t", index=False
    )

    summary_rows: list[dict[str, object]] = []
    for feature_set, feature_df in feature_blocks.items():
        feature_df = feature_df.loc[feature_df["gene_id"].isin(common_gene_ids)].copy()
        feature_df = feature_df.sort_values("gene_id").reset_index(drop=True)
        for target_name, target_df in target_tables.items():
            merged = (
                feature_df.merge(target_df, on="gene_id", how="inner")
                .merge(saluki_reference, on="gene_id", how="inner")
                .sort_values("gene_id")
                .reset_index(drop=True)
            )
            feature_cols = [
                col
                for col in merged.columns
                if col not in {"gene_id", target_name, "saluki_reference_pc1"}
            ]
            pred, best_iterations = _xgb_oof_predictions(merged, feature_cols, target_name)
            target_values = merged[target_name].to_numpy(dtype=np.float32)
            saluki_values = merged["saluki_reference_pc1"].to_numpy(dtype=np.float32)
            target_metrics = _metrics(target_values, pred)
            saluki_metrics = _metrics(saluki_values, pred)
            label_corr = _correlation(target_values, saluki_values)

            oof = merged[["gene_id", target_name, "saluki_reference_pc1"]].copy()
            oof["prediction"] = pred
            oof.to_csv(oof_dir / f"{feature_set}__{target_name}.tsv", sep="\t", index=False)

            summary_rows.append(
                {
                    "feature_set": feature_set,
                    "target": target_name,
                    "n": int(len(merged)),
                    "n_features": int(len(feature_cols)),
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
        ["pearson_vs_saluki_reference", "pearson_vs_target"], ascending=False
    )
    summary_path = results_root / "summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_rows, handle, indent=2)

    note_path = write_benchmark_note(summary, MANUSCRIPT_DIR.parent / "docs" / "moesm3_saluki_label_benchmark.md")
    return summary_path, note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark MOESM3-derived Saluki-matrix labels with human feature blocks.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
        help="Base human sequence feature table used to build both base_sequence and compact_all blocks.",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR / "moesm3_saluki_label_benchmark",
        help="Output directory. Defaults to results/moesm3_saluki_label_benchmark.",
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
    summary_path, note_path = run_moesm3_saluki_label_benchmark(
        feature_path=args.feature_path,
        results_root=args.results_root,
        targets=args.targets,
    )
    print(summary_path)
    print(pd.read_csv(summary_path, sep="\t").to_string(index=False))
    print(note_path)


if __name__ == "__main__":
    main()
