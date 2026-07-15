from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _load_global_table
from mrna_half_life_paper.saluki_naming import HAS_SALUKI_MOUSE_PRIOR, SALUKI_HUMAN_PC1


DEFAULT_OOF_ROOT = RESULTS_DIR / "tenfold_global_prior" / "oof_predictions"


def _metrics(y: pd.Series, pred: pd.Series) -> dict[str, float]:
    y_values = y.to_numpy(dtype=float)
    pred_values = pred.to_numpy(dtype=float)
    return {
        "pearson": float(np.corrcoef(y_values, pred_values)[0, 1]),
        "spearman": float(spearmanr(y_values, pred_values).statistic),
        "r2": float(r2_score(y_values, pred_values)),
    }


def _group_masks(df: pd.DataFrame) -> dict[str, pd.Series]:
    has_mouse = df["has_mouse_pc1"].astype(bool)
    has_saluki = df[HAS_SALUKI_MOUSE_PRIOR].astype(bool)
    return {
        "all_genes": pd.Series(True, index=df.index),
        "has_any_mouse_prior": has_mouse | has_saluki,
        "missing_both_mouse_priors": ~(has_mouse | has_saluki),
        "has_both_mouse_priors": has_mouse & has_saluki,
        "missing_at_least_one_mouse_prior": ~(has_mouse & has_saluki),
        "has_reconstructed_mouse_pc1": has_mouse,
        "missing_reconstructed_mouse_pc1": ~has_mouse,
        "has_saluki_mouse_prior": has_saluki,
        "missing_saluki_mouse_prior": ~has_saluki,
        "only_reconstructed_mouse_pc1": has_mouse & ~has_saluki,
        "only_saluki_mouse_prior": ~has_mouse & has_saluki,
    }


def _group_description(group: str) -> str:
    descriptions = {
        "all_genes": "All genes in the global common human gene universe.",
        "has_any_mouse_prior": "Genes with at least one mouse prior value available.",
        "missing_both_mouse_priors": "Genes with neither reconstructed mouse PC1 nor Saluki mouse prior available.",
        "has_both_mouse_priors": "Genes with both reconstructed mouse PC1 and Saluki mouse prior available.",
        "missing_at_least_one_mouse_prior": "Genes missing one or both mouse prior values.",
        "has_reconstructed_mouse_pc1": "Genes with reconstructed mouse PC1 available.",
        "missing_reconstructed_mouse_pc1": "Genes without reconstructed mouse PC1.",
        "has_saluki_mouse_prior": "Genes with Saluki mouse prior available.",
        "missing_saluki_mouse_prior": "Genes without Saluki mouse prior.",
        "only_reconstructed_mouse_pc1": "Genes with reconstructed mouse PC1 but no Saluki mouse prior.",
        "only_saluki_mouse_prior": "Genes with Saluki mouse prior but no reconstructed mouse PC1.",
    }
    return descriptions[group]


def _setting_label(setting: str) -> str:
    labels = {
        "compact_all_global": "human_only",
        "compact_all_plus_both_mouse_priors_global": "real_both_priors",
        "compact_all_plus_both_mouse_priors_shuffled_global": "shuffled_both_priors",
        "coverage_aware_fallback": "coverage_aware_fallback",
    }
    return labels.get(setting, setting)


def _read_oof(setting: str, seed: int, cv_splits: int, oof_root: Path) -> pd.DataFrame:
    path = oof_root / f"{setting}__folds{cv_splits}__seed{seed}.tsv"
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    return frame.rename(columns={"prediction": f"prediction__{setting}"})


def _aggregate(rows: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["pearson", "spearman", "r2", "delta_vs_human_only_pearson", "delta_vs_human_only_r2"]
    out_rows: list[dict[str, float | int | str]] = []
    group_cols = ["coverage_group", "coverage_description", "setting", "setting_label"]
    for keys, group in rows.groupby(group_cols, sort=False):
        coverage_group, coverage_description, setting, setting_label = keys
        out: dict[str, float | int | str] = {
            "coverage_group": str(coverage_group),
            "coverage_description": str(coverage_description),
            "setting": str(setting),
            "setting_label": str(setting_label),
            "n_random_states": int(group["random_state"].nunique()),
            "cv_splits": int(group["cv_splits"].iloc[0]),
            "n": int(group["n"].iloc[0]),
            "n_with_mouse_pc1": int(group["n_with_mouse_pc1"].iloc[0]),
            "n_with_saluki_mouse_prior": int(group["n_with_saluki_mouse_prior"].iloc[0]),
        }
        for metric in metric_cols:
            out[f"{metric}_mean"] = float(group[metric].mean())
            out[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if group.shape[0] > 1 else 0.0
        out_rows.append(out)
    return pd.DataFrame(out_rows)


def _write_note(summary: pd.DataFrame, note_path: Path) -> Path:
    def row(group: str, label: str) -> pd.Series:
        return summary.loc[(summary["coverage_group"] == group) & (summary["setting_label"] == label)].iloc[0]

    all_real = row("all_genes", "real_both_priors")
    all_fallback = row("all_genes", "coverage_aware_fallback")
    no_prior_real = row("missing_both_mouse_priors", "real_both_priors")
    no_prior_human = row("missing_both_mouse_priors", "human_only")
    both_real = row("has_both_mouse_priors", "real_both_priors")
    both_human = row("has_both_mouse_priors", "human_only")
    lines = [
        "# Mouse-Prior Missingness Analysis",
        "",
        "## Purpose",
        "",
        "评估 global prior-enhanced 模型在测试基因缺少 mouse prior 时的真实 OOF 表现。",
        "该分析不重新训练模型，而是复用 10-fold multi-seed OOF predictions，并按 mouse-prior 覆盖状态分层计算指标。",
        "",
        "## Key Results",
        "",
        f"- 全部 12,916 个基因上，real both-prior 模型 Pearson={all_real['pearson_mean']:.4f}±{all_real['pearson_sd']:.4f}。",
        f"- coverage-aware fallback（有 mouse prior 用 prior-enhanced；两种 prior 都缺失时退回 human-only）"
        f"Pearson={all_fallback['pearson_mean']:.4f}±{all_fallback['pearson_sd']:.4f}。",
        f"- 两种 mouse prior 都缺失的 {int(no_prior_real['n'])} 个基因上，real both-prior 模型 Pearson={no_prior_real['pearson_mean']:.4f}±{no_prior_real['pearson_sd']:.4f}，"
        f"human-only 模型 Pearson={no_prior_human['pearson_mean']:.4f}±{no_prior_human['pearson_sd']:.4f}。",
        f"- 两种 mouse prior 都可用的 {int(both_real['n'])} 个基因上，real both-prior 模型 Pearson={both_real['pearson_mean']:.4f}±{both_real['pearson_sd']:.4f}，"
        f"human-only 模型 Pearson={both_human['pearson_mean']:.4f}±{both_human['pearson_sd']:.4f}。",
        "",
        "## Interpretation",
        "",
        "- 缺少 mouse prior 的基因仍可预测，但不能期待 prior-enhanced 主结果的性能水平。",
        "- 实际部署时可使用 coverage-aware fallback，而不是强行让 missing-prior genes 走 prior-enhanced 模型。",
        "- missing_both_mouse_priors 子集样本量较小，指标方差和选择偏差需要谨慎解释。",
        "- 该结果支持在论文中把模型应用范围写清楚：没有 mouse prior 时可退回 human sequence/regulatory 预测；最高性能依赖可用 ortholog prior。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_prior_missingness_analysis(
    feature_path: Path | None = None,
    oof_root: Path | None = None,
    results_root: Path | None = None,
    random_states: list[int] | None = None,
    cv_splits: int = 10,
) -> dict[str, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    oof_root = oof_root if oof_root is not None else DEFAULT_OOF_ROOT
    results_root = results_root if results_root is not None else RESULTS_DIR / "prior_missingness_analysis"
    random_states = random_states if random_states is not None else [13, 42, 101]
    settings = [
        "compact_all_global",
        "compact_all_plus_both_mouse_priors_global",
        "compact_all_plus_both_mouse_priors_shuffled_global",
    ]
    results_root.mkdir(parents=True, exist_ok=True)

    df, _ = _load_global_table(feature_path)
    availability = df[["gene_id", SALUKI_HUMAN_PC1, "has_mouse_pc1", HAS_SALUKI_MOUSE_PRIOR]].copy()
    masks = _group_masks(availability)

    rows: list[dict[str, float | int | str]] = []
    for seed in random_states:
        seed_predictions = {
            setting: _read_oof(setting, seed, cv_splits, oof_root)
            for setting in settings
        }
        merged = availability.copy()
        for setting, oof in seed_predictions.items():
            merged = merged.merge(oof[["gene_id", f"prediction__{setting}"]], on="gene_id", how="inner")

        human_col = "prediction__compact_all_global"
        prior_col = "prediction__compact_all_plus_both_mouse_priors_global"
        has_any_prior = merged["has_mouse_pc1"].astype(bool) | merged[HAS_SALUKI_MOUSE_PRIOR].astype(bool)
        merged["prediction__coverage_aware_fallback"] = np.where(
            has_any_prior,
            merged[prior_col],
            merged[human_col],
        )
        all_settings = settings + ["coverage_aware_fallback"]
        for coverage_group, mask in masks.items():
            sub = merged.loc[mask.reindex(merged.index).fillna(False)].copy()
            if sub.shape[0] < 3:
                continue
            human_metrics = _metrics(sub[SALUKI_HUMAN_PC1], sub[human_col])
            for setting in all_settings:
                pred_col = f"prediction__{setting}"
                metrics = _metrics(sub[SALUKI_HUMAN_PC1], sub[pred_col])
                rows.append(
                    {
                        "coverage_group": coverage_group,
                        "coverage_description": _group_description(coverage_group),
                        "setting": setting,
                        "setting_label": _setting_label(setting),
                        "random_state": int(seed),
                        "cv_splits": int(cv_splits),
                        "n": int(sub.shape[0]),
                        "n_with_mouse_pc1": int(sub["has_mouse_pc1"].sum()),
                        "n_with_saluki_mouse_prior": int(sub[HAS_SALUKI_MOUSE_PRIOR].sum()),
                        **metrics,
                        "delta_vs_human_only_pearson": float(metrics["pearson"] - human_metrics["pearson"]),
                        "delta_vs_human_only_r2": float(metrics["r2"] - human_metrics["r2"]),
                    }
                )

    by_seed = pd.DataFrame(rows)
    by_seed_path = results_root / "summary_by_seed.tsv"
    by_seed.to_csv(by_seed_path, sep="\t", index=False)
    summary = _aggregate(by_seed)
    sort_order = {
        "all_genes": 0,
        "has_any_mouse_prior": 1,
        "missing_both_mouse_priors": 2,
        "has_both_mouse_priors": 3,
        "missing_at_least_one_mouse_prior": 4,
        "has_reconstructed_mouse_pc1": 5,
        "missing_reconstructed_mouse_pc1": 6,
        "has_saluki_mouse_prior": 7,
        "missing_saluki_mouse_prior": 8,
        "only_reconstructed_mouse_pc1": 9,
        "only_saluki_mouse_prior": 10,
    }
    setting_order = {
        "real_both_priors": 0,
        "coverage_aware_fallback": 1,
        "human_only": 2,
        "shuffled_both_priors": 3,
    }
    summary = summary.sort_values(
        by=["coverage_group", "setting_label"],
        key=lambda s: s.map(sort_order).fillna(99) if s.name == "coverage_group" else s.map(setting_order).fillna(99),
    ).reset_index(drop=True)
    summary_path = results_root / "summary_by_coverage.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    note_path = _write_note(summary, MANUSCRIPT_DIR.parent / "docs" / "prior_missingness_analysis.md")
    return {
        "summary_by_seed": by_seed_path,
        "summary_by_coverage": summary_path,
        "note": note_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate global prior models by mouse-prior availability.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument("--oof-root", type=Path, default=DEFAULT_OOF_ROOT)
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "prior_missingness_analysis")
    parser.add_argument("--random-states", nargs="+", type=int, default=[13, 42, 101])
    parser.add_argument("--cv-splits", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_prior_missingness_analysis(
        feature_path=args.feature_path,
        oof_root=args.oof_root,
        results_root=args.results_root,
        random_states=args.random_states,
        cv_splits=args.cv_splits,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
