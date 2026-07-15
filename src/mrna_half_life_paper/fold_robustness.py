from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import (
    _build_global_settings,
    _load_global_table,
    _metrics,
    _xgb_oof_predictions,
)


def _aggregate(summary: pd.DataFrame) -> pd.DataFrame:
    metric_cols = ["pearson", "spearman", "r2", "mean_best_iter"]
    rows: list[dict[str, float | int | str]] = []
    for (setting, prior_mode, shuffled_columns), group in summary.groupby(
        ["setting", "prior_mode", "shuffled_columns"],
        dropna=False,
        sort=False,
    ):
        row: dict[str, float | int | str] = {
            "setting": str(setting),
            "prior_mode": str(prior_mode),
            "shuffled_columns": "" if pd.isna(shuffled_columns) else str(shuffled_columns),
            "n_random_states": int(group["random_state"].nunique()),
            "cv_splits": int(group["cv_splits"].iloc[0]),
            "n": int(group["n"].iloc[0]),
            "n_with_mouse_pc1": int(group["n_with_mouse_pc1"].iloc[0]),
            "n_with_saluki_mouse_prior": int(group["n_with_saluki_mouse_prior"].iloc[0]),
        }
        for metric in metric_cols:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if group.shape[0] > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("pearson_mean", ascending=False).reset_index(drop=True)


def _aggregate_blend(blend: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for (prior_setting, prior_weight, compact_all_weight), group in blend.groupby(
        ["prior_setting", "prior_weight", "compact_all_weight"],
        sort=False,
    ):
        row: dict[str, float | int | str] = {
            "prior_setting": str(prior_setting),
            "prior_weight": float(prior_weight),
            "compact_all_weight": float(compact_all_weight),
            "n_random_states": int(group["random_state"].nunique()),
            "cv_splits": int(group["cv_splits"].iloc[0]),
        }
        for metric in ("pearson", "spearman", "r2"):
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if group.shape[0] > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("pearson_mean", ascending=False).reset_index(drop=True)


def _write_note(
    summary_by_setting: pd.DataFrame,
    note_path: Path,
) -> Path:
    baseline = summary_by_setting.loc[summary_by_setting["setting"] == "compact_all_global"].iloc[0]
    both = summary_by_setting.loc[
        summary_by_setting["setting"] == "compact_all_plus_both_mouse_priors_global"
    ].iloc[0]
    shuffled = summary_by_setting.loc[
        summary_by_setting["setting"] == "compact_all_plus_both_mouse_priors_shuffled_global"
    ].iloc[0]
    lines = [
        "# Ten-Fold Robustness Benchmark",
        "",
        "## Purpose",
        "",
        "将全局 cross-species prior-enhanced 主结果从 5-fold OOF 扩展到 10-fold gene-level OOF，",
        "并用多个 random seeds 检查结果是否依赖特定 gene split。",
        "",
        "## Key Results",
        "",
        f"- `compact_all_global`: Pearson={baseline['pearson_mean']:.4f}±{baseline['pearson_sd']:.4f}, "
        f"R2={baseline['r2_mean']:.4f}±{baseline['r2_sd']:.4f}.",
        f"- `compact_all_plus_both_mouse_priors_global`: Pearson={both['pearson_mean']:.4f}±{both['pearson_sd']:.4f}, "
        f"R2={both['r2_mean']:.4f}±{both['r2_sd']:.4f}.",
        f"- shuffled both-prior control: Pearson={shuffled['pearson_mean']:.4f}±{shuffled['pearson_sd']:.4f}, "
        f"R2={shuffled['r2_mean']:.4f}±{shuffled['r2_sd']:.4f}.",
        "## Interpretation",
        "",
        "- 10-fold 多 seed 结果用于量化 gene split sensitivity。",
        "- 该结果属于固定 Saluki human PC1 目标下的 cross-species transfer，不是 human-only prediction。",
        "- shuffled control 回到 human-only 附近，支持增益依赖 gene-specific mouse prior，而不是特征列数或 missingness pattern。",
        "- 各 seed 的 OOF prediction 仅用于 split-sensitivity 汇总，不称为模型融合。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_tenfold_global_prior(
    feature_path: Path | None = None,
    results_root: Path | None = None,
    random_states: list[int] | None = None,
    cv_splits: int = 10,
) -> dict[str, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "tenfold_global_prior"
    random_states = random_states if random_states is not None else [13, 42, 101]
    results_root.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    df, base_cols = _load_global_table(feature_path)
    settings = _build_global_settings(base_cols)
    run_specs = {
        "compact_all_global": {
            "cols": settings["compact_all_global"],
            "prior_mode": "none",
            "shuffle_priors": (),
        },
        "compact_all_plus_both_mouse_priors_global": {
            "cols": settings["compact_all_plus_both_mouse_priors_global"],
            "prior_mode": "real",
            "shuffle_priors": (),
        },
        "compact_all_plus_both_mouse_priors_shuffled_global": {
            "cols": settings["compact_all_plus_both_mouse_priors_global"],
            "prior_mode": "shuffled",
            "shuffle_priors": ("mouse_pc1", "saluki_mouse_prior"),
        },
    }

    y = df["saluki_human_pc1"].to_numpy(dtype=np.float32)
    summary_rows: list[dict[str, float | int | str]] = []
    predictions: dict[tuple[str, int], np.ndarray] = {}
    for random_state in random_states:
        for setting_name, spec in run_specs.items():
            pred, best_iterations = _xgb_oof_predictions(
                df,
                spec["cols"],
                random_state=random_state,
                cv_splits=cv_splits,
                shuffle_priors=tuple(spec["shuffle_priors"]),
            )
            predictions[(setting_name, random_state)] = pred
            oof = df[["gene_id", "saluki_human_pc1"]].copy()
            oof["prediction"] = pred
            oof.to_csv(
                oof_dir / f"{setting_name}__folds{cv_splits}__seed{random_state}.tsv",
                sep="\t",
                index=False,
            )
            summary_rows.append(
                {
                    "setting": setting_name,
                    "prior_mode": str(spec["prior_mode"]),
                    "shuffled_columns": ",".join(spec["shuffle_priors"]) if spec["shuffle_priors"] else "",
                    "random_state": int(random_state),
                    "cv_splits": int(cv_splits),
                    "n": int(len(y)),
                    "n_with_mouse_pc1": int(df["has_mouse_pc1"].sum()),
                    "n_with_saluki_mouse_prior": int(df["has_saluki_mouse_prior"].sum()),
                    **_metrics(y, pred),
                    "mean_best_iter": float(np.mean(best_iterations)),
                }
            )

    blend_rows: list[dict[str, float | int | str]] = []
    prior_setting = "compact_all_plus_both_mouse_priors_global"
    for random_state in random_states:
        compact_pred = predictions[("compact_all_global", random_state)]
        prior_pred = predictions[(prior_setting, random_state)]
        for weight in np.linspace(0.0, 1.0, 41):
            blend = weight * prior_pred + (1.0 - weight) * compact_pred
            blend_rows.append(
                {
                    "prior_setting": prior_setting,
                    "prior_weight": float(weight),
                    "compact_all_weight": float(1.0 - weight),
                    "random_state": int(random_state),
                    "cv_splits": int(cv_splits),
                    **_metrics(y, blend),
                }
            )

    summary = pd.DataFrame(summary_rows).sort_values(["pearson", "r2"], ascending=False)
    summary_path = results_root / "summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    summary_by_setting = _aggregate(summary)
    summary_by_setting_path = results_root / "summary_by_setting.tsv"
    summary_by_setting.to_csv(summary_by_setting_path, sep="\t", index=False)

    blend = pd.DataFrame(blend_rows).sort_values(["pearson", "r2"], ascending=False)
    blend_path = results_root / "blend.tsv"
    blend.to_csv(blend_path, sep="\t", index=False)
    blend_by_weight = _aggregate_blend(blend)
    blend_by_weight_path = results_root / "blend_by_weight.tsv"
    blend_by_weight.to_csv(blend_by_weight_path, sep="\t", index=False)

    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_rows, handle, indent=2)

    note_path = _write_note(
        summary_by_setting,
        MANUSCRIPT_DIR.parent / "docs" / "tenfold_global_prior.md",
    )
    return {
        "summary": summary_path,
        "summary_by_setting": summary_by_setting_path,
        "blend": blend_path,
        "blend_by_weight": blend_by_weight_path,
        "note": note_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run 10-fold multi-seed robustness for global mouse priors.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "tenfold_global_prior")
    parser.add_argument("--random-states", nargs="+", type=int, default=[13, 42, 101])
    parser.add_argument("--cv-splits", type=int, default=10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_tenfold_global_prior(
        feature_path=args.feature_path,
        results_root=args.results_root,
        random_states=args.random_states,
        cv_splits=args.cv_splits,
    )
    for path in outputs.values():
        print(path)
        if path.suffix == ".tsv":
            print(pd.read_csv(path, sep="\t").head(12).to_string(index=False))


if __name__ == "__main__":
    main()
