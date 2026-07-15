from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import FIGURES_DIR, MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.label_competition import _load_prior_feature_table, _xgb_oof_predictions
from mrna_half_life_paper.ortholog_regularized_label import (
    _display_label_name,
    _load_saluki_label,
    _load_human_base_label,
    _load_reconstructed_mouse_pc1,
    _mean_sd_text,
    _ortholog_correlation,
    make_ortholog_regularized_label,
)
from mrna_half_life_paper.saluki_naming import HAS_SALUKI_MOUSE_PRIOR, SALUKI_MOUSE_HIGH_CONFIDENCE, SALUKI_MOUSE_PRIOR


SALUKI_PRIOR_BASELINE = 0.836801
SALUKI_PURE_BASELINE = 0.749565


@dataclass(frozen=True)
class FeatureSetting:
    name: str
    description: str
    columns: list[str]


def _source_label(source: str) -> pd.Series:
    if source == "reconstructed_mouse_pc1":
        return _load_reconstructed_mouse_pc1()
    if source == "saluki_mouse_prior":
        return _load_saluki_label("mouse")
    raise ValueError(f"Unknown source: {source}")


def _feature_settings(base_cols: list[str], source: str) -> list[FeatureSetting]:
    mouse_triple = ["mouse_pc1", "has_mouse_pc1", "mouse_pc1_highconf"]
    saluki_mouse_triple = [SALUKI_MOUSE_PRIOR, HAS_SALUKI_MOUSE_PRIOR, SALUKI_MOUSE_HIGH_CONFIDENCE]
    all_masks = [
        "has_mouse_pc1",
        "mouse_pc1_highconf",
        HAS_SALUKI_MOUSE_PRIOR,
        SALUKI_MOUSE_HIGH_CONFIDENCE,
    ]
    if source == "reconstructed_mouse_pc1":
        direct_triple = mouse_triple
        alternate_triple = saluki_mouse_triple
        direct_name = "direct_reconstructed_mouse_pc1"
        alternate_name = "alternate_saluki_mouse_prior"
        no_direct_name = "no_direct_reconstructed_mouse_pc1"
    elif source == "saluki_mouse_prior":
        direct_triple = saluki_mouse_triple
        alternate_triple = mouse_triple
        direct_name = "direct_saluki_mouse_prior"
        alternate_name = "alternate_reconstructed_mouse_pc1"
        no_direct_name = "no_direct_saluki_mouse_prior"
    else:
        raise ValueError(f"Unknown source: {source}")

    return [
        FeatureSetting(
            name="human_only",
            description="Compact human sequence and regulatory features only.",
            columns=base_cols,
        ),
        FeatureSetting(
            name="mask_only",
            description="Human features plus mouse-prior availability/confidence indicators, but no mouse values.",
            columns=base_cols + all_masks,
        ),
        FeatureSetting(
            name=alternate_name,
            description="Human features plus the non-matching mouse prior source.",
            columns=base_cols + alternate_triple,
        ),
        FeatureSetting(
            name=no_direct_name,
            description="Human features plus all indicators and only the non-matching mouse prior value.",
            columns=base_cols + all_masks + [alternate_triple[0]],
        ),
        FeatureSetting(
            name=direct_name,
            description="Human features plus the same mouse prior source used in label regularization.",
            columns=base_cols + direct_triple,
        ),
        FeatureSetting(
            name="both_mouse_priors",
            description="Human features plus both mouse prior values and their indicators.",
            columns=base_cols + mouse_triple + saluki_mouse_triple,
        ),
    ]


def _display_setting_name(name: str) -> str:
    return {
        "human_only": "human_only",
        "mask_only": "mask_only",
        "alternate_saluki_mouse_prior": "alternate_saluki_mouse_prior",
        "alternate_reconstructed_mouse_pc1": "alternate_reconstructed_mouse_pc1",
        "no_direct_reconstructed_mouse_pc1": "no_direct_reconstructed_mouse_pc1",
        "no_direct_saluki_mouse_prior": "no_direct_saluki_mouse_prior",
        "direct_reconstructed_mouse_pc1": "direct_reconstructed_mouse_pc1",
        "direct_saluki_mouse_prior": "direct_saluki_mouse_prior",
        "both_mouse_priors": "both_mouse_priors",
    }.get(name, name)


def _aggregate(summary: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [
        "pearson_vs_target",
        "spearman_vs_target",
        "pearson_vs_saluki",
        "spearman_vs_saluki",
        "delta_vs_saluki_prior_baseline",
        "delta_vs_saluki_pure_baseline",
        "mean_best_iter",
    ]
    rows: list[dict[str, float | int | str]] = []
    group_cols = ["label", "source", "lambda", "setting", "setting_description"]
    for keys, group in summary.groupby(group_cols, sort=False):
        label, source, lambda_, setting, description = keys
        row: dict[str, float | int | str] = {
            "label": str(label),
            "source": str(source),
            "lambda": float(lambda_),
            "setting": str(setting),
            "setting_description": str(description),
            "n_random_states": int(group["random_state"].nunique()),
            "cv_splits": int(group["cv_splits"].iloc[0]) if "cv_splits" in group.columns else 5,
            "n": int(group["n"].iloc[0]),
            "n_features": int(group["n_features"].iloc[0]),
            "target_vs_saluki_pearson": float(group["target_vs_saluki_pearson"].iloc[0]),
            "ortholog_pearson": float(group["ortholog_pearson"].iloc[0]),
        }
        for metric in metric_cols:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if group.shape[0] > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True)


def _plot_ablation(aggregate: pd.DataFrame, figures_dir: Path) -> list[Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    out_paths: list[Path] = []
    labels = {
        "human_only": "Human only",
        "mask_only": "Masks only",
        "alternate_saluki_mouse_prior": "Alternate\nSaluki prior",
        "alternate_reconstructed_mouse_pc1": "Alternate\nreconstructed PC1",
        "no_direct_reconstructed_mouse_pc1": "No direct\nreconstructed PC1",
        "no_direct_saluki_mouse_prior": "No direct\nSaluki prior",
        "direct_reconstructed_mouse_pc1": "Direct\nreconstructed PC1",
        "direct_saluki_mouse_prior": "Direct\nSaluki prior",
        "both_mouse_priors": "Both priors",
    }
    colors = {
        "human_only": "#6c757d",
        "mask_only": "#a7b1b7",
        "alternate_saluki_mouse_prior": "#4c78a8",
        "alternate_reconstructed_mouse_pc1": "#4c78a8",
        "no_direct_reconstructed_mouse_pc1": "#72b7b2",
        "no_direct_saluki_mouse_prior": "#72b7b2",
        "direct_reconstructed_mouse_pc1": "#2f7f6f",
        "direct_saluki_mouse_prior": "#2f7f6f",
        "both_mouse_priors": "#1b5e54",
    }
    for label_name, label_df in aggregate.groupby("label", sort=False):
        plot_df = label_df.copy()
        x = np.arange(plot_df.shape[0], dtype=float)
        fig, ax = plt.subplots(figsize=(8.4, 4.2), constrained_layout=True)
        ax.bar(
            x,
            plot_df["pearson_vs_target_mean"],
            yerr=plot_df["pearson_vs_target_sd"],
            color=[colors.get(setting, "#777777") for setting in plot_df["setting"]],
            capsize=3,
            width=0.72,
        )
        ax.axhline(SALUKI_PRIOR_BASELINE, color="black", linestyle="--", linewidth=1.2, label="Saluki human PC1 prior")
        ax.axhline(SALUKI_PURE_BASELINE, color="#555555", linestyle=":", linewidth=1.2, label="Saluki human PC1 pure")
        ax.set_xticks(x)
        ax.set_xticklabels([labels.get(setting, setting) for setting in plot_df["setting"]], fontsize=8)
        ax.set_ylabel("OOF Pearson to target")
        ax.set_title(f"Prior ablation for {_display_label_name(label_name)}", fontsize=11)
        ax.set_ylim(0.70, max(0.90, float(plot_df["pearson_vs_target_mean"].max()) + 0.03))
        ax.grid(axis="y", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        safe_label = label_name.replace(".", "p")
        png = figures_dir / f"{safe_label}_prior_ablation.png"
        pdf = figures_dir / f"{safe_label}_prior_ablation.pdf"
        fig.savefig(png, dpi=300)
        fig.savefig(pdf)
        plt.close(fig)
        out_paths.extend([png, pdf])
    return out_paths


def _write_note(aggregate: pd.DataFrame, note_path: Path, figure_paths: list[Path]) -> Path:
    lines = [
        "# Ortholog-Informed Target-Shrinkage Prior Ablation",
        "",
        "## Purpose",
        "",
        "测试 ortholog-informed shrinkage 标签的预测优势是否只来自模型直接读取同源 mouse prior 值。",
        "",
        "## Results",
        "",
    ]
    for label_name, label_df in aggregate.groupby("label", sort=False):
        label_df = label_df.set_index("setting")
        lines.extend(
            [
                f"### {_display_label_name(label_name)}",
                "",
            ]
        )
        for setting in label_df.index:
            row = label_df.loc[setting]
            lines.append(
                f"- `{_display_setting_name(setting)}`：Pearson(target)="
                f"{_mean_sd_text(row, 'pearson_vs_target_mean', 'pearson_vs_target_sd')}；"
                f"Pearson(Saluki human PC1)="
                f"{_mean_sd_text(row, 'pearson_vs_saluki_mean', 'pearson_vs_saluki_sd')}。"
            )
        if "no_direct_reconstructed_mouse_pc1" in label_df.index and "both_mouse_priors" in label_df.index:
            no_direct = label_df.loc["no_direct_reconstructed_mouse_pc1"]
            both = label_df.loc["both_mouse_priors"]
            lines.append(
                f"- no-direct setting 保留替代 Saluki mouse PC1 后，target Pearson="
                f"{_mean_sd_text(no_direct, 'pearson_vs_target_mean', 'pearson_vs_target_sd')}；"
                f"完整 both-prior setting 为 "
                f"{_mean_sd_text(both, 'pearson_vs_target_mean', 'pearson_vs_target_sd')}。"
            )
            lines.append(
                f"- no-direct setting 相对 Saluki human PC1 prior baseline 的增益为 "
                f"{_mean_sd_text(no_direct, 'delta_vs_saluki_prior_baseline_mean', 'delta_vs_saluki_prior_baseline_sd')}；"
                f"both-prior 相对 no-direct 只再增加 "
                f"{both['pearson_vs_target_mean'] - no_direct['pearson_vs_target_mean']:.4f}。"
            )
        if "no_direct_reconstructed_mouse_pc1" in label_df.index and "direct_reconstructed_mouse_pc1" in label_df.index:
            no_direct = label_df.loc["no_direct_reconstructed_mouse_pc1"]
            direct = label_df.loc["direct_reconstructed_mouse_pc1"]
            lines.append(
                f"- direct reconstructed mouse PC1 相对 no-direct setting 的增益为 "
                f"{direct['pearson_vs_target_mean'] - no_direct['pearson_vs_target_mean']:.4f}。"
            )
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "- 如果 mask-only 接近 human-only，说明 missingness/confidence 指示本身不能解释提升。",
            "- 如果 no-direct setting 低于 both-prior setting，说明同源 mouse prior 值确实贡献了标签可预测性。",
            "- 该结果用于限定论文主张：ortholog-informed shrinkage target 的优势属于 cross-species prior-enhanced setting，而不是 pure human-only sequence setting。",
            "",
            "## Figures",
            "",
        ]
    )
    for path in figure_paths:
        display_path = path
        try:
            display_path = path.relative_to(MANUSCRIPT_DIR.parent)
        except ValueError:
            pass
        lines.append(f"- `{display_path}`")
    lines.append("")
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_prior_ablation(
    lambdas: list[float],
    source: str,
    random_states: list[int],
    feature_path: Path | None = None,
    results_root: Path | None = None,
    device: str | None = None,
    cv_splits: int = 5,
    saluki_prior_baseline: float = SALUKI_PRIOR_BASELINE,
    saluki_pure_baseline: float = SALUKI_PURE_BASELINE,
    make_figures: bool = True,
) -> tuple[Path, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "ortholog_prior_ablation"
    results_root.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    prior_base, base_cols, _ = _load_prior_feature_table(feature_path)
    settings = _feature_settings(base_cols, source)
    human_base = _load_human_base_label()
    mouse_label = _source_label(source)
    human_saluki = _load_saluki_label("human")

    rows: list[dict[str, float | int | str]] = []
    for lambda_ in lambdas:
        label_name = f"orthoreg_{source}_{lambda_:g}"
        label = make_ortholog_regularized_label(
            human_label=human_base,
            mouse_label=mouse_label,
            lambda_=lambda_,
            saluki_reference=human_saluki,
        ).rename(label_name)
        target_df = label.to_frame().reset_index().rename(columns={"index": "gene_id"})
        saluki_df = human_saluki.rename("saluki_reference_pc1").to_frame().reset_index().rename(columns={"index": "gene_id"})
        merged = (
            prior_base.merge(target_df, on="gene_id", how="inner")
            .merge(saluki_df, on="gene_id", how="inner")
            .sort_values("gene_id")
            .reset_index(drop=True)
        )
        target_saluki_corr = label_correlation(label, human_saluki)
        ortholog_corr = _ortholog_correlation(label, mouse_label)
        for setting in settings:
            missing = [column for column in setting.columns if column not in merged.columns]
            if missing:
                raise ValueError(f"Setting {setting.name} contains missing columns: {missing}")
            for random_state in random_states:
                pred, best_iterations, device_used = _xgb_oof_predictions(
                    merged,
                    setting.columns,
                    label_name,
                    random_state=random_state,
                    cv_splits=cv_splits,
                    device=device,
                )
                oof = merged[["gene_id", label_name, "saluki_reference_pc1"]].copy()
                oof["prediction"] = pred
                safe_label = label_name.replace(".", "p")
                oof.to_csv(
                    oof_dir / f"{safe_label}__{setting.name}__folds{cv_splits}__seed{random_state}.tsv",
                    sep="\t",
                    index=False,
                )
                pearson_vs_target = float(merged[label_name].corr(pd.Series(pred)))
                pearson_vs_saluki = float(merged["saluki_reference_pc1"].corr(pd.Series(pred)))
                rows.append(
                    {
                        "label": label_name,
                        "source": source,
                        "lambda": float(lambda_),
                        "setting": setting.name,
                        "setting_description": setting.description,
                        "random_state": int(random_state),
                        "cv_splits": int(cv_splits),
                        "n": int(merged.shape[0]),
                        "n_features": int(len(setting.columns)),
                        "n_with_mouse_pc1": int(merged["has_mouse_pc1"].sum()),
                        "n_with_saluki_mouse_prior": int(merged["has_saluki_mouse_prior"].sum()),
                        "target_vs_saluki_pearson": target_saluki_corr["pearson"],
                        "target_vs_saluki_spearman": target_saluki_corr["spearman"],
                        "ortholog_n": int(ortholog_corr["n"]),
                        "ortholog_pearson": ortholog_corr["pearson"],
                        "ortholog_spearman": ortholog_corr["spearman"],
                        "pearson_vs_target": pearson_vs_target,
                        "spearman_vs_target": float(merged[label_name].corr(pd.Series(pred), method="spearman")),
                        "pearson_vs_saluki": pearson_vs_saluki,
                        "spearman_vs_saluki": float(
                            merged["saluki_reference_pc1"].corr(pd.Series(pred), method="spearman")
                        ),
                        "delta_vs_saluki_prior_baseline": pearson_vs_target - saluki_prior_baseline,
                        "delta_vs_saluki_pure_baseline": pearson_vs_target - saluki_pure_baseline,
                        "mean_best_iter": float(np.mean(best_iterations)),
                        "device_used": device_used,
                    }
                )

    summary = pd.DataFrame(rows)
    summary_path = results_root / "summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    aggregate = _aggregate(summary)
    aggregate_path = results_root / "summary_by_setting.tsv"
    aggregate.to_csv(aggregate_path, sep="\t", index=False)
    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)
    figure_paths = (
        _plot_ablation(aggregate, FIGURES_DIR / "ortholog_regularized_label" / "prior_ablation")
        if make_figures
        else []
    )
    note_path = _write_note(aggregate, MANUSCRIPT_DIR.parent / "docs" / "ortholog_prior_ablation.md", figure_paths)
    return summary_path, note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run prior ablations for ortholog-informed shrinkage targets.")
    parser.add_argument("--lambdas", nargs="+", type=float, default=[0.1])
    parser.add_argument("--source", default="reconstructed_mouse_pc1")
    parser.add_argument("--random-states", nargs="+", type=int, default=[13, 42, 101])
    parser.add_argument("--cv-splits", type=int, default=5)
    parser.add_argument("--saluki-prior-baseline", type=float, default=SALUKI_PRIOR_BASELINE)
    parser.add_argument("--saluki-pure-baseline", type=float, default=SALUKI_PURE_BASELINE)
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "ortholog_prior_ablation")
    parser.add_argument("--device", default=None)
    parser.add_argument("--no-figures", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary_path, note_path = run_prior_ablation(
        lambdas=args.lambdas,
        source=args.source,
        random_states=args.random_states,
        feature_path=args.feature_path,
        results_root=args.results_root,
        device=args.device,
        cv_splits=args.cv_splits,
        saluki_prior_baseline=args.saluki_prior_baseline,
        saluki_pure_baseline=args.saluki_pure_baseline,
        make_figures=not args.no_figures,
    )
    print(summary_path)
    print(pd.read_csv(summary_path, sep="\t").to_string(index=False))
    aggregate_path = args.results_root / "summary_by_setting.tsv"
    print(aggregate_path)
    print(pd.read_csv(aggregate_path, sep="\t").to_string(index=False))
    print(note_path)


if __name__ == "__main__":
    main()
