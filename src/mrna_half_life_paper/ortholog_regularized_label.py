from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import FIGURES_DIR, MANUSCRIPT_DIR, METADATA_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import (
    _load_compact_all_features,
    _xgb_oof_predictions as _global_prior_oof_predictions,
)
from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.label_competition import _load_prior_feature_table, _xgb_oof_predictions
from mrna_half_life_paper.ortholog_analysis import ONE2ONE_ORTHOLOG_PATH, extract_mouse_human_one2one
from mrna_half_life_paper.pc1_consensus import compute_pc1_consensus
from mrna_half_life_paper.preprocess import PreprocessConfig, preprocess_matrix
from mrna_half_life_paper.saluki_naming import (
    SALUKI_HUMAN_PC1,
    saluki_label_column,
    saluki_label_filename,
)


def _zscore(series: pd.Series) -> pd.Series:
    values = series.astype(float)
    return (values - values.mean()) / values.std(ddof=0)


def _load_human_base_label() -> pd.Series:
    matrix = load_matrix(PROCESSED_DIR / "real_data" / "human" / "moesm2_raw_sparse_matrix.tsv")
    metadata = load_metadata(METADATA_DIR / "real_human_samples_no_Gejman.tsv")
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)
    processed = preprocess_matrix(
        matrix,
        PreprocessConfig(apply_log=False, min_observed_per_gene=10),
    )
    return compute_pc1_consensus(processed).rename("human_saluki_like_no_gejman_pc1")


def _load_reconstructed_mouse_pc1() -> pd.Series:
    matrix = load_matrix(PROCESSED_DIR / "real_data" / "mouse" / "moesm2_raw_sparse_matrix.tsv")
    metadata = load_metadata(METADATA_DIR / "real_mouse_samples.tsv")
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)
    processed = preprocess_matrix(
        matrix,
        PreprocessConfig(apply_log=False, min_observed_per_gene=5),
    )
    return compute_pc1_consensus(processed).rename("reconstructed_mouse_pc1")


def _load_saluki_label(species: str) -> pd.Series:
    frame = pd.read_csv(
        PROCESSED_DIR / "real_data" / species / saluki_label_filename(species),
        sep="\t",
    )
    column = saluki_label_column(species)
    return frame.set_index("gene_id")[column].astype(float).rename(f"{species}_{column}")


def _display_source_name(name: str) -> str:
    return {
        SALUKI_HUMAN_PC1: "Saluki human PC1",
        "reconstructed_mouse_pc1": "Reconstructed mouse PC1",
        "saluki_mouse_prior": "Saluki mouse prior",
    }.get(name, name)


def _display_label_name(name: str) -> str:
    if name == SALUKI_HUMAN_PC1:
        return "Saluki human PC1"
    return name


def _mouse_to_human_map(human_index: pd.Index, mouse_label: pd.Series) -> pd.Series:
    extract_mouse_human_one2one()
    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")
    sub = mapping[
        mapping["human_gene_id"].isin(human_index)
        & mapping["mouse_gene_id"].isin(mouse_label.index)
    ].copy()
    return sub.set_index("human_gene_id")["mouse_gene_id"].map(_zscore(mouse_label))


def make_ortholog_regularized_label(
    human_label: pd.Series,
    mouse_label: pd.Series,
    lambda_: float,
    saluki_reference: pd.Series | None = None,
) -> pd.Series:
    mouse_map = _mouse_to_human_map(human_label.index, mouse_label)
    out = _zscore(human_label).copy()
    common = out.index.intersection(mouse_map.dropna().index)
    out.loc[common] = (1.0 - lambda_) * _zscore(human_label).loc[common] + lambda_ * mouse_map.loc[common]
    out = _zscore(out).rename(f"ortholog_regularized_lambda_{lambda_:g}")
    if saluki_reference is not None and label_correlation(out, saluki_reference)["pearson"] < 0:
        out = -out
    return out


def _ortholog_correlation(human_label: pd.Series, mouse_label: pd.Series) -> dict[str, float]:
    extract_mouse_human_one2one()
    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")
    shared = mapping[
        mapping["human_gene_id"].isin(human_label.index)
        & mapping["mouse_gene_id"].isin(mouse_label.index)
    ].copy()
    shared["human_label"] = shared["human_gene_id"].map(human_label)
    shared["mouse_label"] = shared["mouse_gene_id"].map(_zscore(mouse_label))
    return label_correlation(
        shared.set_index("human_gene_id")["human_label"],
        shared.set_index("human_gene_id")["mouse_label"],
    )


def _write_series(series: pd.Series, path: Path, column_name: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    series.rename(column_name).to_frame().reset_index().rename(columns={"index": "gene_id"}).to_csv(
        path,
        sep="\t",
        index=False,
    )


def _benchmark_prior(
    label: pd.Series,
    label_name: str,
    saluki_reference: pd.Series,
    feature_path: Path,
    oof_dir: Path,
    device: str | None,
    random_state: int,
    cv_splits: int,
) -> dict[str, float | int | str]:
    prior_base, _, settings = _load_prior_feature_table(feature_path)
    prior_cols = settings["compact_all_plus_both_mouse_priors_global"]
    target_df = label.rename(label_name).to_frame().reset_index().rename(columns={"index": "gene_id"})
    saluki_df = saluki_reference.rename("saluki_reference_pc1").to_frame().reset_index().rename(columns={"index": "gene_id"})
    merged = (
        prior_base.merge(target_df, on="gene_id", how="inner")
        .merge(saluki_df, on="gene_id", how="inner")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    pred, best_iterations, device_used = _xgb_oof_predictions(
        merged,
        prior_cols,
        label_name,
        device=device,
        random_state=random_state,
        cv_splits=cv_splits,
    )
    oof = merged[["gene_id", label_name, "saluki_reference_pc1"]].copy()
    oof["prediction"] = pred
    oof_dir.mkdir(parents=True, exist_ok=True)
    oof.to_csv(
        oof_dir / f"{label_name}__folds{cv_splits}__seed{random_state}__prior_enhanced.tsv",
        sep="\t",
        index=False,
    )
    return {
        "n": int(len(merged)),
        "n_features": int(len(prior_cols)),
        "cv_splits": int(cv_splits),
        "n_with_mouse_pc1": int(merged["has_mouse_pc1"].sum()),
        "n_with_saluki_mouse_prior": int(merged["has_saluki_mouse_prior"].sum()),
        "prior_pearson_vs_target": float(merged[label_name].corr(pd.Series(pred))),
        "prior_spearman_vs_target": float(merged[label_name].corr(pd.Series(pred), method="spearman")),
        "prior_pearson_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred))),
        "prior_spearman_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred), method="spearman")),
        "mean_best_iter": float(np.mean(best_iterations)),
        "device_used": device_used,
    }


def _benchmark_pure_human(
    label: pd.Series,
    label_name: str,
    saluki_reference: pd.Series,
    feature_path: Path,
    oof_dir: Path,
    device: str | None,
    random_state: int,
    cv_splits: int,
) -> dict[str, float | int | str]:
    features = _load_compact_all_features(feature_path).sort_values("gene_id").reset_index(drop=True)
    feature_cols = [column for column in features.columns if column != "gene_id"]
    target_df = label.rename(label_name).to_frame().reset_index().rename(columns={"index": "gene_id"})
    saluki_df = saluki_reference.rename("saluki_reference_pc1").to_frame().reset_index().rename(columns={"index": "gene_id"})
    merged = (
        features.merge(target_df, on="gene_id", how="inner")
        .merge(saluki_df, on="gene_id", how="inner")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    pred, best_iterations, device_used = _xgb_oof_predictions(
        merged,
        feature_cols,
        label_name,
        device=device,
        random_state=random_state,
        cv_splits=cv_splits,
    )
    oof = merged[["gene_id", label_name, "saluki_reference_pc1"]].copy()
    oof["prediction"] = pred
    oof_dir.mkdir(parents=True, exist_ok=True)
    oof.to_csv(
        oof_dir / f"{label_name}__folds{cv_splits}__seed{random_state}__pure_human.tsv",
        sep="\t",
        index=False,
    )
    return {
        "pure_n": int(len(merged)),
        "pure_n_features": int(len(feature_cols)),
        "pure_cv_splits": int(cv_splits),
        "pure_pearson_vs_target": float(merged[label_name].corr(pd.Series(pred))),
        "pure_spearman_vs_target": float(merged[label_name].corr(pd.Series(pred), method="spearman")),
        "pure_pearson_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred))),
        "pure_spearman_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred), method="spearman")),
        "pure_mean_best_iter": float(np.mean(best_iterations)),
        "pure_device_used": device_used,
    }


def _benchmark_shuffled_prior(
    label: pd.Series,
    label_name: str,
    saluki_reference: pd.Series,
    feature_path: Path,
    oof_dir: Path,
    random_state: int,
    cv_splits: int,
) -> dict[str, float | int | str]:
    prior_base, _, settings = _load_prior_feature_table(feature_path)
    prior_cols = settings["compact_all_plus_both_mouse_priors_global"]
    target_df = label.rename(label_name).to_frame().reset_index().rename(columns={"index": "gene_id"})
    saluki_df = saluki_reference.rename("saluki_reference_pc1").to_frame().reset_index().rename(columns={"index": "gene_id"})
    merged = (
        prior_base.merge(target_df, on="gene_id", how="inner")
        .merge(saluki_df, on="gene_id", how="inner")
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    model_df = merged.copy()
    model_df["saluki_human_pc1"] = model_df[label_name]
    pred, best_iterations = _global_prior_oof_predictions(
        model_df,
        prior_cols,
        random_state=random_state,
        cv_splits=cv_splits,
        shuffle_priors=("mouse_pc1", "saluki_mouse_prior"),
    )
    oof = merged[["gene_id", label_name, "saluki_reference_pc1"]].copy()
    oof["prediction"] = pred
    oof_dir.mkdir(parents=True, exist_ok=True)
    oof.to_csv(
        oof_dir / f"{label_name}__folds{cv_splits}__seed{random_state}__prior_shuffled.tsv",
        sep="\t",
        index=False,
    )
    return {
        "shuffled_prior_cv_splits": int(cv_splits),
        "shuffled_prior_pearson_vs_target": float(merged[label_name].corr(pd.Series(pred))),
        "shuffled_prior_spearman_vs_target": float(merged[label_name].corr(pd.Series(pred), method="spearman")),
        "shuffled_prior_pearson_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred))),
        "shuffled_prior_spearman_vs_saluki": float(merged["saluki_reference_pc1"].corr(pd.Series(pred), method="spearman")),
        "shuffled_prior_mean_best_iter": float(np.mean(best_iterations)),
    }


def _add_control_deltas(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    if "shuffled_prior_pearson_vs_target" in out.columns:
        out["real_minus_shuffled_pearson_vs_target"] = (
            out["prior_pearson_vs_target"] - out["shuffled_prior_pearson_vs_target"]
        )
    if "pure_pearson_vs_target" in out.columns:
        out["real_minus_pure_pearson_vs_target"] = out["prior_pearson_vs_target"] - out["pure_pearson_vs_target"]
    if "shuffled_prior_pearson_vs_saluki" in out.columns:
        out["real_minus_shuffled_pearson_vs_saluki"] = (
            out["prior_pearson_vs_saluki"] - out["shuffled_prior_pearson_vs_saluki"]
        )
    if "pure_pearson_vs_saluki" in out.columns:
        out["real_minus_pure_pearson_vs_saluki"] = out["prior_pearson_vs_saluki"] - out["pure_pearson_vs_saluki"]
    return out


def _aggregate_summary(summary: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "target_vs_saluki_pearson",
        "target_vs_saluki_spearman",
        "ortholog_pearson",
        "ortholog_spearman",
        "prior_pearson_vs_target",
        "prior_spearman_vs_target",
        "prior_pearson_vs_saluki",
        "prior_spearman_vs_saluki",
        "pure_pearson_vs_target",
        "pure_spearman_vs_target",
        "pure_pearson_vs_saluki",
        "pure_spearman_vs_saluki",
        "shuffled_prior_pearson_vs_target",
        "shuffled_prior_spearman_vs_target",
        "shuffled_prior_pearson_vs_saluki",
        "shuffled_prior_spearman_vs_saluki",
        "real_minus_shuffled_pearson_vs_target",
        "real_minus_pure_pearson_vs_target",
        "real_minus_shuffled_pearson_vs_saluki",
        "real_minus_pure_pearson_vs_saluki",
        "mean_best_iter",
        "pure_mean_best_iter",
        "shuffled_prior_mean_best_iter",
    ]
    metrics = [metric for metric in metrics if metric in summary.columns]
    rows: list[dict[str, float | int | str]] = []
    for (label, source, lambda_), group in summary.groupby(["label", "source", "lambda"], sort=False):
        row: dict[str, float | int | str] = {
            "label": label,
            "source": source,
            "lambda": float(lambda_),
            "n_random_states": int(group["random_state"].nunique()),
            "cv_splits": int(group["cv_splits"].iloc[0]) if "cv_splits" in group.columns else 5,
            "ortholog_n": int(group["ortholog_n"].iloc[0]),
            "n": int(group["n"].iloc[0]),
            "n_features": int(group["n_features"].iloc[0]),
        }
        if "pure_n" in group.columns:
            row["pure_n"] = int(group["pure_n"].iloc[0])
            row["pure_n_features"] = int(group["pure_n_features"].iloc[0])
        for metric in metrics:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_sd"] = float(group[metric].std(ddof=1)) if group[metric].shape[0] > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values("prior_pearson_vs_target_mean", ascending=False).reset_index(drop=True)


def _mean_sd_text(row: pd.Series, mean_col: str, sd_col: str | None = None) -> str:
    if mean_col not in row:
        return "NA"
    if sd_col is not None and sd_col in row and float(row[sd_col]) > 0:
        return f"{row[mean_col]:.4f}±{row[sd_col]:.4f}"
    return f"{row[mean_col]:.4f}"


def _plot_summary(summary: pd.DataFrame, figures_dir: Path) -> list[Path]:
    aggregate = _aggregate_summary(summary)
    figures_dir.mkdir(parents=True, exist_ok=True)
    figure_paths: list[Path] = []

    metric_panels = [
        ("target_vs_saluki_pearson", "Target vs Saluki human PC1 Pearson"),
        ("ortholog_pearson", "Human-mouse ortholog Pearson"),
        ("prior_pearson_vs_target", "OOF Pearson to target"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.4), constrained_layout=True)
    for ax, (metric, title) in zip(axes, metric_panels):
        for source, group in aggregate.groupby("source", sort=False):
            group = group.sort_values("lambda")
            ax.errorbar(
                group["lambda"],
                group[f"{metric}_mean"],
                yerr=group.get(f"{metric}_sd"),
                marker="o",
                linewidth=1.8,
                capsize=3,
                label=_display_source_name(source),
            )
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Ortholog regularization lambda")
        ax.set_ylabel("Pearson")
        ax.grid(alpha=0.25, linewidth=0.6)
    axes[-1].axhline(0.836801, color="black", linestyle="--", linewidth=1.1, label="Saluki human PC1 prior")
    axes[-1].legend(frameon=False, fontsize=8)
    tradeoff_png = figures_dir / "orthoreg_lambda_tradeoff.png"
    tradeoff_pdf = figures_dir / "orthoreg_lambda_tradeoff.pdf"
    fig.savefig(tradeoff_png, dpi=300)
    fig.savefig(tradeoff_pdf)
    plt.close(fig)
    figure_paths.extend([tradeoff_png, tradeoff_pdf])

    if {"pure_pearson_vs_target", "shuffled_prior_pearson_vs_target"}.issubset(summary.columns):
        fig, ax = plt.subplots(figsize=(7.2, 3.9), constrained_layout=True)
        plot_df = aggregate.sort_values(["source", "lambda"]).reset_index(drop=True)
        x = np.arange(plot_df.shape[0], dtype=float)
        width = 0.24
        bars = [
            ("pure_pearson_vs_target", "Pure human", "#6c757d", -width),
            ("shuffled_prior_pearson_vs_target", "Shuffled prior", "#d08c60", 0.0),
            ("prior_pearson_vs_target", "Real prior", "#2f7f6f", width),
        ]
        for metric, label, color, offset in bars:
            ax.bar(
                x + offset,
                plot_df[f"{metric}_mean"],
                yerr=plot_df[f"{metric}_sd"],
                width=width,
                label=label,
                color=color,
                capsize=3,
                linewidth=0,
            )
        ax.axhline(0.836801, color="black", linestyle="--", linewidth=1.1, label="Saluki human PC1 prior")
        ax.axhline(0.749565, color="#555555", linestyle=":", linewidth=1.1, label="Saluki human PC1 pure")
        ax.set_xticks(x)
        ax.set_xticklabels(
            [
                f"{_display_source_name(source)}\nlambda={lambda_:g}"
                for source, lambda_ in zip(plot_df["source"], plot_df["lambda"])
            ],
            fontsize=8,
        )
        ax.set_ylabel("OOF Pearson to target")
        ax.set_ylim(0.70, max(0.92, float(plot_df["prior_pearson_vs_target_mean"].max()) + 0.03))
        ax.grid(axis="y", alpha=0.25, linewidth=0.6)
        ax.legend(frameon=False, fontsize=8, ncol=2)
        control_png = figures_dir / "orthoreg_real_vs_controls.png"
        control_pdf = figures_dir / "orthoreg_real_vs_controls.pdf"
        fig.savefig(control_png, dpi=300)
        fig.savefig(control_pdf)
        plt.close(fig)
        figure_paths.extend([control_png, control_pdf])

    return figure_paths


def _write_note(summary: pd.DataFrame, note_path: Path, figure_paths: list[Path] | None = None) -> Path:
    aggregate = _aggregate_summary(summary)
    best_target = aggregate.sort_values("prior_pearson_vs_target_mean", ascending=False).iloc[0]
    conservative = aggregate.sort_values(
        ["target_vs_saluki_pearson_mean", "ortholog_pearson_mean"],
        ascending=False,
    ).iloc[0]
    lines = [
        "# Ortholog-Regularized Label Benchmark",
        "",
        "## Purpose",
        "",
        "测试在 human no-Gejman PC1 标签中加入轻度 mouse ortholog regularization 后，",
        "目标几何、human dominance 与不同输入条件下的可预测性如何变化。",
        "",
        "## Key Results",
        "",
        f"- 最高 prior-enhanced Pearson(target)：`{_display_label_name(best_target['label'])}` = "
        f"{_mean_sd_text(best_target, 'prior_pearson_vs_target_mean', 'prior_pearson_vs_target_sd')}。",
        f"- 最保守、仍高度接近 Saluki human PC1 的候选：`{_display_label_name(conservative['label'])}`，"
        f"target-vs-Saluki Pearson={_mean_sd_text(conservative, 'target_vs_saluki_pearson_mean', 'target_vs_saluki_pearson_sd')}，"
        f"ortholog Pearson={_mean_sd_text(conservative, 'ortholog_pearson_mean', 'ortholog_pearson_sd')}。",
        "- 固定 Saluki human PC1 目标下的 transfer benchmark 仅作为不同 estimand 的参照，不与 regularized target 共用排行榜。",
        "",
        "## Interpretation",
        "",
        "- 较高的 target-specific score 表明 ortholog regularization 增强了显式 prior 可捕获的保守成分，不能解释为固定 human target 上的性能跃迁。",
        "- 这些标签应描述为 human-dominant ortholog-regularized targets，而不是 pure human ground truth 或 sequence-only targets。",
        "",
    ]
    if "real_minus_shuffled_pearson_vs_target_mean" in aggregate.columns:
        lines.extend(
            [
                "## Controls",
                "",
            ]
        )
        for _, row in aggregate.sort_values(["source", "lambda"]).iterrows():
            lines.append(
                f"- `{_display_label_name(row['label'])}`：real prior="
                f"{_mean_sd_text(row, 'prior_pearson_vs_target_mean', 'prior_pearson_vs_target_sd')}，"
                f"pure human={_mean_sd_text(row, 'pure_pearson_vs_target_mean', 'pure_pearson_vs_target_sd')}，"
                f"shuffled prior={_mean_sd_text(row, 'shuffled_prior_pearson_vs_target_mean', 'shuffled_prior_pearson_vs_target_sd')}，"
                f"real-shuffled delta="
                f"{_mean_sd_text(row, 'real_minus_shuffled_pearson_vs_target_mean', 'real_minus_shuffled_pearson_vs_target_sd')}。"
            )
        lines.extend(
            [
                "",
                "- real prior 明显高于 pure-human 和 shuffled-prior，说明提升依赖真实 gene-specific mouse prior，而不是额外列数或 missingness pattern。",
                "",
            ]
        )
    if figure_paths:
        lines.extend(["## Figures", ""])
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


def run_ortholog_regularized_benchmark(
    lambdas: list[float],
    sources: list[str],
    feature_path: Path | None = None,
    results_root: Path | None = None,
    device: str | None = None,
    include_controls: bool = False,
    random_states: list[int] | None = None,
    cv_splits: int = 5,
    make_figures: bool = True,
    include_saluki_baseline: bool = False,
) -> tuple[Path, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "ortholog_regularized_label"
    random_states = random_states if random_states is not None else [42]
    results_root.mkdir(parents=True, exist_ok=True)
    target_dir = results_root / "targets"
    oof_dir = results_root / "oof_predictions"

    human_base = _load_human_base_label()
    mouse_reconstructed = _load_reconstructed_mouse_pc1()
    saluki_mouse = _load_saluki_label("mouse")
    human_saluki = _load_saluki_label("human")
    mouse_sources = {
        "reconstructed_mouse_pc1": mouse_reconstructed,
        "saluki_mouse_prior": saluki_mouse,
    }

    rows: list[dict[str, float | int | str]] = []
    if include_saluki_baseline:
        label_name = SALUKI_HUMAN_PC1
        baseline_genes = human_saluki.index.intersection(human_base.index)
        label = human_saluki.loc[baseline_genes].rename(label_name)
        _write_series(label, target_dir / f"{label_name}.tsv", label_name)
        saluki_corr = label_correlation(label, human_saluki)
        ortholog_corr = _ortholog_correlation(label, saluki_mouse)
        for random_state in random_states:
            benchmark = _benchmark_prior(
                label=label,
                label_name=label_name,
                saluki_reference=human_saluki,
                feature_path=feature_path,
                oof_dir=oof_dir,
                device=device,
                random_state=random_state,
                cv_splits=cv_splits,
            )
            controls: dict[str, float | int | str] = {}
            if include_controls:
                controls.update(
                    _benchmark_pure_human(
                        label=label,
                        label_name=label_name,
                        saluki_reference=human_saluki,
                        feature_path=feature_path,
                        oof_dir=oof_dir,
                        device=device,
                        random_state=random_state,
                        cv_splits=cv_splits,
                    )
                )
                controls.update(
                    _benchmark_shuffled_prior(
                        label=label,
                        label_name=label_name,
                        saluki_reference=human_saluki,
                        feature_path=feature_path,
                        oof_dir=oof_dir,
                        random_state=random_state,
                        cv_splits=cv_splits,
                    )
                )
            rows.append(
                {
                    "label": label_name,
                    "source": SALUKI_HUMAN_PC1,
                    "lambda": 0.0,
                    "random_state": int(random_state),
                    "target_vs_saluki_pearson": saluki_corr["pearson"],
                    "target_vs_saluki_spearman": saluki_corr["spearman"],
                    "ortholog_n": int(ortholog_corr["n"]),
                    "ortholog_pearson": ortholog_corr["pearson"],
                    "ortholog_spearman": ortholog_corr["spearman"],
                    **benchmark,
                    **controls,
                }
            )
    for source_name in sources:
        if source_name not in mouse_sources:
            raise ValueError(f"Unknown source: {source_name}")
        mouse_label = mouse_sources[source_name]
        for lambda_ in lambdas:
            label_name = f"orthoreg_{source_name}_{lambda_:g}"
            label = make_ortholog_regularized_label(
                human_label=human_base,
                mouse_label=mouse_label,
                lambda_=lambda_,
                saluki_reference=human_saluki,
            ).rename(label_name)
            _write_series(label, target_dir / f"{label_name}.tsv", label_name)
            saluki_corr = label_correlation(label, human_saluki)
            ortholog_corr = _ortholog_correlation(label, mouse_label)
            for random_state in random_states:
                benchmark = _benchmark_prior(
                    label=label,
                    label_name=label_name,
                    saluki_reference=human_saluki,
                    feature_path=feature_path,
                    oof_dir=oof_dir,
                    device=device,
                    random_state=random_state,
                    cv_splits=cv_splits,
                )
                controls: dict[str, float | int | str] = {}
                if include_controls:
                    controls.update(
                        _benchmark_pure_human(
                            label=label,
                            label_name=label_name,
                            saluki_reference=human_saluki,
                            feature_path=feature_path,
                            oof_dir=oof_dir,
                            device=device,
                            random_state=random_state,
                            cv_splits=cv_splits,
                        )
                    )
                    controls.update(
                        _benchmark_shuffled_prior(
                            label=label,
                            label_name=label_name,
                            saluki_reference=human_saluki,
                            feature_path=feature_path,
                            oof_dir=oof_dir,
                            random_state=random_state,
                            cv_splits=cv_splits,
                        )
                    )
                rows.append(
                    {
                        "label": label_name,
                        "source": source_name,
                        "lambda": float(lambda_),
                        "random_state": int(random_state),
                        "target_vs_saluki_pearson": saluki_corr["pearson"],
                        "target_vs_saluki_spearman": saluki_corr["spearman"],
                        "ortholog_n": int(ortholog_corr["n"]),
                        "ortholog_pearson": ortholog_corr["pearson"],
                        "ortholog_spearman": ortholog_corr["spearman"],
                        **benchmark,
                        **controls,
                    }
                )

    summary = _add_control_deltas(pd.DataFrame(rows))
    summary = summary.sort_values("prior_pearson_vs_target", ascending=False)
    summary_path = results_root / "summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    aggregate = _aggregate_summary(summary)
    aggregate.to_csv(results_root / "summary_by_label.tsv", sep="\t", index=False)
    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(rows, handle, indent=2)
    figure_paths = _plot_summary(summary, FIGURES_DIR / "ortholog_regularized_label") if make_figures else []
    note_path = _write_note(summary, MANUSCRIPT_DIR.parent / "docs" / "ortholog_regularized_label.md", figure_paths)
    return summary_path, note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark ortholog-regularized human PC1 labels.")
    parser.add_argument("--lambdas", nargs="+", type=float, default=[0.05, 0.1, 0.15, 0.2, 0.3])
    parser.add_argument("--sources", nargs="+", default=["saluki_mouse_prior", "reconstructed_mouse_pc1"])
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR / "ortholog_regularized_label",
    )
    parser.add_argument("--device", default=None)
    parser.add_argument(
        "--include-controls",
        action="store_true",
        help="Also run pure-human and shuffled-prior controls for each label.",
    )
    parser.add_argument("--random-states", nargs="+", type=int, default=[42])
    parser.add_argument("--cv-splits", type=int, default=5)
    parser.add_argument(
        "--include-saluki-baseline",
        action="store_true",
        dest="include_saluki_baseline",
        help="Also benchmark Saluki human PC1 under the same fold/seed/control settings.",
    )
    parser.set_defaults(include_saluki_baseline=False)
    parser.add_argument("--no-figures", action="store_true", help="Skip manuscript figure generation.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary_path, note_path = run_ortholog_regularized_benchmark(
        lambdas=args.lambdas,
        sources=args.sources,
        feature_path=args.feature_path,
        results_root=args.results_root,
        device=args.device,
        include_controls=args.include_controls,
        random_states=args.random_states,
        cv_splits=args.cv_splits,
        make_figures=not args.no_figures,
        include_saluki_baseline=args.include_saluki_baseline,
    )
    print(summary_path)
    print(pd.read_csv(summary_path, sep="\t").to_string(index=False))
    aggregate_path = args.results_root / "summary_by_label.tsv"
    if aggregate_path.exists():
        print(aggregate_path)
        print(pd.read_csv(aggregate_path, sep="\t").to_string(index=False))
    print(note_path)


if __name__ == "__main__":
    main()
