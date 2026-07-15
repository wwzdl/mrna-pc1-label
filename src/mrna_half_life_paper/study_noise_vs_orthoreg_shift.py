from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mrna_half_life_paper.config import FIGURES_DIR, MANUSCRIPT_DIR, METADATA_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.ortholog_label_distance import run_ortholog_label_distance
from mrna_half_life_paper.ortholog_regularized_label import _load_saluki_label, _load_human_base_label, _zscore
from mrna_half_life_paper.preprocess import PreprocessConfig, preprocess_matrix, zscore_by_sample


def _metrics(
    left: pd.Series,
    right: pd.Series,
    *,
    comparison_type: str,
    cohort: str,
    left_name: str,
    right_name: str,
) -> dict[str, float | int | str]:
    aligned = pd.concat([left.rename("left"), right.rename("right")], axis=1).dropna()
    diff = aligned["left"] - aligned["right"]
    abs_diff = diff.abs()
    return {
        "comparison_type": comparison_type,
        "cohort": cohort,
        "left": left_name,
        "right": right_name,
        "n": int(aligned.shape[0]),
        "pearson": float(aligned["left"].corr(aligned["right"], method="pearson")),
        "spearman": float(aligned["left"].corr(aligned["right"], method="spearman")),
        "mse": float(np.mean(np.square(diff))),
        "rmse": float(np.sqrt(np.mean(np.square(diff)))),
        "mae": float(np.mean(abs_diff)),
        "median_abs_error": float(np.median(abs_diff)),
        "p90_abs_error": float(np.quantile(abs_diff, 0.90)),
        "p95_abs_error": float(np.quantile(abs_diff, 0.95)),
        "mean_signed_diff": float(np.mean(diff)),
        "sd_signed_diff": float(np.std(diff, ddof=0)),
    }


def _safe_zscore(series: pd.Series) -> pd.Series:
    values = series.astype(float)
    std = values.std(ddof=0)
    if float(std) == 0.0 or np.isnan(std):
        return values - values.mean()
    return (values - values.mean()) / std


def _align_sign(label: pd.Series, reference: pd.Series) -> pd.Series:
    aligned = pd.concat([label.rename("label"), reference.rename("reference")], axis=1).dropna()
    if aligned.shape[0] < 3:
        return label
    corr = aligned["label"].corr(aligned["reference"], method="pearson")
    if pd.notna(corr) and corr < 0:
        return -label
    return label


def _load_human_moesm2() -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = load_matrix(PROCESSED_DIR / "real_data" / "human" / "moesm2_raw_sparse_matrix.tsv")
    metadata = load_metadata(METADATA_DIR / "real_human_samples.tsv")
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)
    return matrix, metadata


def _study_mean_labels(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    reference: pd.Series,
    gene_universe: pd.Index,
    min_observed_per_study: int,
) -> dict[str, pd.Series]:
    matrix = matrix.loc[matrix.index.intersection(gene_universe)].copy()
    z_matrix = zscore_by_sample(matrix)
    labels: dict[str, pd.Series] = {}
    for study in sorted(metadata["study"].unique().tolist()):
        samples = metadata.loc[metadata["study"] == study, "sample_id"].tolist()
        study_values = z_matrix.loc[:, samples]
        observed = study_values.notna().sum(axis=1)
        label = study_values.mean(axis=1, skipna=True).where(observed >= min_observed_per_study).dropna()
        label = _safe_zscore(label).rename(study)
        label = _align_sign(label, reference)
        labels[study] = label
    return labels


def _preprocessed_study_mean_labels(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    reference: pd.Series,
    gene_universe: pd.Index,
) -> dict[str, pd.Series]:
    processed = preprocess_matrix(
        matrix,
        PreprocessConfig(apply_log=False, min_observed_per_gene=10),
    )
    processed = processed.loc[processed.index.intersection(gene_universe)].copy()
    labels: dict[str, pd.Series] = {}
    for study in sorted(metadata["study"].unique().tolist()):
        samples = metadata.loc[metadata["study"] == study, "sample_id"].tolist()
        samples = [sample for sample in samples if sample in processed.columns]
        if not samples:
            continue
        label = processed.loc[:, samples].mean(axis=1).dropna()
        label = _safe_zscore(label).rename(study)
        label = _align_sign(label, reference)
        labels[study] = label
    return labels


def _pairwise_label_metrics(
    labels: dict[str, pd.Series],
    *,
    cohort: str,
    min_common_genes: int,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for left_name, right_name in combinations(sorted(labels.keys()), 2):
        left = labels[left_name]
        right = labels[right_name]
        aligned = pd.concat([left, right], axis=1).dropna()
        if aligned.shape[0] < min_common_genes:
            continue
        rows.append(
            _metrics(
                aligned.iloc[:, 0],
                aligned.iloc[:, 1],
                comparison_type="study_mean_label_pair",
                cohort=cohort,
                left_name=left_name,
                right_name=right_name,
            )
        )
    return pd.DataFrame(rows)


def _study_residuals(
    labels: dict[str, pd.Series],
    reference: pd.Series,
    *,
    cohort: str,
    min_common_genes: int,
) -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    rows: list[dict[str, float | int | str]] = []
    residuals: dict[str, pd.Series] = {}
    for study, label in sorted(labels.items()):
        aligned = pd.concat([label.rename("study_label"), reference.rename("reference")], axis=1).dropna()
        if aligned.shape[0] < min_common_genes:
            continue
        residual = (aligned["study_label"] - aligned["reference"]).rename(study)
        residuals[study] = residual
        rows.append(
            _metrics(
                aligned["study_label"],
                aligned["reference"],
                comparison_type="study_label_vs_human_no_gejman_pc1",
                cohort=cohort,
                left_name=study,
                right_name="human_no_gejman_pc1",
            )
        )
    return pd.DataFrame(rows), residuals


def _pairwise_residual_metrics(
    residuals: dict[str, pd.Series],
    *,
    cohort: str,
    min_common_genes: int,
) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for left_name, right_name in combinations(sorted(residuals.keys()), 2):
        left = residuals[left_name]
        right = residuals[right_name]
        aligned = pd.concat([left, right], axis=1).dropna()
        if aligned.shape[0] < min_common_genes:
            continue
        rows.append(
            _metrics(
                aligned.iloc[:, 0],
                aligned.iloc[:, 1],
                comparison_type="study_residual_pair",
                cohort=cohort,
                left_name=left_name,
                right_name=right_name,
            )
        )
    return pd.DataFrame(rows)


def _summarize_distribution(
    frame: pd.DataFrame,
    *,
    comparison_type: str,
    cohort: str,
    source: str,
) -> list[dict[str, float | int | str]]:
    sub = frame.loc[(frame["comparison_type"] == comparison_type) & (frame["cohort"] == cohort)].copy()
    rows: list[dict[str, float | int | str]] = []
    for metric in ("pearson", "spearman", "rmse", "mae", "median_abs_error", "p90_abs_error"):
        values = sub[metric].dropna().astype(float)
        if values.empty:
            continue
        rows.append(
            {
                "source": source,
                "comparison_type": comparison_type,
                "cohort": cohort,
                "metric": metric,
                "n_comparisons": int(values.shape[0]),
                "mean": float(values.mean()),
                "median": float(values.median()),
                "q25": float(values.quantile(0.25)),
                "q75": float(values.quantile(0.75)),
                "min": float(values.min()),
                "max": float(values.max()),
            }
        )
    return rows


def _load_lambda_shift(results_root: Path, lambda_: float) -> pd.Series:
    distance_root = RESULTS_DIR / "ortholog_label_distance"
    if not (distance_root / "distance_summary.tsv").exists():
        run_ortholog_label_distance(results_root=distance_root, lambda_=lambda_)
    summary = pd.read_csv(distance_root / "distance_summary.tsv", sep="\t")
    row = summary.loc[
        (summary["comparison"] == "orthoreg_lambda0p10_vs_human_no_gejman_pc1")
        & (summary["subset"] == "all_one2one")
    ].iloc[0]
    out = row.copy()
    out["comparison_type"] = "orthoreg_lambda0p10_shift"
    return out


def _load_gene_universe() -> pd.Index:
    path = RESULTS_DIR / "ortholog_label_distance" / "ortholog_label_values.tsv"
    if not path.exists():
        run_ortholog_label_distance(results_root=RESULTS_DIR / "ortholog_label_distance", lambda_=0.10)
    frame = pd.read_csv(path, sep="\t")
    return pd.Index(frame["human_gene_id"].astype(str).unique())


def _write_note(
    summary: pd.DataFrame,
    lambda_shift: pd.Series,
    outpath: Path,
) -> Path:
    def pick(source: str, comparison_type: str, cohort: str, metric: str) -> pd.Series:
        return summary.loc[
            (summary["source"] == source)
            & (summary["comparison_type"] == comparison_type)
            & (summary["cohort"] == cohort)
            & (summary["metric"] == metric)
        ].iloc[0]

    pair_pearson = pick("study_mean_labels", "study_mean_label_pair", "no_gejman", "pearson")
    pair_rmse = pick("study_mean_labels", "study_mean_label_pair", "no_gejman", "rmse")
    pair_mae = pick("study_mean_labels", "study_mean_label_pair", "no_gejman", "mae")
    residual_pearson = pick("study_residuals", "study_residual_pair", "no_gejman", "pearson")
    residual_rmse = pick("study_residuals", "study_label_vs_human_no_gejman_pc1", "no_gejman", "rmse")
    residual_mae = pick("study_residuals", "study_label_vs_human_no_gejman_pc1", "no_gejman", "mae")
    multistudy_pair_rmse = pick(
        "study_mean_labels",
        "study_mean_label_pair",
        "no_gejman_min2samples",
        "rmse",
    )
    multistudy_pair_mae = pick(
        "study_mean_labels",
        "study_mean_label_pair",
        "no_gejman_min2samples",
        "mae",
    )
    multistudy_residual_pearson = pick(
        "study_residuals",
        "study_residual_pair",
        "no_gejman_min2samples",
        "pearson",
    )
    full_preprocess_pair_pearson = pick(
        "study_mean_labels",
        "study_mean_label_pair",
        "no_gejman_full_preprocess",
        "pearson",
    )
    full_preprocess_pair_rmse = pick(
        "study_mean_labels",
        "study_mean_label_pair",
        "no_gejman_full_preprocess",
        "rmse",
    )
    full_preprocess_pair_mae = pick(
        "study_mean_labels",
        "study_mean_label_pair",
        "no_gejman_full_preprocess",
        "mae",
    )

    shift_rmse = float(lambda_shift["rmse"])
    shift_mae = float(lambda_shift["mae"])
    rmse_fold = float(pair_rmse["median"]) / shift_rmse
    mae_fold = float(pair_mae["median"]) / shift_mae
    residual_rmse_fold = float(residual_rmse["median"]) / shift_rmse
    residual_mae_fold = float(residual_mae["median"]) / shift_mae

    lines = [
        "# Study Noise Versus Ortholog-Regularized Label Shift",
        "",
        "## Purpose",
        "",
        "比较 `lambda=0.10` ortholog-regularized target 相对原 human no-Gejman PC1 的微小改变量，",
        "与 human compendium 中不同 study 之间的实验差异是否处在同一量级。",
        "",
        "## Key Results",
        "",
        f"- `lambda=0.10` 混合前后标签几乎不改变 human 标签主体结构：N={int(lambda_shift['n'])}，"
        f"Pearson={float(lambda_shift['pearson']):.4f}，RMSE={shift_rmse:.4f}，MAE={shift_mae:.4f}。",
        f"- 在排除 Gejman 后的 normal-study 对照中，study mean label 两两比较的中位 Pearson 为 "
        f"{float(pair_pearson['median']):.4f}（IQR {float(pair_pearson['q25']):.4f}-{float(pair_pearson['q75']):.4f}），"
        f"中位 RMSE={float(pair_rmse['median']):.4f}，中位 MAE={float(pair_mae['median']):.4f}。",
        f"- 因此，normal studies 之间的中位 label RMSE 约为 `lambda=0.10` 标签 shift 的 {rmse_fold:.1f} 倍，"
        f"中位 MAE 约为 {mae_fold:.1f} 倍。",
        f"- 若把每个 study label 相对 human no-Gejman PC1 的差异视为 study-specific residual，"
        f"不同 study residual 两两相关的中位 Pearson 为 {float(residual_pearson['median']):.4f}，"
        f"说明误差结构并不高度一致。",
        f"- 单个 normal study 相对 human no-Gejman PC1 的中位 RMSE={float(residual_rmse['median']):.4f}，"
        f"MAE={float(residual_mae['median']):.4f}；分别约为 `lambda=0.10` shift 的 "
        f"{residual_rmse_fold:.1f} 倍和 {residual_mae_fold:.1f} 倍。",
        f"- 更保守地只保留 no-Gejman 且样本数不少于 2 的 7 个 study 时，study-pair 中位 "
        f"RMSE={float(multistudy_pair_rmse['median']):.4f}，MAE={float(multistudy_pair_mae['median']):.4f}，"
        f"残差两两相关中位 Pearson={float(multistudy_residual_pearson['median']):.4f}；"
        "结论仍不依赖单样本 study。",
        f"- 作为预处理敏感性检查，先对 no-Gejman 样本执行完整 Saluki-like preprocessing 后再按 study 取均值，"
        f"study-pair 中位 Pearson 升至 {float(full_preprocess_pair_pearson['median']):.4f}，"
        f"但中位 RMSE={float(full_preprocess_pair_rmse['median']):.4f}，"
        f"MAE={float(full_preprocess_pair_mae['median']):.4f}，仍远大于 `lambda=0.10` shift。",
        "",
        "## Interpretation",
        "",
        "- 可以说：`lambda=0.10` 是轻度 target regularization，而不是重写 human 标签。",
        "- 更准确的写法不是“完全可以忽略”，而是“其幅度远小于 normal study-to-study variability，"
        "因此不太可能改变 human 标签的主要排序结构”。",
        "- 这组结果可作为补充材料防守：ortholog regularization 的收益来自对 noisy consensus target 的轻度跨物种收缩，"
        "而不是把 human target 换成 mouse target。",
        "",
    ]
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text("\n".join(lines), encoding="utf-8")
    return outpath


def _plot_summary(summary: pd.DataFrame, lambda_shift: pd.Series, path: Path) -> Path:
    no_gejman = summary.loc[summary["cohort"] == "no_gejman"].copy()
    rows = [
        {
            "comparison": "lambda=0.10 shift",
            "rmse": float(lambda_shift["rmse"]),
            "mae": float(lambda_shift["mae"]),
        }
    ]
    for comparison_type, source, label in [
        ("study_mean_label_pair", "study_mean_labels", "study-study label pair"),
        ("study_label_vs_human_no_gejman_pc1", "study_residuals", "study vs human no-Gejman PC1"),
    ]:
        for metric in ("rmse", "mae"):
            row = no_gejman.loc[
                (no_gejman["source"] == source)
                & (no_gejman["comparison_type"] == comparison_type)
                & (no_gejman["metric"] == metric)
            ].iloc[0]
            existing = next((item for item in rows if item["comparison"] == label), None)
            if existing is None:
                existing = {"comparison": label}
                rows.append(existing)
            existing[metric] = float(row["median"])

    plot_df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(plot_df.shape[0])
    width = 0.36
    ax.bar(x - width / 2, plot_df["rmse"], width=width, label="RMSE", color="#4C7A9F")
    ax.bar(x + width / 2, plot_df["mae"], width=width, label="MAE", color="#C67B48")
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["comparison"], rotation=20, ha="right")
    ax.set_ylabel("z-scored label units")
    ax.set_title("Ortholog-Regularized Shift Versus Study-Level Noise")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)
    return path


def run_study_noise_vs_orthoreg_shift(
    results_root: Path | None = None,
    figures_root: Path | None = None,
    lambda_: float = 0.10,
    min_common_genes: int = 500,
    min_observed_per_study: int = 1,
) -> dict[str, Path]:
    results_root = results_root if results_root is not None else RESULTS_DIR / "study_noise_vs_orthoreg_shift"
    figures_root = figures_root if figures_root is not None else FIGURES_DIR / "study_noise_vs_orthoreg_shift"
    results_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)

    gene_universe = _load_gene_universe()
    matrix, metadata = _load_human_moesm2()
    human_base = _zscore(_load_human_base_label()).rename("human_no_gejman_pc1")
    saluki_label = _zscore(_load_saluki_label("human")).rename("saluki_human_pc1")
    reference = human_base.loc[human_base.index.intersection(gene_universe)].copy()
    saluki_reference = saluki_label.loc[saluki_label.index.intersection(gene_universe)].copy()

    all_labels = _study_mean_labels(
        matrix=matrix,
        metadata=metadata,
        reference=reference,
        gene_universe=gene_universe,
        min_observed_per_study=min_observed_per_study,
    )
    no_gejman_labels = {key: value for key, value in all_labels.items() if key != "Gejman"}
    sample_counts = metadata.groupby("study").size()
    no_gejman_min2_labels = {
        key: value
        for key, value in no_gejman_labels.items()
        if int(sample_counts.get(key, 0)) >= 2
    }
    no_gejman_metadata = metadata.loc[metadata["study"] != "Gejman"].copy()
    no_gejman_samples = no_gejman_metadata["sample_id"].tolist()
    no_gejman_preprocessed_labels = _preprocessed_study_mean_labels(
        matrix=matrix.loc[:, no_gejman_samples],
        metadata=no_gejman_metadata,
        reference=reference,
        gene_universe=gene_universe,
    )

    pair_frames = [
        _pairwise_label_metrics(all_labels, cohort="all", min_common_genes=min_common_genes),
        _pairwise_label_metrics(no_gejman_labels, cohort="no_gejman", min_common_genes=min_common_genes),
        _pairwise_label_metrics(
            no_gejman_min2_labels,
            cohort="no_gejman_min2samples",
            min_common_genes=min_common_genes,
        ),
        _pairwise_label_metrics(
            no_gejman_preprocessed_labels,
            cohort="no_gejman_full_preprocess",
            min_common_genes=min_common_genes,
        ),
    ]
    pair_metrics = pd.concat(pair_frames, ignore_index=True)
    pair_path = results_root / "study_pair_label_metrics.tsv"
    pair_metrics.to_csv(pair_path, sep="\t", index=False)

    residual_frames: list[pd.DataFrame] = []
    residual_pair_frames: list[pd.DataFrame] = []
    residual_label_tables: list[pd.DataFrame] = []
    for cohort, labels in (
        ("all", all_labels),
        ("no_gejman", no_gejman_labels),
        ("no_gejman_min2samples", no_gejman_min2_labels),
        ("no_gejman_full_preprocess", no_gejman_preprocessed_labels),
    ):
        residual_metrics, residuals = _study_residuals(
            labels,
            reference,
            cohort=cohort,
            min_common_genes=min_common_genes,
        )
        residual_frames.append(residual_metrics)
        residual_pair_frames.append(
            _pairwise_residual_metrics(residuals, cohort=cohort, min_common_genes=min_common_genes)
        )
        for study, residual in residuals.items():
            residual_label_tables.append(
                residual.rename("residual")
                .to_frame()
                .reset_index()
                .rename(columns={"index": "gene_id"})
                .assign(study=study, cohort=cohort)
            )

    residual_metrics = pd.concat(residual_frames, ignore_index=True)
    residual_pair_metrics = pd.concat(residual_pair_frames, ignore_index=True)
    residual_path = results_root / "study_residual_vs_reference_metrics.tsv"
    residual_pair_path = results_root / "study_residual_pair_metrics.tsv"
    residual_metrics.to_csv(residual_path, sep="\t", index=False)
    residual_pair_metrics.to_csv(residual_pair_path, sep="\t", index=False)
    pd.concat(residual_label_tables, ignore_index=True).to_csv(
        results_root / "study_residual_values.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    lambda_shift = _load_lambda_shift(results_root, lambda_)
    summary_rows: list[dict[str, float | int | str]] = []
    for cohort in ("all", "no_gejman", "no_gejman_min2samples", "no_gejman_full_preprocess"):
        summary_rows.extend(
            _summarize_distribution(
                pair_metrics,
                comparison_type="study_mean_label_pair",
                cohort=cohort,
                source="study_mean_labels",
            )
        )
        summary_rows.extend(
            _summarize_distribution(
                residual_metrics,
                comparison_type="study_label_vs_human_no_gejman_pc1",
                cohort=cohort,
                source="study_residuals",
            )
        )
        summary_rows.extend(
            _summarize_distribution(
                residual_pair_metrics,
                comparison_type="study_residual_pair",
                cohort=cohort,
                source="study_residuals",
            )
        )
    summary = pd.DataFrame(summary_rows)
    summary_path = results_root / "noise_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    lambda_path = results_root / "lambda0p10_shift_metrics.tsv"
    lambda_shift.to_frame().T.to_csv(lambda_path, sep="\t", index=False)

    reference_alignment = pd.DataFrame(
        [
            _metrics(
                human_base,
                saluki_reference,
                comparison_type="human_no_gejman_pc1_vs_saluki_human_pc1",
                cohort="ortholog_gene_universe",
                left_name="human_no_gejman_pc1",
                right_name="saluki_human_pc1",
            )
        ]
    )
    reference_alignment_path = results_root / "reference_alignment_metrics.tsv"
    reference_alignment.to_csv(reference_alignment_path, sep="\t", index=False)

    note_path = _write_note(summary, lambda_shift, MANUSCRIPT_DIR.parent / "docs" / "study_noise_vs_orthoreg_shift.md")
    figure_path = _plot_summary(summary, lambda_shift, figures_root / "study_noise_vs_orthoreg_shift.png")

    metadata_out = {
        "lambda": lambda_,
        "min_common_genes": min_common_genes,
        "min_observed_per_study": min_observed_per_study,
        "gene_universe": "one-to-one ortholog human genes used in ortholog_label_distance",
        "n_gene_universe": int(len(gene_universe)),
        "study_label_definition": "mean of sample-wise z-scored MOESM2 values within each study, z-scored and sign-aligned to human no-Gejman PC1",
        "reference": "human no-Gejman Saluki-like PC1, z-scored",
    }
    metadata_path = results_root / "analysis_metadata.json"
    metadata_path.write_text(json.dumps(metadata_out, indent=2), encoding="utf-8")

    return {
        "pair_metrics": pair_path,
        "residual_metrics": residual_path,
        "residual_pair_metrics": residual_pair_path,
        "summary": summary_path,
        "lambda_shift": lambda_path,
        "reference_alignment": reference_alignment_path,
        "note": note_path,
        "figure": figure_path,
        "metadata": metadata_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare study-to-study noise with lambda=0.10 label shift.")
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "study_noise_vs_orthoreg_shift")
    parser.add_argument("--figures-root", type=Path, default=FIGURES_DIR / "study_noise_vs_orthoreg_shift")
    parser.add_argument("--lambda", dest="lambda_", type=float, default=0.10)
    parser.add_argument("--min-common-genes", type=int, default=500)
    parser.add_argument("--min-observed-per-study", type=int, default=1)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_study_noise_vs_orthoreg_shift(
        results_root=args.results_root,
        figures_root=args.figures_root,
        lambda_=args.lambda_,
        min_common_genes=args.min_common_genes,
        min_observed_per_study=args.min_observed_per_study,
    )
    for path in outputs.values():
        print(path)
        if path.suffix == ".tsv":
            print(pd.read_csv(path, sep="\t").head(20).to_string(index=False))


if __name__ == "__main__":
    main()
