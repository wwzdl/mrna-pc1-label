from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import FIGURES_DIR, MANUSCRIPT_DIR, METADATA_DIR, RESULTS_DIR
from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.pc1_consensus import compute_pc1_consensus
from mrna_half_life_paper.preprocess import PreprocessConfig, iterative_pca_impute, preprocess_matrix
from mrna_half_life_paper.real_data import extract_real_data
from mrna_half_life_paper.run_workflow import run_workflow
from mrna_half_life_paper.saluki_naming import saluki_label_column, saluki_label_filename


def _load_species_inputs(species: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    matrix = load_matrix(RESULTS_DIR.parents[0] / "data" / "processed" / "real_data" / species / "moesm2_raw_sparse_matrix.tsv")
    metadata = load_metadata(METADATA_DIR / f"real_{species}_samples.tsv")
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)
    matrix = matrix.loc[matrix.notna().sum(axis=1) >= 3].copy()

    saluki_frame = pd.read_csv(
        RESULTS_DIR.parents[0] / "data" / "processed" / "real_data" / species / saluki_label_filename(species),
        sep="\t",
    )
    saluki_col = saluki_label_column(species)
    saluki_label = saluki_frame.set_index("gene_id")[saluki_col].astype(float)
    return matrix, metadata, saluki_label


def _compute_sample_pca(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    apply_preprocess: bool,
) -> pd.DataFrame:
    if apply_preprocess:
        plot_matrix = preprocess_matrix(matrix, PreprocessConfig(apply_log=False))
    else:
        plot_matrix = iterative_pca_impute(matrix, n_components=5, max_iter=10, tol=1e-5)

    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(plot_matrix.T.to_numpy(dtype=float))
    pca_df = pd.DataFrame(
        {
            "sample_id": plot_matrix.columns,
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
        }
    )
    return pca_df.merge(metadata, on="sample_id", how="left")


def _plot_sample_pca(df: pd.DataFrame, hue_col: str, title: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, 6.5))
    sns.scatterplot(data=df, x="PC1", y="PC2", hue=hue_col, s=55, edgecolor="none")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def _plot_influence_bar(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    path: Path,
    color: str,
) -> None:
    plot_df = df.sort_values(value_col, ascending=True).copy()
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8.5, max(4.5, 0.33 * len(plot_df))))
    sns.barplot(data=plot_df, x=value_col, y="study", color=color)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def compute_study_influence(species: str) -> dict[str, object]:
    matrix, metadata, saluki_reference_pc1 = _load_species_inputs(species)
    processed_full = preprocess_matrix(matrix, PreprocessConfig(apply_log=False))
    full_pc1 = compute_pc1_consensus(processed_full)
    full_saluki_corr = label_correlation(full_pc1, saluki_reference_pc1)

    rows: list[dict[str, object]] = []
    for study in sorted(metadata["study"].unique().tolist()):
        subset_meta = metadata.loc[metadata["study"] != study].copy()
        subset_samples = subset_meta["sample_id"].tolist()
        subset_matrix = matrix.loc[:, subset_samples]
        subset_processed = preprocess_matrix(subset_matrix, PreprocessConfig(apply_log=False))
        subset_pc1 = compute_pc1_consensus(subset_processed)
        stability = label_correlation(full_pc1, subset_pc1)
        saluki_corr = label_correlation(subset_pc1, saluki_reference_pc1)
        rows.append(
            {
                "study": study,
                "n_samples_removed": int(metadata.loc[metadata["study"] == study].shape[0]),
                "pc1_stability_pearson": stability["pearson"],
                "pc1_stability_spearman": stability["spearman"],
                "saluki_reference_pc1_pearson": saluki_corr["pearson"],
                "saluki_reference_pc1_spearman": saluki_corr["spearman"],
                "delta_saluki_reference_pc1_pearson": saluki_corr["pearson"] - full_saluki_corr["pearson"],
                "delta_saluki_reference_pc1_spearman": saluki_corr["spearman"] - full_saluki_corr["spearman"],
            }
        )

    influence_df = pd.DataFrame(rows).sort_values("pc1_stability_pearson", ascending=True)
    top_influential_study = influence_df.iloc[0]["study"]
    top_influential_delta = float(influence_df.iloc[0]["delta_saluki_reference_pc1_pearson"])

    sample_pca_raw = _compute_sample_pca(matrix, metadata, apply_preprocess=False)
    sample_pca_processed = _compute_sample_pca(matrix, metadata, apply_preprocess=True)

    return {
        "species": species,
        "matrix": matrix,
        "metadata": metadata,
        "saluki_reference_pc1": saluki_reference_pc1,
        "full_pc1": full_pc1,
        "full_saluki_corr": full_saluki_corr,
        "study_influence": influence_df,
        "sample_pca_raw": sample_pca_raw,
        "sample_pca_processed": sample_pca_processed,
        "top_influential_study": top_influential_study,
        "top_influential_delta_saluki_reference_pc1_pearson": top_influential_delta,
    }


def run_study_influence_analysis(
    results_root: Path | None = None,
    figures_root: Path | None = None,
) -> dict[str, dict[str, object]]:
    extract_real_data()

    results_root = results_root if results_root is not None else RESULTS_DIR / "study_influence"
    figures_root = figures_root if figures_root is not None else FIGURES_DIR / "study_influence"
    results_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)

    all_results: dict[str, dict[str, object]] = {}
    summary_rows: list[dict[str, object]] = []

    for species in ("human", "mouse"):
        result = compute_study_influence(species)
        all_results[species] = result

        species_result_dir = results_root / species
        species_figure_dir = figures_root / species
        species_result_dir.mkdir(parents=True, exist_ok=True)
        species_figure_dir.mkdir(parents=True, exist_ok=True)

        result["study_influence"].to_csv(species_result_dir / "study_influence.tsv", sep="\t", index=False)
        result["sample_pca_raw"].to_csv(species_result_dir / "sample_pca_raw.tsv", sep="\t", index=False)
        result["sample_pca_processed"].to_csv(species_result_dir / "sample_pca_processed.tsv", sep="\t", index=False)

        _plot_sample_pca(
            result["sample_pca_raw"],
            hue_col="method",
            title=f"{species.capitalize()} MOESM2 Raw Sample PCA Colored By Method",
            path=species_figure_dir / "raw_sample_pca_by_method.png",
        )
        _plot_sample_pca(
            result["sample_pca_processed"],
            hue_col="method",
            title=f"{species.capitalize()} Processed Sample PCA Colored By Method",
            path=species_figure_dir / "processed_sample_pca_by_method.png",
        )
        _plot_influence_bar(
            result["study_influence"],
            value_col="pc1_stability_pearson",
            title=f"{species.capitalize()} Study Influence On PC1 Stability",
            path=species_figure_dir / "study_influence_pc1_stability.png",
            color="#4C7A9F",
        )
        _plot_influence_bar(
            result["study_influence"],
            value_col="delta_saluki_reference_pc1_pearson",
            title=f"{species.capitalize()} Gain In Saluki-PC1 Correlation After Study Removal",
            path=species_figure_dir / "study_influence_saluki_gain.png",
            color="#C67B48",
        )

        full_corr = result["full_saluki_corr"]
        summary_rows.append(
            {
                "species": species,
                "full_pc1_pearson_vs_saluki": full_corr["pearson"],
                "full_pc1_spearman_vs_saluki": full_corr["spearman"],
                "top_influential_study": result["top_influential_study"],
                "top_influential_delta_saluki_reference_pc1_pearson": result[
                    "top_influential_delta_saluki_reference_pc1_pearson"
                ],
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(results_root / "study_influence_summary.tsv", sep="\t", index=False)
    with open(results_root / "study_influence_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_rows, handle, indent=2)
    return all_results


def run_human_without_gejman(results_root: Path | None = None) -> Path:
    results_root = results_root if results_root is not None else RESULTS_DIR
    matrix, metadata, _saluki_reference_pc1 = _load_species_inputs("human")
    filtered_meta = metadata.loc[metadata["study"] != "Gejman"].copy()
    filtered_meta_path = METADATA_DIR / "real_human_samples_no_Gejman.tsv"
    filtered_meta.to_csv(filtered_meta_path, sep="\t", index=False)

    outdir = results_root / "real_human_moesm2_no_gejman"
    run_workflow(
        matrix_path=RESULTS_DIR.parents[0] / "data" / "processed" / "real_data" / "human" / "moesm2_raw_sparse_matrix.tsv",
        metadata_path=filtered_meta_path,
        outdir=outdir,
        apply_log=False,
    )
    return outdir


def write_cn_report(
    analysis_results: dict[str, dict[str, object]],
    outpath: Path | None = None,
) -> Path:
    outpath = outpath if outpath is not None else MANUSCRIPT_DIR / "paper2_outlier_study_note_cn.md"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    human = analysis_results["human"]
    mouse = analysis_results["mouse"]
    human_top = human["study_influence"].iloc[0]
    mouse_top = mouse["study_influence"].iloc[0]

    text = f"""# Study-aware 审计量化高影响 study 对 mRNA 半衰期共识标签的作用

## 一句话结论

基于 Saluki 公开补充数据，Saluki-label-independent leave-one-study-out PC1 stability 将 human 中的 `{human['top_influential_study']}` 排为对共识标签几何影响最大的 study。随后使用 Saluki human PC1 进行 reference-informed validation；移除该 study 后，重建标签与 Saluki human PC1 的 Pearson 从 {human['full_saluki_corr']['pearson']:.3f} 提升到 {human_top['saluki_reference_pc1_pearson']:.3f}。该排序不使用下游模型分数；它重现并量化 Saluki 数据说明中记录的 study influence。Mouse 数据没有出现同等级的多证据一致性。

## 目前已完成的真实结果

- human：54 个样本，19 个 study，5 种方法。
- mouse：27 个样本，17 个 study，3 种方法。
- human 全量数据的重建 PC1 与 Saluki human PC1 的 Pearson 为 {human['full_saluki_corr']['pearson']:.3f}。
- mouse 全量数据的重建 PC1 与 Saluki mouse PC1 的 Pearson 为 {mouse['full_saluki_corr']['pearson']:.3f}。
- human 中影响最大的 study 是 `{human['top_influential_study']}`，它被剔除后与 Saluki human PC1 的对齐增益为 {human_top['delta_saluki_reference_pc1_pearson']:.3f}。
- mouse 中最强的 study 影响者是 `{mouse['top_influential_study']}`，但其增益只有 {mouse_top['delta_saluki_reference_pc1_pearson']:.3f}，远小于 human。

## 这个题目为什么适合快速出文章

- 不需要训练深度模型，全部基于公开补充数据即可完成。
- 不需要再采集实验数据，工作量主要是统计分析、稳定性验证和图表整理。
- 结论很聚焦：不是重新发明 Saluki，而是说明 study influence 如何改变半衰期共识标签，以及如何用稳定性筛查和外部验证分层审计。
- human 和 mouse 可以形成天然对照：human 存在主导 study influence，mouse 没有同等级的多证据一致性。

## 建议论文标题

`Study-aware auditing of consensus labels in mammalian mRNA half-life datasets`

中文可写为：

`哺乳动物 mRNA 半衰期数据中的 study-aware 共识标签审计`

## 核心图建议

1. raw sample PCA，按 method 着色，显示原始 MOESM2 数据存在明显方法/研究来源聚类。
2. processed sample PCA，显示标准化后方法聚类显著减弱。
3. leave-one-study-out 稳定性柱状图，展示 human 中 `{human['top_influential_study']}` 明显最不稳定。
4. 剔除不同 study 后，与 Saluki PC1 的相关性增益柱状图，突出 human 中 `{human['top_influential_study']}` 的改善最大。

## 文章主体结构

### 1. Introduction

说明 mammalian mRNA half-life 数据来自多个研究，实验策略异质性强。原论文已经意识到 bias 问题，但异常 study 的识别仍主要依赖人工经验。我们提出一个更简单的、可复现的 study influence 分析框架。

### 2. Methods

- 使用论文官方 MOESM2 和 MOESM3 补充数据。
- 对 MOESM2 进行作者论文一致的预处理思路：过滤、插补、标准化、PC1 重建。
- 对每个 study 做 leave-one-study-out。
- 用两类指标评价：
  - 与全量 PC1 的稳定性。
  - 与 Saluki PC1 的一致性变化。

### 3. Results

- human 中存在主导影响 study。
- 去掉该 study 后，人类 PC1 共识标签更接近 Saluki human PC1。
- mouse 中不存在同等级别的多证据一致性，说明高 influence 不自动等于应删除。
- 因此，study influence analysis 可以作为半衰期 meta-analysis 的常规质控步骤。

### 4. Discussion

- 这篇文章不需要否定原文，而是强调原文作者的人工剔除是有数据基础的。
- 我们的贡献是把这一过程形式化、自动化、可复现化。
- 后续可扩展到其它 RNA stability compendium、translation rate compendium、protein half-life compendium。

## 当前可直接使用的结果文件

- `results/study_influence/human/study_influence.tsv`
- `results/study_influence/mouse/study_influence.tsv`
- `figures/study_influence/human/study_influence_pc1_stability.png`
- `figures/study_influence/human/study_influence_saluki_gain.png`
- `figures/study_influence/mouse/study_influence_pc1_stability.png`
- `results/real_human_moesm2_no_gejman/`

## 下一步最省力的投稿路线

- 先把这篇做成 short report / brief communication。
- 不扩展深度学习部分，只做元分析稳定性和质控方法。
- 若要再补一点工作量，增加一个“自动 outlier 剔除前后，和 mouse ortholog label 一致性”的结果即可。
"""
    outpath.write_text(text, encoding="utf-8")
    return outpath


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run study influence analysis for the Saluki half-life compendium.")
    parser.add_argument(
        "--run-pruned-human",
        action="store_true",
        help="Also rerun the human workflow after removing the dominant influential study Gejman.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_results = run_study_influence_analysis()
    report_path = write_cn_report(analysis_results)
    print(report_path)
    if args.run_pruned_human:
        outdir = run_human_without_gejman()
        print(outdir)


if __name__ == "__main__":
    main()
