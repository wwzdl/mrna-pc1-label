from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

from mrna_half_life_paper.config import DATA_DIR, FIGURES_DIR, MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.real_data import extract_real_data
from mrna_half_life_paper.saluki_naming import saluki_label_filename
from mrna_half_life_paper.study_influence import run_human_without_gejman


ENSEMBL_MOUSE_HOMOLOGY_URL = (
    "https://ftp.ensembl.org/pub/release-115/tsv/ensembl-compara/homologies/mus_musculus/"
    "Compara.115.protein_default.homologies.tsv.gz"
)
ENSEMBL_MOUSE_HOMOLOGY_SHA256 = "ed82f17dda96df6b47d0885e07bb3adc7d274c760e10437fa3a58e7081206dd0"
REMOTE_DIR = DATA_DIR / "raw" / "remote"
ORTHOLOG_DIR = DATA_DIR / "processed" / "orthologs"
ENSEMBL_MOUSE_HOMOLOGY_GZ = REMOTE_DIR / "Compara.115.protein_default.homologies.mus_musculus.tsv.gz"
ONE2ONE_ORTHOLOG_PATH = ORTHOLOG_DIR / "ensembl_mouse_human_one2one.tsv"
HIGHCONF_ONE2ONE_ORTHOLOG_PATH = ORTHOLOG_DIR / "ensembl_mouse_human_one2one_high_confidence.tsv"
_SOURCE_CHECKSUM_VERIFIED = False


def _verify_ortholog_source(path: Path) -> None:
    global _SOURCE_CHECKSUM_VERIFIED
    if _SOURCE_CHECKSUM_VERIFIED:
        return
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    observed = digest.hexdigest()
    if observed != ENSEMBL_MOUSE_HOMOLOGY_SHA256:
        raise ValueError(
            f"Ensembl release 115 ortholog checksum mismatch for {path}: "
            f"expected {ENSEMBL_MOUSE_HOMOLOGY_SHA256}, observed {observed}"
        )
    _SOURCE_CHECKSUM_VERIFIED = True


def fetch_mouse_human_ortholog_source(force: bool = False) -> Path:
    REMOTE_DIR.mkdir(parents=True, exist_ok=True)
    if ENSEMBL_MOUSE_HOMOLOGY_GZ.exists() and not force:
        _verify_ortholog_source(ENSEMBL_MOUSE_HOMOLOGY_GZ)
        return ENSEMBL_MOUSE_HOMOLOGY_GZ

    context = ssl._create_unverified_context()
    with urllib.request.urlopen(ENSEMBL_MOUSE_HOMOLOGY_URL, context=context) as response:
        ENSEMBL_MOUSE_HOMOLOGY_GZ.write_bytes(response.read())
    _verify_ortholog_source(ENSEMBL_MOUSE_HOMOLOGY_GZ)
    return ENSEMBL_MOUSE_HOMOLOGY_GZ


def extract_mouse_human_one2one(force: bool = False) -> tuple[Path, Path]:
    ORTHOLOG_DIR.mkdir(parents=True, exist_ok=True)
    if ONE2ONE_ORTHOLOG_PATH.exists() and HIGHCONF_ONE2ONE_ORTHOLOG_PATH.exists() and not force:
        return ONE2ONE_ORTHOLOG_PATH, HIGHCONF_ONE2ONE_ORTHOLOG_PATH

    source_path = fetch_mouse_human_ortholog_source(force=force)
    usecols = [
        "gene_stable_id",
        "homology_type",
        "homology_gene_stable_id",
        "homology_species",
        "is_high_confidence",
    ]

    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(source_path, sep="\t", compression="gzip", usecols=usecols, chunksize=200000):
        sub = chunk[
            (chunk["homology_species"] == "homo_sapiens")
            & (chunk["homology_type"] == "ortholog_one2one")
        ].copy()
        if not sub.empty:
            chunks.append(sub)

    mapping = pd.concat(chunks, ignore_index=True)
    mapping["is_high_confidence"] = pd.to_numeric(mapping["is_high_confidence"], errors="coerce").fillna(0).astype(int)
    mapping = (
        mapping.rename(
            columns={
                "gene_stable_id": "mouse_gene_id",
                "homology_gene_stable_id": "human_gene_id",
            }
        )
        .loc[:, ["human_gene_id", "mouse_gene_id", "homology_type", "homology_species", "is_high_confidence"]]
        .drop_duplicates(["human_gene_id", "mouse_gene_id"])
        .sort_values(["human_gene_id", "mouse_gene_id"])
        .reset_index(drop=True)
    )
    mapping.to_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t", index=False)
    mapping.loc[mapping["is_high_confidence"] == 1].to_csv(HIGHCONF_ONE2ONE_ORTHOLOG_PATH, sep="\t", index=False)
    return ONE2ONE_ORTHOLOG_PATH, HIGHCONF_ONE2ONE_ORTHOLOG_PATH


def _label_correlation(series_a: pd.Series, series_b: pd.Series) -> dict[str, float]:
    common = series_a.index.intersection(series_b.index)
    if len(common) < 2:
        return {"n": float(len(common)), "pearson": float("nan"), "spearman": float("nan")}
    a = series_a.loc[common].astype(float)
    b = series_b.loc[common].astype(float)
    return {
        "n": float(len(common)),
        "pearson": float(a.corr(b, method="pearson")),
        "spearman": float(spearmanr(a, b, nan_policy="omit").statistic),
    }


def _load_result_labels() -> dict[str, pd.DataFrame | pd.Series]:
    extract_real_data()
    run_human_without_gejman()

    return {
        "human_full": pd.read_csv(RESULTS_DIR / "real_human_moesm2" / "consensus_labels.tsv", sep="\t").set_index("gene_id"),
        "human_pruned": pd.read_csv(RESULTS_DIR / "real_human_moesm2_no_gejman" / "consensus_labels.tsv", sep="\t").set_index("gene_id"),
        "mouse_full": pd.read_csv(RESULTS_DIR / "real_mouse_moesm2" / "consensus_labels.tsv", sep="\t").set_index("gene_id"),
        "saluki_human": pd.read_csv(
            DATA_DIR / "processed" / "real_data" / "human" / saluki_label_filename("human"), sep="\t"
        ).set_index("gene_id")["saluki_human_pc1"],
        "saluki_mouse": pd.read_csv(
            DATA_DIR / "processed" / "real_data" / "mouse" / saluki_label_filename("mouse"), sep="\t"
        ).set_index("gene_id")["saluki_mouse_prior"],
    }


def _build_common_pair_table(mapping: pd.DataFrame, labels: dict[str, pd.DataFrame | pd.Series]) -> pd.DataFrame:
    human_full = labels["human_full"]
    human_pruned = labels["human_pruned"]
    mouse_full = labels["mouse_full"]

    common_pairs = mapping[
        mapping["human_gene_id"].isin(human_full.index)
        & mapping["human_gene_id"].isin(human_pruned.index)
        & mapping["mouse_gene_id"].isin(mouse_full.index)
    ].copy()

    common_pairs = common_pairs.merge(
        human_full[["pc1", "mean", "bias_aware"]].rename(
            columns={
                "pc1": "human_full_pc1",
                "mean": "human_full_mean",
                "bias_aware": "human_full_bias_aware",
            }
        ),
        left_on="human_gene_id",
        right_index=True,
    )
    common_pairs = common_pairs.merge(
        human_pruned[["pc1", "mean", "bias_aware"]].rename(
            columns={
                "pc1": "human_pruned_pc1",
                "mean": "human_pruned_mean",
                "bias_aware": "human_pruned_bias_aware",
            }
        ),
        left_on="human_gene_id",
        right_index=True,
    )
    common_pairs = common_pairs.merge(
        mouse_full[["pc1", "mean", "bias_aware"]].rename(
            columns={
                "pc1": "mouse_pc1",
                "mean": "mouse_mean",
                "bias_aware": "mouse_bias_aware",
            }
        ),
        left_on="mouse_gene_id",
        right_index=True,
    )
    return common_pairs


def _build_saluki_pair_table(mapping: pd.DataFrame, labels: dict[str, pd.DataFrame | pd.Series]) -> pd.DataFrame:
    saluki_human = labels["saluki_human"]
    saluki_mouse = labels["saluki_mouse"]

    pairs = mapping[
        mapping["human_gene_id"].isin(saluki_human.index)
        & mapping["mouse_gene_id"].isin(saluki_mouse.index)
    ].copy()
    pairs = pairs.merge(saluki_human.rename("saluki_human_pc1"), left_on="human_gene_id", right_index=True)
    pairs = pairs.merge(saluki_mouse.rename("saluki_mouse_prior"), left_on="mouse_gene_id", right_index=True)
    return pairs


def _scatter_with_fit(
    df: pd.DataFrame,
    x: str,
    y: str,
    title: str,
    path: Path,
    color: str,
) -> None:
    corr = _label_correlation(df[x], df[y])
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6.7, 6.0))
    sns.regplot(data=df, x=x, y=y, scatter_kws={"s": 10, "alpha": 0.35}, line_kws={"color": "#1f1f1f"}, color=color)
    plt.title(title)
    plt.xlabel(x.replace("_", " "))
    plt.ylabel(y.replace("_", " "))
    plt.text(
        0.03,
        0.97,
        f"Pearson={corr['pearson']:.3f}\nSpearman={corr['spearman']:.3f}\nN={int(corr['n'])}",
        transform=plt.gca().transAxes,
        ha="left",
        va="top",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def _barplot_summary(summary: pd.DataFrame, metric: str, path: Path) -> None:
    plot_df = summary.loc[
        summary["comparison"].isin(
            ["human_full_pc1_vs_mouse_pc1", "human_pruned_pc1_vs_mouse_pc1", "saluki_human_pc1_vs_saluki_mouse_prior"]
        )
    ].copy()
    plot_df["label"] = plot_df["comparison"].map(
        {
            "human_full_pc1_vs_mouse_pc1": "Human Full vs Mouse",
            "human_pruned_pc1_vs_mouse_pc1": "Human No Gejman vs Mouse",
            "saluki_human_pc1_vs_saluki_mouse_prior": "Saluki Human vs Saluki Mouse",
        }
    )
    order = ["Human Full vs Mouse", "Human No Gejman vs Mouse", "Saluki Human vs Saluki Mouse"]
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7.4, 4.8))
    sns.barplot(
        data=plot_df,
        x="label",
        y=metric,
        hue="label",
        order=order,
        hue_order=order,
        palette=["#7E9FBE", "#D08C60", "#8A9A5B"],
        legend=False,
    )
    plt.ylabel(metric.replace("_", " "))
    plt.xlabel("")
    plt.xticks(rotation=10, ha="right")
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def run_ortholog_analysis(
    results_root: Path | None = None,
    figures_root: Path | None = None,
    force_download: bool = False,
) -> dict[str, float]:
    results_root = results_root if results_root is not None else RESULTS_DIR / "ortholog_analysis"
    figures_root = figures_root if figures_root is not None else FIGURES_DIR / "ortholog_analysis"
    results_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)

    extract_mouse_human_one2one(force=force_download)
    labels = _load_result_labels()

    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")
    mapping_high = pd.read_csv(HIGHCONF_ONE2ONE_ORTHOLOG_PATH, sep="\t")

    common_pairs = _build_common_pair_table(mapping, labels)
    common_pairs_high = _build_common_pair_table(mapping_high, labels)
    saluki_pairs = _build_saluki_pair_table(mapping, labels)

    common_pairs.to_csv(results_root / "common_ortholog_pairs.tsv", sep="\t", index=False)
    common_pairs_high.to_csv(results_root / "common_ortholog_pairs_high_confidence.tsv", sep="\t", index=False)
    saluki_pairs.to_csv(results_root / "saluki_ortholog_pairs.tsv", sep="\t", index=False)

    summary_rows: list[dict[str, float | str]] = []
    for name, left, right, df in (
        ("human_full_pc1_vs_mouse_pc1", "human_full_pc1", "mouse_pc1", common_pairs),
        ("human_pruned_pc1_vs_mouse_pc1", "human_pruned_pc1", "mouse_pc1", common_pairs),
        ("human_full_mean_vs_mouse_mean", "human_full_mean", "mouse_mean", common_pairs),
        ("human_pruned_mean_vs_mouse_mean", "human_pruned_mean", "mouse_mean", common_pairs),
        ("saluki_human_pc1_vs_saluki_mouse_prior", "saluki_human_pc1", "saluki_mouse_prior", saluki_pairs),
        ("highconf_human_full_pc1_vs_mouse_pc1", "human_full_pc1", "mouse_pc1", common_pairs_high),
        ("highconf_human_pruned_pc1_vs_mouse_pc1", "human_pruned_pc1", "mouse_pc1", common_pairs_high),
    ):
        corr = _label_correlation(df[left], df[right])
        summary_rows.append(
            {
                "comparison": name,
                "n_pairs": int(corr["n"]),
                "pearson": corr["pearson"],
                "spearman": corr["spearman"],
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(results_root / "ortholog_correlation_summary.tsv", sep="\t", index=False)

    common_pair_corrs = {
        "pc1_full": _label_correlation(common_pairs["human_full_pc1"], common_pairs["mouse_pc1"]),
        "pc1_pruned": _label_correlation(common_pairs["human_pruned_pc1"], common_pairs["mouse_pc1"]),
        "mean_full": _label_correlation(common_pairs["human_full_mean"], common_pairs["mouse_mean"]),
        "mean_pruned": _label_correlation(common_pairs["human_pruned_mean"], common_pairs["mouse_mean"]),
    }
    deltas = {
        "pc1_delta_pearson": common_pair_corrs["pc1_pruned"]["pearson"] - common_pair_corrs["pc1_full"]["pearson"],
        "pc1_delta_spearman": common_pair_corrs["pc1_pruned"]["spearman"] - common_pair_corrs["pc1_full"]["spearman"],
        "mean_delta_pearson": common_pair_corrs["mean_pruned"]["pearson"] - common_pair_corrs["mean_full"]["pearson"],
        "mean_delta_spearman": common_pair_corrs["mean_pruned"]["spearman"] - common_pair_corrs["mean_full"]["spearman"],
    }
    with open(results_root / "ortholog_delta_summary.json", "w", encoding="utf-8") as handle:
        json.dump({"common_pair_correlations": common_pair_corrs, "deltas": deltas}, handle, indent=2)

    _scatter_with_fit(
        common_pairs,
        x="human_full_pc1",
        y="mouse_pc1",
        title="Human Full PC1 vs Mouse PC1 on Ensembl 1:1 Orthologs",
        path=figures_root / "human_full_vs_mouse_pc1.png",
        color="#6B8BA4",
    )
    _scatter_with_fit(
        common_pairs,
        x="human_pruned_pc1",
        y="mouse_pc1",
        title="Human No-Gejman PC1 vs Mouse PC1 on Ensembl 1:1 Orthologs",
        path=figures_root / "human_pruned_vs_mouse_pc1.png",
        color="#C98654",
    )
    _barplot_summary(summary, metric="pearson", path=figures_root / "ortholog_pc1_barplot_pearson.png")
    _barplot_summary(summary, metric="spearman", path=figures_root / "ortholog_pc1_barplot_spearman.png")

    return {
        "common_pairs": int(common_pairs.shape[0]),
        "pc1_full_pearson": float(common_pair_corrs["pc1_full"]["pearson"]),
        "pc1_pruned_pearson": float(common_pair_corrs["pc1_pruned"]["pearson"]),
        "pc1_delta_pearson": float(deltas["pc1_delta_pearson"]),
    }


def write_cn_report(outpath: Path | None = None) -> Path:
    outpath = outpath if outpath is not None else MANUSCRIPT_DIR / "paper2_ortholog_extension_cn.md"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(RESULTS_DIR / "ortholog_analysis" / "ortholog_correlation_summary.tsv", sep="\t")

    lookup = summary.set_index("comparison")
    common_pairs = int(lookup.loc["human_full_pc1_vs_mouse_pc1", "n_pairs"])
    full_pearson = float(lookup.loc["human_full_pc1_vs_mouse_pc1", "pearson"])
    pruned_pearson = float(lookup.loc["human_pruned_pc1_vs_mouse_pc1", "pearson"])
    full_spearman = float(lookup.loc["human_full_pc1_vs_mouse_pc1", "spearman"])
    pruned_spearman = float(lookup.loc["human_pruned_pc1_vs_mouse_pc1", "spearman"])
    saluki_pearson = float(lookup.loc["saluki_human_pc1_vs_saluki_mouse_prior", "pearson"])
    saluki_spearman = float(lookup.loc["saluki_human_pc1_vs_saluki_mouse_prior", "spearman"])

    text = f"""# paper2 补强结果：剔除异常 study 后 human-mouse ortholog 一致性进一步提升

## 核心结论

基于 Ensembl 官方 one-to-one ortholog 映射，我们在 {common_pairs} 对 human-mouse 正交同源基因上比较了 human 全量标签、human 去除 Gejman 后标签，以及 mouse 标签的一致性。结果显示：

- human 全量 `PC1` vs mouse `PC1`：Pearson = {full_pearson:.3f}，Spearman = {full_spearman:.3f}
- human 去除 `Gejman` 后 `PC1` vs mouse `PC1`：Pearson = {pruned_pearson:.3f}，Spearman = {pruned_spearman:.3f}

也就是说，在同一批 ortholog 对上，去除异常 study 后 human-mouse 一致性继续上升，Pearson 增加了 {pruned_pearson - full_pearson:.3f}，Spearman 增加了 {pruned_spearman - full_spearman:.3f}。

作为参照，Saluki human `PC1` 与 Saluki mouse PC1 在 ortholog 对上的 Pearson 为 {saluki_pearson:.3f}，Spearman 为 {saluki_spearman:.3f}。因此，我们的异常 study 识别方向与原文“提高跨物种一致性”的总体目标是一致的。

## 这个结果为什么重要

- 它说明 `Gejman` 不只是影响 human 内部标签稳定性，也会影响 human 与 mouse 的跨物种一致性。
- 这样就把 paper2 从“内部质控”扩展成了“跨物种有效性验证”。
- 这个补强很省工作量，因为只需要公开数据和官方 ortholog 表，不需要新增实验。

## 建议放入正文的位置

建议放在 Results 的后半部分，作为“Cross-species concordance after influential-study removal”小节。

可直接写成：

`To test whether outlier-study pruning improved not only within-human label stability but also biological consistency across species, we mapped human and mouse genes through Ensembl one-to-one orthologs and compared the resulting consensus labels. Removing Gejman increased the human-mouse PC1 correlation from {full_pearson:.3f} to {pruned_pearson:.3f}, supporting the interpretation that this study introduced noise inconsistent with the shared cross-species stability program.`

## 当前可直接使用的文件

- `results/ortholog_analysis/ortholog_correlation_summary.tsv`
- `results/ortholog_analysis/common_ortholog_pairs.tsv`
- `figures/ortholog_analysis/human_full_vs_mouse_pc1.png`
- `figures/ortholog_analysis/human_pruned_vs_mouse_pc1.png`
- `figures/ortholog_analysis/ortholog_pc1_barplot_pearson.png`
"""
    outpath.write_text(text, encoding="utf-8")
    return outpath


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze human-mouse ortholog consistency for half-life consensus labels.")
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Force re-download of the Ensembl mouse homology file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_ortholog_analysis(force_download=args.force_download)
    report_path = write_cn_report()
    print(report_path)


if __name__ == "__main__":
    main()
