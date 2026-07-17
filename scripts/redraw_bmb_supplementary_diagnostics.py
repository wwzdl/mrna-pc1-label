#!/usr/bin/env python3
"""Regenerate the quantitative supplementary diagnostic figures."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "manuscript" / "bmb_submission" / "figures" / "supplement"
DPI = 600

BLUE = "#315fb0"
TEAL = "#16858f"
ORANGE = "#ef7d00"
GRAY = "#7b8794"
GRID = "#d8dee6"
INK = "#1f2933"


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 9.5,
            "axes.titlesize": 11.2,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.6,
            "ytick.labelsize": 8.6,
            "legend.fontsize": 7.7,
            "axes.linewidth": 0.8,
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
        }
    )


def save_all(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    if os.environ.get("BMB_SKIP_PDF") != "1":
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str, x: float = -0.13) -> None:
    ax.text(x, 1.04, label, transform=ax.transAxes, ha="left", va="top", fontsize=13, fontweight="bold", color=INK)


def draw_mouse_pca() -> None:
    data_dir = RESULTS / "study_influence" / "mouse"
    raw = pd.read_csv(data_dir / "sample_pca_raw.tsv", sep="\t")
    processed = pd.read_csv(data_dir / "sample_pca_processed.tsv", sep="\t")
    methods = sorted(pd.concat([raw["method"], processed["method"]]).dropna().unique())
    studies = sorted(pd.concat([raw["study"], processed["study"]]).dropna().unique())

    method_colors = dict(zip(methods, sns.color_palette("Set2", n_colors=len(methods))))
    marker_values = ["o", "s", "^", "v", "D", "P", "X", "*", "<", ">", "h", "H", "p", "8", "d", "1", "2", "3"]
    study_markers = {study: marker_values[i] for i, study in enumerate(studies)}
    line_markers = {"1", "2", "3", "4", "+", "x"}

    legend_rows = math.ceil(len(studies) / 5)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 5.7 + 0.18 * max(0, legend_rows - 3)))
    for ax, frame, title in zip(axes, (raw, processed), ("Mouse raw matrix", "Mouse processed matrix")):
        for (method, study), group in frame.groupby(["method", "study"], sort=False):
            marker = study_markers[study]
            kwargs = {
                "s": 58,
                "marker": marker,
                "color": method_colors[method],
                "alpha": 0.96,
                "linewidths": 0.9 if marker in line_markers else 0.3,
            }
            if marker not in line_markers:
                kwargs["edgecolors"] = "black"
            ax.scatter(group["PC1"], group["PC2"], **kwargs)
        ax.set_title(title, pad=7, color=INK)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.axhline(0, color="#dddddd", lw=0.8)
        ax.axvline(0, color="#dddddd", lw=0.8)
        ax.grid(color=GRID, linewidth=0.7)

    panel_label(axes[0], "a")
    panel_label(axes[1], "b")
    method_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="white", markerfacecolor=method_colors[m],
               markeredgecolor="black", markeredgewidth=0.4, markersize=6.6, label=m)
        for m in methods
    ]
    study_handles = [
        Line2D([0], [0], marker=study_markers[s], linestyle="None", color="black", markerfacecolor="white",
               markeredgecolor="black", markeredgewidth=0.5, markersize=6.1, label=s)
        for s in studies
    ]
    fig.legend(method_handles, methods, title="Method (color)", loc="lower center", bbox_to_anchor=(0.5, 0.19),
               ncol=len(methods), frameon=False, fontsize=8.8, title_fontsize=9.1)
    fig.legend(study_handles, studies, title="Study (marker)", loc="lower center", bbox_to_anchor=(0.5, 0.01),
               ncol=5, frameon=False, fontsize=7.7, title_fontsize=8.6, handletextpad=0.35, columnspacing=0.9)
    fig.tight_layout(rect=[0, 0.35, 1, 1], w_pad=2.0)
    save_all(fig, OUT / "FigS_mouse_pca_panel")


def draw_mouse_influence() -> None:
    frame = pd.read_csv(RESULTS / "study_influence" / "mouse" / "study_influence.tsv", sep="\t")
    stability = frame.sort_values("pc1_stability_pearson", ascending=True).reset_index(drop=True)
    agreement = frame.sort_values("delta_saluki_reference_pc1_pearson", ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 6.3))

    y = np.arange(len(stability))
    colors = [ORANGE if study == "Hanna" else TEAL for study in stability["study"]]
    for i, row in stability.iterrows():
        value = float(row["pc1_stability_pearson"])
        axes[0].scatter(value, i, s=42, color=colors[i], zorder=3)
    axes[0].axvline(1.0, color=GRAY, linestyle=(0, (4, 3)), lw=1.0)
    axes[0].set_xlim(float(stability["pc1_stability_pearson"].min()) - 0.004, 1.002)
    axes[0].set_yticks(y, stability["study"])
    axes[0].invert_yaxis()
    axes[0].set_title("Full-label geometry")
    axes[0].set_xlabel("Correlation with full mouse PC1")
    axes[0].set_ylabel("Study removed")
    axes[0].grid(axis="x", color=GRID, lw=0.65)
    axes[0].grid(axis="y", visible=False)

    y = np.arange(len(agreement))
    values = agreement["delta_saluki_reference_pc1_pearson"].astype(float)
    colors = [ORANGE if study == "Hanna" else BLUE for study in agreement["study"]]
    axes[1].barh(y, values, color=colors, height=0.62)
    axes[1].axvline(0, color=GRAY, linestyle=(0, (4, 3)), lw=1.0)
    margin = max(0.0015, 0.08 * (values.max() - values.min()))
    axes[1].set_xlim(values.min() - margin, values.max() + margin)
    axes[1].set_yticks(y, agreement["study"])
    axes[1].invert_yaxis()
    axes[1].set_title("Agreement change vs Saluki mouse PC1")
    axes[1].set_xlabel("Delta Pearson after removing one study")
    axes[1].set_ylabel("Study removed")
    axes[1].grid(axis="x", color=GRID, lw=0.65)
    axes[1].grid(axis="y", visible=False)

    panel_label(axes[0], "a", x=-0.16)
    panel_label(axes[1], "b", x=-0.18)
    fig.tight_layout(w_pad=2.2)
    save_all(fig, OUT / "FigS_mouse_study_influence")


def draw_study_influence_sensitivity() -> None:
    data_dir = RESULTS / "study_influence_sensitivity"
    comparators = pd.read_csv(data_dir / "conditional_same_size_deletions.tsv", sep="\t")
    comparator_summary = pd.read_csv(data_dir / "conditional_same_size_summary.tsv", sep="\t")
    sensitivity = pd.read_csv(data_dir / "pipeline_sensitivity_summary.tsv", sep="\t")
    weighting = pd.read_csv(data_dir / "study_weighted_summary.tsv", sep="\t")

    observed = comparator_summary.set_index("metric")["observed_gejman"]
    tail_proportions = comparator_summary.set_index("metric")["empirical_tail_proportion"]
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.2))

    ax = axes[0, 0]
    values = comparators["pc1_stability_pearson"].astype(float)
    ax.hist(values, bins=24, color="#aab7c7", edgecolor="white", linewidth=0.5)
    observed_stability = float(observed["pc1_stability_pearson"])
    ax.axvline(observed_stability, color=ORANGE, lw=2.0)
    ax.text(
        observed_stability - 0.0005,
        ax.get_ylim()[1] * 0.93,
        f"Gejman {observed_stability:.3f}\nempirical tail = {tail_proportions['pc1_stability_pearson']:.3f}",
        ha="right",
        va="top",
        color=ORANGE,
        fontsize=8.5,
    )
    ax.set_title("Conditional same-size deletions", color=INK)
    ax.set_xlabel("Correlation with full human PC1")
    ax.set_ylabel("Comparator deletions (n = 500)")
    ax.grid(axis="x", color=GRID, lw=0.65)
    ax.grid(axis="y", visible=False)

    ax = axes[0, 1]
    gain_columns = ["delta_saluki_pearson", "delta_ortholog_pearson"]
    gain_labels = ["Saluki agreement", "Ortholog concordance"]
    long_comparators = comparators[gain_columns].rename(columns=dict(zip(gain_columns, gain_labels))).melt(
        var_name="Comparison axis", value_name="Delta Pearson"
    )
    sns.violinplot(
        data=long_comparators,
        x="Comparison axis",
        y="Delta Pearson",
        color="#aab7c7",
        inner=None,
        cut=0,
        linewidth=0.7,
        ax=ax,
    )
    sns.boxplot(
        data=long_comparators,
        x="Comparison axis",
        y="Delta Pearson",
        width=0.20,
        color="white",
        showfliers=False,
        linewidth=0.8,
        ax=ax,
    )
    for x, (column, label) in enumerate(zip(gain_columns, gain_labels)):
        value = float(observed[column])
        ax.scatter(x, value, s=48, color=ORANGE, edgecolor="white", linewidth=0.7, zorder=5)
        ax.text(
            x,
            value + 0.006,
            f"Gejman {value:+.3f}\nempirical tail = {tail_proportions[column]:.3f}",
            ha="center",
            va="bottom",
            color=ORANGE,
            fontsize=8.2,
        )
    ax.axhline(0, color=GRAY, linestyle=(0, (4, 3)), lw=0.9)
    ax.set_ylim(-0.135, 0.060)
    ax.set_title("Saluki and ortholog comparison gains", color=INK, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("Delta Pearson after removal")
    ax.tick_params(axis="x", labelrotation=8)
    ax.grid(axis="y", color=GRID, lw=0.65)
    ax.grid(axis="x", visible=False)

    ax = axes[1, 0]
    palette = {"iterative_pca": BLUE, "sample_median": ORANGE}
    markers = {"dynamic": "o", "fixed": "s"}
    for (method, universe), frame in sensitivity.groupby(["imputation_method", "universe_mode"]):
        ax.scatter(
            frame["gejman_delta_saluki_pearson"],
            frame["gejman_delta_ortholog_pearson"],
            s=56,
            color=palette[method],
            marker=markers[universe],
            edgecolor="white",
            linewidth=0.6,
            label=f"{method.replace('_', ' ')}; {universe}",
        )
    ax.axhline(0, color=GRAY, linestyle=(0, (4, 3)), lw=0.8)
    ax.axvline(0, color=GRAY, linestyle=(0, (4, 3)), lw=0.8)
    ax.set_xlim(-0.006, 0.130)
    ax.set_ylim(-0.003, 0.092)
    ax.set_title("Pipeline sensitivity", color=INK)
    ax.set_xlabel("Delta Saluki agreement")
    ax.set_ylabel("Delta ortholog concordance")
    ax.text(
        0.03,
        0.95,
        "Gejman ranked 1/19\nin all 12 settings",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.4,
        color=INK,
    )
    ax.legend(loc="lower right", frameon=False, fontsize=7.3)
    ax.grid(color=GRID, lw=0.65)

    ax = axes[1, 1]
    weighting = weighting.copy()
    weighting["label"] = weighting["estimator"].map(
        {
            "sample_weighted": "Samples equally weighted",
            "equal_study_weighted": "Studies equally weighted",
            "study_mean_collapsed": "Study means",
        }
    )
    y = np.arange(len(weighting))
    for i, row in weighting.iterrows():
        worst = float(row["worst_stability_pearson"])
        gejman = float(row["gejman_stability_pearson"])
        ax.plot([worst, gejman], [i, i], color="#c7ced8", lw=2.0, zorder=1)
        ax.scatter(worst, i, s=42, color=BLUE, zorder=3)
        ax.scatter(gejman, i, s=52, color=ORANGE, edgecolor="white", linewidth=0.6, zorder=4)
        ax.text(
            gejman + 0.0008,
            i,
            f"rank {int(row['gejman_rank'])}",
            ha="left",
            va="center",
            fontsize=8.1,
            color=ORANGE,
        )
    ax.set_yticks(y, weighting["label"])
    ax.invert_yaxis()
    ax.set_xlim(0.950, 0.988)
    ax.set_title("Sample-count adjustment", color=INK)
    ax.set_xlabel("Leave-one-study-out PC1 stability")
    ax.set_ylabel("")
    ax.grid(axis="x", color=GRID, lw=0.65)
    ax.grid(axis="y", visible=False)
    handles = [
        Line2D([0], [0], marker="o", linestyle="None", color=BLUE, label="Most influential study"),
        Line2D([0], [0], marker="o", linestyle="None", color=ORANGE, label="Gejman"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7.4)

    for label, ax in zip("abcd", axes.ravel()):
        panel_label(ax, label, x=-0.15)
    fig.tight_layout(h_pad=2.4, w_pad=2.2)
    save_all(fig, OUT / "FigS_study_influence_sensitivity")

    panel_data = pd.concat(
        [
            comparators.assign(panel="a-b", source="conditional_same_size_deletion"),
            sensitivity.assign(panel="c", source="pipeline_sensitivity"),
            weighting.assign(panel="d", source="study_weighting"),
        ],
        ignore_index=True,
        sort=False,
    )
    data_out = ROOT / "manuscript" / "bmb_submission" / "figures" / "data"
    data_out.mkdir(parents=True, exist_ok=True)
    panel_data.to_csv(
        data_out / "FigS04_panel_data.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )


def main() -> None:
    setup_style()
    draw_mouse_pca()
    draw_mouse_influence()
    draw_study_influence_sensitivity()
    print(f"Wrote supplementary diagnostics to {OUT}")


if __name__ == "__main__":
    main()
