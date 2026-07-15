#!/usr/bin/env python3
"""Regenerate the three quantitative supplementary diagnostic figures."""

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
        axes[0].hlines(i, value, 1.0, color=colors[i], lw=1.8)
        axes[0].scatter(value, i, s=38, color=colors[i], zorder=3)
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
    axes[1].set_title("Agreement change vs Saluki mouse prior")
    axes[1].set_xlabel("Delta Pearson after removing one study")
    axes[1].set_ylabel("Study removed")
    axes[1].grid(axis="x", color=GRID, lw=0.65)
    axes[1].grid(axis="y", visible=False)

    panel_label(axes[0], "a", x=-0.16)
    panel_label(axes[1], "b", x=-0.18)
    fig.tight_layout(w_pad=2.2)
    save_all(fig, OUT / "FigS_mouse_study_influence")


def draw_tuning_scan() -> None:
    frame = pd.read_csv(RESULTS / "human_compact_all_tuning_scan.tsv", sep="\t").sort_values("pearson")
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    labels = [f"Config {int(value)}" for value in frame["config_id"]]
    colors = [ORANGE if int(value) == 6 else BLUE for value in frame["config_id"]]
    bars = ax.barh(labels, frame["pearson"], color=colors, height=0.66)
    xmin = float(frame["pearson"].min()) - 0.0010
    xmax = float(frame["pearson"].max()) + 0.0014
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("Parameter configuration")
    ax.set_title("Compact human-only XGBoost tuning scan", pad=8, color=INK)
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.grid(axis="y", visible=False)
    for bar, value in zip(bars, frame["pearson"]):
        ax.text(float(value) + 0.00012, bar.get_y() + bar.get_height() / 2, f"{value:.5f}",
                va="center", ha="left", fontsize=8.3, color=INK)
    ax.text(0.99, 0.02, "Orange: fixed main-setting configuration", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.2, color=ORANGE)
    fig.tight_layout()
    save_all(fig, OUT / "FigS_tuning_scan")


def main() -> None:
    setup_style()
    draw_mouse_pca()
    draw_mouse_influence()
    draw_tuning_scan()
    print(f"Wrote supplementary diagnostics to {OUT}")


if __name__ == "__main__":
    main()
