#!/usr/bin/env python3
"""Redraw the main-text human sample PCA panel with larger readable text."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "study_influence" / "human"
FIG_MAIN = ROOT / "manuscript" / "bmb_submission" / "figures" / "main"
DPI = 600

COLOR = {
    "ink": "#1f2933",
    "grid": "#d8dee6",
}


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 9.5,
            "axes.titlesize": 11.5,
            "axes.labelsize": 10.0,
            "xtick.labelsize": 8.8,
            "ytick.labelsize": 8.8,
            "legend.fontsize": 7.5,
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
            "axes.linewidth": 0.8,
        }
    )


def save_all(fig: plt.Figure, out_no_ext: Path) -> None:
    out_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_no_ext.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(
        out_no_ext.with_suffix(".tiff"),
        dpi=DPI,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    if os.environ.get("BMB_SKIP_PDF") != "1":
        fig.savefig(out_no_ext.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_no_ext.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def load_pca_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(RESULTS / "sample_pca_raw.tsv", sep="\t")
    processed = pd.read_csv(RESULTS / "sample_pca_processed.tsv", sep="\t")
    return raw, processed


def panel_label(ax: plt.Axes, label: str, x: float = -0.13) -> None:
    ax.text(
        x,
        1.04,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=COLOR["ink"],
    )


def main() -> None:
    setup_style()
    raw, processed = load_pca_tables()

    methods = sorted(pd.concat([raw["method"], processed["method"]]).dropna().unique().tolist())
    studies = sorted(pd.concat([raw["study"], processed["study"]]).dropna().unique().tolist())

    palette_values = sns.color_palette("Set2", n_colors=max(len(methods), 3))
    method_palette = {method: palette_values[i] for i, method in enumerate(methods)}
    marker_values = [
        "o",
        "s",
        "^",
        "v",
        "D",
        "P",
        "X",
        "*",
        "<",
        ">",
        "h",
        "H",
        "p",
        "8",
        "d",
        "1",
        "2",
        "3",
        "4",
        "+",
        "x",
    ]
    study_markers = {study: marker_values[i % len(marker_values)] for i, study in enumerate(studies)}
    line_art_markers = {"+", "x", "1", "2", "3", "4"}

    legend_rows = math.ceil(len(studies) / 4)
    fig_height = 5.45 + 0.34 * max(0, legend_rows - 3)
    fig, axes = plt.subplots(1, 2, figsize=(7.8, fig_height), sharex=False, sharey=False)

    for ax, df, title in zip(axes, [raw, processed], ["Human Raw matrix", "Human Processed matrix"]):
        for (method, study), group in df.groupby(["method", "study"], sort=False):
            marker = study_markers[study]
            scatter_kwargs = {
                "s": 52,
                "marker": marker,
                "color": method_palette[method],
                "alpha": 0.96,
                "linewidths": 0.9 if marker in line_art_markers else 0.28,
            }
            if marker not in line_art_markers:
                scatter_kwargs["edgecolors"] = "black"
            ax.scatter(group["PC1"], group["PC2"], **scatter_kwargs)

        ax.set_title(title, pad=8, color=COLOR["ink"])
        ax.axhline(0, color="#dddddd", lw=0.85)
        ax.axvline(0, color="#dddddd", lw=0.85)
        ax.tick_params(axis="both", labelsize=8.8)
        ax.grid(color=COLOR["grid"], linewidth=0.8)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    panel_label(axes[0], "a")
    panel_label(axes[1], "b")

    method_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="white",
            label=method,
            markerfacecolor=method_palette[method],
            markeredgecolor="black",
            markeredgewidth=0.4,
            markersize=6.8,
            linestyle="None",
        )
        for method in methods
    ]
    study_handles = [
        Line2D(
            [0],
            [0],
            marker=study_markers[study],
            color="black",
            label=study,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=6.4,
            linestyle="None",
        )
        for study in studies
    ]

    fig.legend(
        handles=method_handles,
        title="Method (color)",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.19),
        ncol=min(len(method_handles), 6),
        frameon=False,
        fontsize=8.8,
        title_fontsize=9.2,
        handletextpad=0.55,
        columnspacing=1.15,
    )
    fig.legend(
        handles=study_handles,
        title="Study (marker)",
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=4,
        frameon=False,
        fontsize=8.1,
        title_fontsize=8.8,
        handletextpad=0.45,
        columnspacing=1.0,
    )

    fig.tight_layout(rect=[0, 0.38, 1, 1], w_pad=2.0)
    out = FIG_MAIN / "Fig03_human_pca_panel"
    save_all(fig, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
