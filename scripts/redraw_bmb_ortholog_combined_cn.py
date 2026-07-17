#!/usr/bin/env python3
"""Draw the four-panel ortholog concordance figure used in the BMB manuscript."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "ortholog_analysis"
PAIRS = RESULTS / "common_ortholog_pairs.tsv"
SUMMARY = RESULTS / "ortholog_correlation_summary.tsv"
OUT = ROOT / "figures" / "paper" / "main" / "Fig05_ortholog_summary_cn"
DPI = 600

BLUE = "#4C72C9"
ORANGE = "#B9743A"
RED = "#C7544D"
INK = "#1F2933"
GRID = "#D6DEE8"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
        }
    )


def save_all(fig: plt.Figure) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT.with_suffix(".png"), dpi=DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(
        OUT.with_suffix(".tiff"),
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    if os.environ.get("BMB_SKIP_PDF") != "1":
        fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    svg_path = OUT.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor="white")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n",
        encoding="utf-8",
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )


def style_axis(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("#BFC5CC")


def scatter_panel(
    ax: plt.Axes,
    data: pd.DataFrame,
    x_col: str,
    title: str,
    r_value: float,
) -> None:
    x = data[x_col].to_numpy(dtype=float)
    y = data["mouse_pc1"].to_numpy(dtype=float)
    ax.scatter(x, y, s=9, color=BLUE, alpha=0.10, edgecolors="none", rasterized=True)
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.linspace(np.nanpercentile(x, 0.2), np.nanpercentile(x, 99.8), 200)
    ax.plot(x_line, slope * x_line + intercept, color=RED, linewidth=1.6)
    ax.set_title(f"{title}\nr = {r_value:.3f}", pad=5)
    ax.set_xlabel("Human consensus PC1")
    ax.set_ylabel("Mouse consensus PC1")
    style_axis(ax)


def main() -> None:
    setup_style()
    pairs = pd.read_csv(PAIRS, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t").set_index("comparison")

    keys = {
        "full": "human_full_pc1_vs_mouse_pc1",
        "pruned": "human_pruned_pc1_vs_mouse_pc1",
        "high_full": "highconf_human_full_pc1_vs_mouse_pc1",
        "high_pruned": "highconf_human_pruned_pc1_vs_mouse_pc1",
    }
    metrics = {
        name: {
            "pearson": float(summary.loc[key, "pearson"]),
            "spearman": float(summary.loc[key, "spearman"]),
        }
        for name, key in keys.items()
    }

    fig, axes = plt.subplots(2, 2, figsize=(7.35, 7.05), constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    scatter_panel(ax_a, pairs, "human_full_pc1", "Full human PC1", metrics["full"]["pearson"])
    scatter_panel(
        ax_b,
        pairs,
        "human_pruned_pc1",
        "No-Gejman human PC1",
        metrics["pruned"]["pearson"],
    )

    labels = ["Full\nPC1", "No-Gejman\nPC1", "High-conf\nfull", "High-conf\nno-Gejman"]
    x = np.arange(len(labels))
    width = 0.38
    pearson = [metrics[name]["pearson"] for name in ("full", "pruned", "high_full", "high_pruned")]
    spearman = [metrics[name]["spearman"] for name in ("full", "pruned", "high_full", "high_pruned")]
    bars_p = ax_c.bar(x - width / 2, pearson, width, color=BLUE, label="Pearson")
    bars_s = ax_c.bar(x + width / 2, spearman, width, color=ORANGE, label="Spearman")
    ax_c.set_title("Correlation summary", pad=5)
    ax_c.set_ylabel("Correlation")
    ax_c.set_xticks(x, labels)
    ax_c.set_ylim(0.70, 0.80)
    ax_c.legend(frameon=False, loc="upper left")
    ax_c.bar_label(bars_p, fmt="%.3f", padding=2, fontsize=8)
    ax_c.bar_label(bars_s, fmt="%.3f", padding=2, fontsize=8)
    style_axis(ax_c)

    delta_pearson = [
        metrics["pruned"]["pearson"] - metrics["full"]["pearson"],
        metrics["high_pruned"]["pearson"] - metrics["high_full"]["pearson"],
    ]
    delta_spearman = [
        metrics["pruned"]["spearman"] - metrics["full"]["spearman"],
        metrics["high_pruned"]["spearman"] - metrics["high_full"]["spearman"],
    ]
    x_d = np.arange(2)
    bars_dp = ax_d.bar(x_d - width / 2, delta_pearson, width, color=BLUE, label="Pearson")
    bars_ds = ax_d.bar(x_d + width / 2, delta_spearman, width, color=ORANGE, label="Spearman")
    ax_d.set_title("Gain after Gejman removal", pad=5)
    ax_d.set_ylabel("Delta correlation")
    ax_d.set_xticks(x_d, ["All orthologs", "High-conf orthologs"])
    ax_d.set_ylim(0, 0.0175)
    ax_d.bar_label(bars_dp, fmt="%.3f", padding=2, fontsize=8)
    ax_d.bar_label(bars_ds, fmt="%.3f", padding=2, fontsize=8)
    style_axis(ax_d)

    for ax, label in zip(axes.ravel(), "abcd", strict=True):
        panel_label(ax, label)

    save_all(fig)
    plt.close(fig)
    formats = "png|tiff|svg" if os.environ.get("BMB_SKIP_PDF") == "1" else "png|tiff|pdf|svg"
    print(f"Wrote {OUT}.[{formats}]")


if __name__ == "__main__":
    main()
