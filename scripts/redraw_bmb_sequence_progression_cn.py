#!/usr/bin/env python3
"""Draw a richer Chinese-main-text human-only benchmark summary figure.

This figure replaces the separate progression and ablation plots in the
Chinese manuscript with a 2x2 panel:
A) algorithm and feature comparison
B) regulatory-block comparison
C) label-definition sensitivity
D) subset sensitivity
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIG_MAIN = ROOT / "figures" / "paper" / "main"
DPI = 600

COLOR = {
    "gray": "#7A828A",
    "teal": "#2F8C6B",
    "blue": "#4A74D9",
    "purple": "#7A78B8",
    "orange": "#C9772F",
    "green": "#5D9A6F",
    "red": "#C45A5A",
    "ink": "#1f2933",
    "light_grid": "#d8dee6",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 8.8,
            "axes.linewidth": 0.8,
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
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
    svg_path = out_no_ext.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.14,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        color=COLOR["ink"],
        va="bottom",
        ha="left",
    )


def annotate_vertical_bars(ax: plt.Axes, fmt: str = "{:.3f}", dy: float = 0.0012) -> None:
    for patch in ax.patches:
        height = patch.get_height()
        x = patch.get_x() + patch.get_width() / 2
        ax.text(
            x,
            height + dy,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=7.2,
            color=COLOR["ink"],
        )


def annotate_horizontal_bars(ax: plt.Axes, fmt: str = "{:.3f}", dx: float = 0.0015) -> None:
    for patch in ax.patches:
        width = patch.get_width()
        y = patch.get_y() + patch.get_height() / 2
        ax.text(
            width + dx,
            y,
            fmt.format(width),
            ha="left",
            va="center",
            fontsize=7.2,
            color=COLOR["ink"],
        )


def draw_vertical_points(
    ax: plt.Axes,
    labels: list[str],
    values: list[float],
    colors: list[str],
    *,
    value_offset: float,
) -> None:
    x = np.arange(len(labels))
    ax.scatter(x, values, s=52, c=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_xticks(x, labels)
    for xpos, value in zip(x, values):
        ax.text(
            xpos,
            value + value_offset,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.2,
            color=COLOR["ink"],
        )


def draw_panel_a(ax: plt.Axes) -> None:
    seq = pd.read_csv(RESULTS / "human_sequence_prediction_common_benchmark.tsv", sep="\t")
    compact = pd.read_csv(RESULTS / "human_compact_all_common_targets.tsv", sep="\t")
    rows = [
        ("Ridge", float(seq.loc[(seq["target"] == "saluki_human_pc1") & (seq["model"] == "ridge"), "pearson"].iloc[0])),
        ("HGB", float(seq.loc[(seq["target"] == "saluki_human_pc1") & (seq["model"] == "hgb"), "pearson"].iloc[0])),
        ("XGBoost/CUDA", float(seq.loc[(seq["target"] == "saluki_human_pc1") & (seq["model"] == "xgb_cuda"), "pearson"].iloc[0])),
        ("Compact all", float(compact.loc[(compact["target"] == "saluki_human_pc1") & (compact["model"] == "compact_all_xgb_cuda"), "pearson"].iloc[0])),
    ]
    df = pd.DataFrame(rows, columns=["Model", "Pearson"])
    palette = [COLOR["gray"], COLOR["gray"], COLOR["teal"], COLOR["blue"]]

    draw_vertical_points(
        ax,
        df["Model"].tolist(),
        df["Pearson"].astype(float).tolist(),
        palette,
        value_offset=0.0022,
    )
    ax.set_title(
        "Algorithm and feature comparison\nN = 12,916; matched 5-fold",
        fontsize=9.4,
        color=COLOR["ink"],
        pad=6,
        linespacing=1.25,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.67, 0.75)
    ax.tick_params(axis="x", rotation=18, labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)


def draw_panel_b(ax: plt.Axes) -> None:
    df = pd.read_csv(RESULTS / "human_compact_regulatory_benchmark.tsv", sep="\t")
    order = ["compact_base", "compact_cwcs", "compact_seqweaver", "compact_deepripe", "compact_all"]
    label_map = {
        "compact_base": "Base",
        "compact_cwcs": "+ CWCS",
        "compact_seqweaver": "+ SeqWeaver",
        "compact_deepripe": "+ DeepRiPe",
        "compact_all": "Compact all",
    }
    keep = df.loc[df["setting"].isin(order)].copy()
    keep["setting"] = pd.Categorical(keep["setting"], categories=order, ordered=True)
    keep = keep.sort_values("setting")
    keep["label"] = keep["setting"].map(label_map)

    palette = [COLOR["gray"], COLOR["teal"], COLOR["purple"], COLOR["green"], COLOR["blue"]]
    draw_vertical_points(
        ax,
        keep["label"].astype(str).tolist(),
        keep["pearson"].astype(float).tolist(),
        palette,
        value_offset=0.0008,
    )
    ax.set_title(
        "Regulatory-block comparison\nN = 13,532; 5-fold",
        fontsize=9.4,
        color=COLOR["ink"],
        pad=6,
        linespacing=1.25,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.705, 0.73)
    ax.tick_params(axis="x", rotation=18, labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)


def draw_panel_c(ax: plt.Axes) -> None:
    df = pd.read_csv(RESULTS / "human_compact_all_common_targets.tsv", sep="\t")
    order = ["human_full_pc1", "human_pruned_pc1", "saluki_human_pc1"]
    label_map = {
        "human_full_pc1": "Full human\nPC1",
        "human_pruned_pc1": "No-Gejman\nPC1",
        "saluki_human_pc1": "Saluki\nhuman PC1",
    }
    keep = df.loc[df["target"].isin(order)].copy()
    keep["target"] = pd.Categorical(keep["target"], categories=order, ordered=True)
    keep = keep.sort_values("target")
    keep["label"] = keep["target"].map(label_map)
    palette = [COLOR["gray"], COLOR["blue"], COLOR["green"]]

    draw_vertical_points(
        ax,
        keep["label"].astype(str).tolist(),
        keep["pearson"].astype(float).tolist(),
        palette,
        value_offset=0.0015,
    )
    ax.set_title(
        "Label-definition sensitivity\nN = 12,916; matched 5-fold",
        fontsize=9.4,
        color=COLOR["ink"],
        pad=6,
        linespacing=1.25,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.70, 0.75)
    ax.tick_params(axis="x", rotation=0, labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)


def draw_panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(RESULTS / "human_compact_all_subset_scan.tsv", sep="\t")
    order = [
        "all_common",
        "ortholog_all",
        "ortholog_highconf",
        "reg_top33",
    ]
    label_map = {
        "all_common": "All common",
        "ortholog_all": "All orthologs",
        "ortholog_highconf": "High-conf orthologs",
        "reg_top33": "Top regulatory 33%",
    }
    keep = df.loc[df["subset"].isin(order)].copy()
    keep["subset"] = pd.Categorical(keep["subset"], categories=order, ordered=True)
    keep = keep.sort_values("subset", ascending=False)
    keep["label"] = keep["subset"].map(label_map)

    y = np.arange(len(keep))
    values = keep["pearson"].astype(float).to_numpy()
    ax.scatter(values, y, s=52, color=COLOR["purple"], edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_yticks(y, keep["label"].astype(str).tolist())
    ax.set_title(
        "Gene-subset sensitivity\n3-fold",
        fontsize=9.4,
        color=COLOR["ink"],
        pad=6,
        linespacing=1.25,
    )
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("")
    ax.set_xlim(0.70, 0.745)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=7.6)
    for xpos, ypos in zip(values, y):
        ax.text(
            xpos + 0.0013,
            ypos,
            f"{xpos:.3f}",
            ha="left",
            va="center",
            fontsize=7.2,
            color=COLOR["ink"],
        )


def main() -> None:
    setup_style()
    sns.set_style("whitegrid", {"grid.color": COLOR["light_grid"], "grid.linewidth": 0.6})

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 6.35))
    draw_panel_a(axes[0, 0])
    draw_panel_b(axes[0, 1])
    draw_panel_c(axes[1, 0])
    draw_panel_d(axes[1, 1])

    for label, ax in zip(["a", "b", "c", "d"], axes.flatten()):
        add_panel_label(ax, label)

    fig.tight_layout(w_pad=1.5, h_pad=3.2)
    out = FIG_MAIN / "Fig06_sequence_progression_cn"
    save_all(fig, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
