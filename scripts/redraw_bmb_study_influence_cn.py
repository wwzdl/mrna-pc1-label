#!/usr/bin/env python3
"""Redraw the main-text human study-influence figure with a cleaner layout."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as pe


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIG_MAIN = ROOT / "figures" / "paper" / "main"
DPI = 600

COLOR = {
    "blue": "#315fb0",
    "blue_light": "#dbe7ff",
    "orange": "#ef7d00",
    "orange_light": "#fff1df",
    "gray": "#7b8794",
    "grid": "#d8dee6",
    "ink": "#1f2933",
}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 10.0,
            "axes.linewidth": 0.8,
            "figure.dpi": DPI,
            "savefig.dpi": DPI,
        }
    )


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=13.0,
        fontweight="bold",
        color=COLOR["ink"],
        va="top",
        ha="left",
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


def load_data() -> pd.DataFrame:
    df = pd.read_csv(RESULTS / "study_influence" / "human" / "study_influence.tsv", sep="\t")
    return df


def draw_panel_a(ax: plt.Axes, df: pd.DataFrame) -> None:
    keep = df.sort_values("pc1_stability_pearson", ascending=True).head(8).copy()
    keep["color"] = [COLOR["orange"] if s == "Gejman" else COLOR["blue"] for s in keep["study"]]
    y = np.arange(len(keep))
    label_overrides = {
        "Gejman": {"dx": 0.0022, "dy": 0.24, "ha": "left"},
        "Rinn": {"dx": -0.0042, "dy": -0.18, "ha": "right"},
        "Dieterich": {"dx": -0.0040, "dy": 0.12, "ha": "right"},
        "Rissland2": {"dx": -0.0042, "dy": -0.18, "ha": "right"},
    }

    for i, (_, row) in enumerate(keep.iterrows()):
        ax.scatter(row["pc1_stability_pearson"], i, s=58, color=row["color"], zorder=3)
        value = row["pc1_stability_pearson"]
        override = label_overrides.get(row["study"])
        if override is not None:
            x_text = value + override["dx"]
            y_text = i + override["dy"]
            ha = override["ha"]
        elif value >= 0.982:
            x_text = value - 0.0030
            y_text = i + (0.15 if i % 2 == 0 else -0.15)
            ha = "right"
        else:
            x_text = value + 0.0022
            y_text = i + (0.15 if i % 2 == 0 else -0.15)
            ha = "left"
        label = ax.text(
            x_text,
            y_text,
            f"{value:.3f}",
            ha=ha,
            va="center",
            fontsize=8.6,
            color=row["color"],
            zorder=4,
            bbox=dict(boxstyle="round,pad=0.14,rounding_size=0.08", facecolor="white", edgecolor="none", alpha=0.92),
        )
        label.set_path_effects([pe.withStroke(linewidth=1.0, foreground="white")])

    ax.axvline(1.0, color=COLOR["gray"], linestyle=(0, (4, 3)), linewidth=1.0)
    ax.set_title("Full-label geometry", fontsize=11.8, color=COLOR["ink"], pad=6)
    ax.set_xlabel("Correlation with full human PC1", fontsize=10.5)
    ax.set_ylabel("Study removed", fontsize=10.5)
    ax.set_xlim(0.952, 1.003)
    ax.set_yticks(y)
    ax.set_yticklabels(keep["study"], fontsize=9.2)
    ax.tick_params(axis="x", labelsize=9.2)
    ax.invert_yaxis()
    ax.grid(axis="x", color=COLOR["grid"], linewidth=0.6)
    ax.grid(axis="y", visible=False)


def draw_panel_b(ax: plt.Axes, df: pd.DataFrame) -> None:
    keep = df.sort_values("delta_saluki_reference_pc1_pearson", ascending=False).head(8).copy()
    colors = [COLOR["orange"] if s == "Gejman" else COLOR["blue"] for s in keep["study"]]
    y = np.arange(len(keep))
    ax.barh(y, keep["delta_saluki_reference_pc1_pearson"], color=colors, height=0.62)
    ax.axvline(0, color=COLOR["gray"], linestyle=(0, (4, 3)), linewidth=1.0)
    ax.set_title("Agreement gain vs Saluki human PC1", fontsize=11.0, color=COLOR["ink"], pad=8, x=0.58)
    ax.set_xlabel("Delta Pearson after removing one study", fontsize=10.5)
    ax.set_ylabel("Study removed", fontsize=10.5)
    ax.set_yticks(y)
    ax.set_yticklabels(keep["study"], fontsize=9.2)
    ax.tick_params(axis="x", labelsize=9.2)
    ax.invert_yaxis()
    ax.grid(axis="x", color=COLOR["grid"], linewidth=0.6)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(-0.02, 0.04)
    ax.set_xticks(np.arange(-0.02, 0.041, 0.01))

    for i, (_, row) in enumerate(keep.iterrows()):
        value = row["delta_saluki_reference_pc1_pearson"]
        if row["study"] == "Gejman":
            ha = "right"
            x = value - 0.0038
            y_text = i + 0.74
        else:
            ha = "right"
            x = value - 0.0018
            y_text = i
        ax.text(
            x,
            y_text,
            f"{value:.3f}",
            ha=ha,
            va="center",
            fontsize=8.4,
            color=colors[i],
            bbox=dict(boxstyle="round,pad=0.10,rounding_size=0.08", facecolor="white", edgecolor="none", alpha=0.9),
        )

def main() -> None:
    setup_style()
    df = load_data()

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0))
    draw_panel_a(axes[0], df)
    draw_panel_b(axes[1], df)
    add_panel_label(axes[0], "a")
    add_panel_label(axes[1], "b", x=-0.20, y=1.07)

    fig.tight_layout(w_pad=2.2)
    out = FIG_MAIN / "Fig04_human_study_influence"
    save_all(fig, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
