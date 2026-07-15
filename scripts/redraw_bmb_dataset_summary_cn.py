#!/usr/bin/env python3
"""Draw a richer Chinese-main-text dataset summary figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
FIG_MAIN = ROOT / "manuscript" / "bmb_submission" / "figures" / "main"
DPI = 600

COLOR = {
    "gray": "#7A828A",
    "blue": "#4A74D9",
    "teal": "#2F8C6B",
    "green": "#5D9A6F",
    "orange": "#C9772F",
    "purple": "#7A78B8",
    "ink": "#1f2933",
    "light_grid": "#d8dee6",
    "panel_bg": "#f7f9fc",
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
    fig.savefig(out_no_ext.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_no_ext.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str, *, x: float = -0.14, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        color=COLOR["ink"],
        va="top",
        ha="left",
    )


def annotate_vertical_bars(ax: plt.Axes, dy: float) -> None:
    for patch in ax.patches:
        height = patch.get_height()
        if height <= 0:
            continue
        x = patch.get_x() + patch.get_width() / 2
        label = f"{int(round(height)):,}" if height >= 100 else f"{int(round(height))}"
        ax.text(x, height + dy, label, ha="center", va="bottom", fontsize=6.9, color=COLOR["ink"])


def annotate_horizontal_bars(ax: plt.Axes, dx: float = 120) -> None:
    for patch in ax.patches:
        width = patch.get_width()
        if width <= 0:
            continue
        y = patch.get_y() + patch.get_height() / 2
        ax.text(width + dx, y, f"{int(round(width)):,}", ha="left", va="center", fontsize=7.0, color=COLOR["ink"])


def draw_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    w: float,
    h: float,
    title: str,
    body: str,
    edge: str,
    *,
    facecolor: str = "white",
    title_fs: float = 8.0,
    body_fs: float = 6.8,
    title_y: float = 0.66,
    body_y: float = 0.32,
) -> None:
    patch = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.0,
        edgecolor=edge,
        facecolor=facecolor,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h * title_y,
        title,
        ha="center",
        va="center",
        fontsize=title_fs,
        color=COLOR["ink"],
        fontweight="bold",
    )
    ax.text(
        xy[0] + w / 2,
        xy[1] + h * body_y,
        body,
        ha="center",
        va="center",
        fontsize=body_fs,
        color=COLOR["ink"],
    )


def draw_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str) -> None:
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.2, color=color))


def panel_a(ax: plt.Axes) -> None:
    data = pd.DataFrame(
        {
            "Species": ["Human", "Mouse", "No-Gejman\nhuman"],
            "Samples": [54, 27, 39],
            "Studies": [19, 17, 18],
            "Methods": [5, 3, 5],
        }
    )
    plot = data.melt(id_vars="Species", value_vars=["Samples", "Studies", "Methods"], var_name="Metric", value_name="Count")
    palette = {"Samples": COLOR["blue"], "Studies": COLOR["orange"], "Methods": COLOR["purple"]}
    sns.barplot(
        data=plot,
        x="Species",
        y="Count",
        hue="Metric",
        palette=palette,
        saturation=0.92,
        edgecolor="white",
        linewidth=0.9,
        ax=ax,
    )
    ax.set_title("Experimental design size", fontsize=10.2, color=COLOR["ink"], pad=5)
    ax.set_xlabel("")
    ax.set_ylabel("Count")
    ax.tick_params(axis="x", labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(title="", fontsize=7.0, frameon=False, loc="upper right")
    annotate_vertical_bars(ax, dy=0.7)


def panel_b(ax: plt.Axes) -> None:
    data = pd.DataFrame(
        {
            "Set": ["Human", "Mouse", "No-Gejman\nhuman"],
            "Processed genes": [16444, 16951, 15788],
        }
    )
    sns.barplot(
        data=data,
        x="Set",
        y="Processed genes",
        color=COLOR["purple"],
        saturation=0.92,
        edgecolor="white",
        linewidth=0.9,
        ax=ax,
    )
    ax.set_title("Consensus-label gene coverage", fontsize=10.2, color=COLOR["ink"], pad=5)
    ax.set_xlabel("")
    ax.set_ylabel("Processed genes", labelpad=3.0, fontsize=8.4)
    ax.tick_params(axis="x", labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)
    annotate_vertical_bars(ax, dy=140)


def panel_c(ax: plt.Axes) -> None:
    data = pd.DataFrame(
        [
            ("MOESM3 control\ngenes", 13532, "Gene universe"),
            ("Global prior\ncommon genes", 12916, "Gene universe"),
            ("Ortholog validation\npairs", 12592, "Ortholog pairs"),
            ("0.10 regularized target\ngenes", 12307, "Gene universe"),
            ("0.10 target ortholog\npairs", 10768, "Ortholog pairs"),
        ],
        columns=["Universe", "Count", "Type"],
    )
    palette = {"Gene universe": COLOR["purple"], "Ortholog pairs": COLOR["orange"]}

    colors = [palette[label] for label in data["Type"]]
    y_positions = list(range(len(data)))
    bars = ax.barh(
        y_positions,
        data["Count"],
        color=colors,
        edgecolor="white",
        linewidth=0.9,
        height=0.72,
        alpha=0.92,
    )

    ax.set_title("Downstream analysis universes", fontsize=10.2, color=COLOR["ink"], pad=5)
    ax.set_xlabel("Items retained downstream")
    ax.set_ylabel("")
    ax.set_yticks(y_positions, data["Universe"])
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=7.1, pad=2)
    ax.set_xlim(0, 16200)
    ax.grid(axis="x", color=COLOR["light_grid"], linewidth=0.6)
    ax.grid(axis="y", visible=False)

    short_type = {"Gene universe": "genes", "Ortholog pairs": "pairs"}
    for bar, count, item_type in zip(bars, data["Count"], data["Type"]):
        y = bar.get_y() + bar.get_height() / 2
        ax.text(count + 70, y, f"{count:,}", ha="left", va="center", fontsize=7.4, color=COLOR["ink"])
        ax.text(
            430,
            y,
            short_type[item_type],
            ha="left",
            va="center",
            fontsize=6.9,
            color=palette[item_type],
            bbox=dict(
                boxstyle="round,pad=0.18,rounding_size=0.12",
                fc="white",
                ec=palette[item_type],
                lw=0.8,
            ),
        )


def panel_d(ax: plt.Axes) -> None:
    ax.set_title("Data source roles", fontsize=10.5, color=COLOR["ink"], pad=5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.0)
    ax.axis("off")
    ax.set_facecolor("white")

    lanes = [
        {
            "y": 0.79,
            "color": COLOR["blue"],
            "source_title": "MOESM2",
            "source_body": "transformed matrix",
            "role_title": "Study audit",
            "role_body": "PC1 rebuild + audit",
        },
        {
            "y": 0.50,
            "color": COLOR["green"],
            "source_title": "MOESM3",
            "source_body": "Saluki labels + matrix",
            "role_title": "Reference / control",
            "role_body": "Saluki target + controls",
        },
        {
            "y": 0.21,
            "color": COLOR["orange"],
            "source_title": "1:1 orthologs",
            "source_body": "Ensembl pairs",
            "role_title": "Ortholog checks",
            "role_body": "concordance +\nregularization",
        },
    ]

    for lane in lanes:
        draw_box(
            ax,
            (0.05, lane["y"] - 0.085),
            0.35,
            0.17,
            lane["source_title"],
            lane["source_body"],
            lane["color"],
            facecolor=mpl.colors.to_rgba(lane["color"], 0.05),
            title_fs=10.0,
            body_fs=8.0,
            title_y=0.68,
            body_y=0.27,
        )
        draw_box(
            ax,
            (0.54, lane["y"] - 0.085),
            0.35,
            0.17,
            lane["role_title"],
            lane["role_body"],
            lane["color"],
            facecolor=mpl.colors.to_rgba(lane["color"], 0.05),
            title_fs=9.2,
            body_fs=7.6,
            title_y=0.68,
            body_y=0.27,
        )
        ax.text(0.465, lane["y"], "→", fontsize=18, color=lane["color"], ha="center", va="center")


def main() -> None:
    setup_style()
    sns.set_style("whitegrid", {"grid.color": COLOR["light_grid"], "grid.linewidth": 0.6})

    fig = plt.figure(figsize=(8.3, 7.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[0.98, 1.22])
    top = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.24)
    bottom = outer[1].subgridspec(1, 2, width_ratios=[0.92, 1.08], wspace=0.16)
    axes = [
        fig.add_subplot(top[0, 0]),
        fig.add_subplot(top[0, 1]),
        fig.add_subplot(bottom[0, 0]),
        fig.add_subplot(bottom[0, 1]),
    ]

    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    panel_d(axes[3])

    for label, ax in zip(["a", "b", "c"], axes[:3]):
        add_panel_label(ax, label)
    add_panel_label(axes[3], "d", x=-0.06, y=1.08)

    fig.tight_layout(w_pad=2.2, h_pad=1.95)
    out = FIG_MAIN / "Fig02_dataset_summary_cn"
    save_all(fig, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
