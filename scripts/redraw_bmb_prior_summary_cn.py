#!/usr/bin/env python3
"""Draw a richer Chinese-main-text prior-enhanced benchmark summary figure.

This figure replaces the separate prior benchmark, coverage, and ortholog
conditional plots in the Chinese manuscript with a 2x2 panel:
A) global prior benchmark
B) repeated 10-fold controls
C) prior availability
D) ortholog-conditioned prior gains
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIG_MAIN = ROOT / "manuscript" / "bmb_submission" / "figures" / "main"
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
    fig.savefig(out_no_ext.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_no_ext.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        color=COLOR["ink"],
        va="top",
        ha="left",
    )


def annotate_vertical_bars(ax: plt.Axes, fmt: str = "{:.3f}", dy: float = 0.0018) -> None:
    for patch in ax.patches:
        height = patch.get_height()
        if height <= 0:
            continue
        x = patch.get_x() + patch.get_width() / 2
        ax.text(
            x,
            height + dy,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=7.1,
            color=COLOR["ink"],
        )


def annotate_horizontal_counts(ax: plt.Axes, total: int) -> None:
    for patch in ax.patches:
        width = patch.get_width()
        if width <= 0:
            continue
        y = patch.get_y() + patch.get_height() / 2
        ax.text(
            width + 130,
            y,
            f"{int(width):,} ({width / total * 100:.1f}%)",
            ha="left",
            va="center",
            fontsize=7.1,
            color=COLOR["ink"],
        )


def draw_panel_a(ax: plt.Axes) -> None:
    df = pd.read_csv(RESULTS / "human_global_mouse_prior_benchmark.tsv", sep="\t")
    order = [
        "compact_all_global",
        "compact_all_plus_mouse_pc1_global",
        "compact_all_plus_saluki_mouse_prior_global",
        "compact_all_plus_both_mouse_priors_global",
    ]
    label_map = {
        "compact_all_global": "Human-only",
        "compact_all_plus_mouse_pc1_global": "+ mouse PC1",
        "compact_all_plus_saluki_mouse_prior_global": "+ Saluki prior",
        "compact_all_plus_both_mouse_priors_global": "+ both priors",
    }
    keep = df.loc[df["setting"].isin(order)].copy()
    keep["setting"] = pd.Categorical(keep["setting"], categories=order, ordered=True)
    keep = keep.sort_values("setting")
    keep["label"] = keep["setting"].map(label_map)
    palette = [COLOR["gray"], COLOR["teal"], COLOR["green"], COLOR["blue"]]

    sns.barplot(data=keep, x="label", y="pearson", palette=palette, hue="label", legend=False, ax=ax)
    ax.set_title("Global prior benchmark", fontsize=10.5, color=COLOR["ink"], pad=6)
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.72, 0.845)
    ax.tick_params(axis="x", rotation=18, labelsize=7.6)
    ax.tick_params(axis="y", labelsize=8)
    annotate_vertical_bars(ax)


def draw_panel_b(ax: plt.Axes) -> None:
    summary = pd.read_csv(RESULTS / "tenfold_global_prior" / "summary_by_setting.tsv", sep="\t")
    coverage = pd.read_csv(RESULTS / "prior_missingness_analysis" / "summary_by_coverage.tsv", sep="\t")

    def row_from_summary(setting: str) -> pd.Series:
        return summary.loc[summary["setting"] == setting].iloc[0]

    def row_from_coverage(group: str, label: str) -> pd.Series:
        return coverage.loc[
            (coverage["coverage_group"] == group) & (coverage["setting_label"] == label)
        ].iloc[0]

    rows = [
        ("Human-only", row_from_summary("compact_all_global")["pearson_mean"], row_from_summary("compact_all_global")["pearson_sd"]),
        (
            "Shuffled priors",
            row_from_summary("compact_all_plus_both_mouse_priors_shuffled_global")["pearson_mean"],
            row_from_summary("compact_all_plus_both_mouse_priors_shuffled_global")["pearson_sd"],
        ),
        (
            "Real both priors",
            row_from_summary("compact_all_plus_both_mouse_priors_global")["pearson_mean"],
            row_from_summary("compact_all_plus_both_mouse_priors_global")["pearson_sd"],
        ),
        (
            "Fallback",
            row_from_coverage("all_genes", "coverage_aware_fallback")["pearson_mean"],
            row_from_coverage("all_genes", "coverage_aware_fallback")["pearson_sd"],
        ),
    ]
    plot = pd.DataFrame(rows, columns=["label", "pearson_mean", "pearson_sd"])
    colors = [COLOR["gray"], COLOR["orange"], COLOR["blue"], COLOR["purple"]]
    x = range(len(plot))

    ax.bar(x, plot["pearson_mean"], yerr=plot["pearson_sd"], color=colors, capsize=3, width=0.72)
    ax.set_xticks(list(x), plot["label"])
    ax.set_title("Repeated 10-fold controls", fontsize=10.5, color=COLOR["ink"], pad=6)
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.72, 0.845)
    ax.tick_params(axis="x", rotation=18, labelsize=7.6)
    ax.tick_params(axis="y", labelsize=8)
    for i, row in plot.iterrows():
        ax.text(
            i,
            row["pearson_mean"] + 0.0023,
            f"{row['pearson_mean']:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.1,
            color=COLOR["ink"],
        )


def draw_panel_c(ax: plt.Axes) -> None:
    coverage = pd.read_csv(RESULTS / "prior_missingness_analysis" / "summary_by_coverage.tsv", sep="\t")

    def n_for(group: str) -> int:
        return int(coverage.loc[coverage["coverage_group"] == group, "n"].iloc[0])

    total = n_for("all_genes")
    rows = [
        ("Both priors", n_for("has_both_mouse_priors")),
        ("Only reconstructed\nmouse PC1", n_for("only_reconstructed_mouse_pc1")),
        ("Neither prior", n_for("missing_both_mouse_priors")),
    ]
    plot = pd.DataFrame(rows, columns=["label", "n"])
    colors = [COLOR["blue"], COLOR["teal"], COLOR["gray"]]
    bars = ax.barh(plot["label"], plot["n"], color=colors, height=0.62)
    ax.invert_yaxis()
    ax.set_title("Prior availability", fontsize=10.5, color=COLOR["ink"], pad=6)
    ax.set_xlabel("Genes in the common-gene universe")
    ax.set_ylabel("")
    ax.set_xlim(0, 13_900)
    ax.tick_params(axis="x", labelsize=8)
    ax.tick_params(axis="y", labelsize=7.8)
    for bar, (_, row) in zip(bars, plot.iterrows()):
        ax.text(
            bar.get_width() + 180,
            bar.get_y() + bar.get_height() / 2,
            f"{int(row['n']):,} ({row['n'] / total * 100:.1f}%)",
            ha="left",
            va="center",
            fontsize=7.2,
            color=COLOR["ink"],
        )
    ax.text(
        0.98,
        0.04,
        f"Total N = {total:,}; Saluki-only = 0",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.0,
        color=COLOR["ink"],
    )


def draw_panel_d(ax: plt.Axes) -> None:
    df = pd.read_csv(RESULTS / "human_cross_species_compact_all_benchmark.tsv", sep="\t")
    order = [
        "compact_all__ortholog_all",
        "compact_all_plus_mouse_pc1__ortholog_all",
        "compact_all_plus_saluki_mouse_prior__ortholog_all",
        "compact_all__ortholog_highconf",
        "compact_all_plus_mouse_pc1__ortholog_highconf",
        "compact_all_plus_saluki_mouse_prior__ortholog_highconf",
    ]
    subset_map = {
        "compact_all__ortholog_all": "All orthologs",
        "compact_all_plus_mouse_pc1__ortholog_all": "All orthologs",
        "compact_all_plus_saluki_mouse_prior__ortholog_all": "All orthologs",
        "compact_all__ortholog_highconf": "High-conf orthologs",
        "compact_all_plus_mouse_pc1__ortholog_highconf": "High-conf orthologs",
        "compact_all_plus_saluki_mouse_prior__ortholog_highconf": "High-conf orthologs",
    }
    model_map = {
        "compact_all__ortholog_all": "Human-only",
        "compact_all_plus_mouse_pc1__ortholog_all": "+ mouse PC1",
        "compact_all_plus_saluki_mouse_prior__ortholog_all": "+ Saluki prior",
        "compact_all__ortholog_highconf": "Human-only",
        "compact_all_plus_mouse_pc1__ortholog_highconf": "+ mouse PC1",
        "compact_all_plus_saluki_mouse_prior__ortholog_highconf": "+ Saluki prior",
    }
    keep = df.loc[df["setting"].isin(order)].copy()
    keep["setting"] = pd.Categorical(keep["setting"], categories=order, ordered=True)
    keep = keep.sort_values("setting")
    keep["Subset"] = keep["setting"].map(subset_map)
    keep["Model"] = keep["setting"].map(model_map)
    palette = [COLOR["gray"], COLOR["blue"], COLOR["green"]]

    sns.barplot(data=keep, x="Subset", y="pearson", hue="Model", palette=palette, ax=ax)
    ax.set_title("Ortholog-conditioned gains", fontsize=10.5, color=COLOR["ink"], pad=6)
    ax.set_xlabel("")
    ax.set_ylabel("Pearson r")
    ax.set_ylim(0.72, 0.845)
    ax.tick_params(axis="x", labelsize=7.8)
    ax.tick_params(axis="y", labelsize=8)
    ax.legend(
        title="",
        fontsize=6.8,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=3,
        frameon=False,
    )
    for container in ax.containers:
        ax.bar_label(container, fmt="%.3f", padding=2, fontsize=6.8, color=COLOR["ink"])


def main() -> None:
    setup_style()
    sns.set_style("whitegrid", {"grid.color": COLOR["light_grid"], "grid.linewidth": 0.6})

    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.15))
    draw_panel_a(axes[0, 0])
    draw_panel_b(axes[0, 1])
    draw_panel_c(axes[1, 0])
    draw_panel_d(axes[1, 1])

    for label, ax in zip(["a", "b", "c", "d"], axes.flatten()):
        add_panel_label(ax, label)

    fig.tight_layout(w_pad=1.2, h_pad=1.45)
    out = FIG_MAIN / "Fig08_prior_summary_cn"
    save_all(fig, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
