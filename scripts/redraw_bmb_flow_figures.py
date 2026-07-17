#!/usr/bin/env python3
"""Draw submission-safe BMB workflow figures with native matplotlib vectors."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
SUBMISSION = ROOT / "figures" / "paper"
MAIN_OUT = SUBMISSION / "main" / "Fig01_workflow"
SUPP_OUT = SUBMISSION / "supplement" / "FigS01_prediction_boundary_flow"
MAIN_DATA = SUBMISSION / "data" / "Fig01_panel_data.tsv"

BLUE = "#356AC3"
BLUE_LIGHT = "#F3F7FE"
TEAL = "#178B8D"
TEAL_LIGHT = "#F1FAFA"
ORANGE = "#D9771E"
ORANGE_LIGHT = "#FFF8F0"
PURPLE = "#7257B6"
PURPLE_LIGHT = "#F7F4FC"
INK = "#202936"
MUTED = "#586575"
GRID = "#CBD4DF"
WHITE = "#FFFFFF"


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.2,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def _box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str = "",
    *,
    edge: str = BLUE,
    face: str = WHITE,
    title_size: float = 8.0,
    body_size: float = 6.7,
    title_color: str = INK,
    linewidth: float = 1.0,
) -> None:
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=linewidth,
        edgecolor=edge,
        facecolor=face,
        zorder=1,
    )
    ax.add_patch(patch)
    title_y = y + height * (0.60 if body else 0.50)
    ax.text(
        x + width / 2,
        title_y,
        title,
        ha="center",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=title_color,
        linespacing=1.05,
        zorder=2,
    )
    if body:
        ax.text(
            x + width / 2,
            y + height * 0.27,
            body,
            ha="center",
            va="center",
            fontsize=body_size,
            color=MUTED,
            linespacing=1.08,
            zorder=2,
        )


def _arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = INK,
    linewidth: float = 1.2,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=linewidth,
            color=color,
            # Keep connector endpoints visibly clear of rounded box borders.
            shrinkA=5,
            shrinkB=5,
            zorder=3,
        )
    )


def _band(
    ax: plt.Axes,
    y: float,
    height: float,
    label: str,
    heading: str,
    *,
    edge: str,
    face: str,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (0.015, y),
            0.97,
            height,
            boxstyle="round,pad=0.006,rounding_size=0.012",
            linewidth=1.2,
            edgecolor=edge,
            facecolor=face,
            zorder=0,
        )
    )
    ax.text(
        0.035,
        y + height - 0.036,
        label,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color=WHITE,
        bbox=dict(boxstyle="round,pad=0.28", facecolor=edge, edgecolor=edge),
    )
    ax.text(
        0.066,
        y + height - 0.036,
        heading,
        ha="left",
        va="center",
        fontsize=10.2,
        fontweight="bold",
        color=edge,
    )


def _save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    svg_path = stem.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight", facecolor=WHITE)
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    if os.environ.get("BMB_SKIP_PDF") != "1":
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor=WHITE)
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor=WHITE)
    fig.savefig(
        stem.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        facecolor=WHITE,
        pil_kwargs={"compression": "tiff_lzw"},
    )


def _build_main_panel_data() -> pd.DataFrame:
    conditional_path = RESULTS / "study_influence_sensitivity" / "conditional_same_size_summary.tsv"
    conditional = pd.read_csv(conditional_path, sep="\t").set_index("metric")
    global_path = RESULTS / "tenfold_global_prior" / "summary_by_setting.tsv"
    global_summary = pd.read_csv(global_path, sep="\t").set_index("prior_mode")
    distance_path = RESULTS / "ortholog_label_distance" / "distance_summary.tsv"
    distance = pd.read_csv(distance_path, sep="\t").set_index(["comparison", "subset"])
    target_path = RESULTS / "ortholog_regularized_label_10fold_main" / "summary_by_label.tsv"
    target = pd.read_csv(target_path, sep="\t").set_index("label")

    values = [
        (
            "A",
            "study_audit",
            "saluki_delta_r",
            conditional.loc["delta_saluki_pearson", "observed_gejman"],
            "Saluki delta r",
            conditional_path,
            "Fixed 13,265-gene agreement gain after Gejman removal",
        ),
        (
            "A",
            "study_audit",
            "ortholog_delta_r",
            conditional.loc["delta_ortholog_pearson", "observed_gejman"],
            "ortholog delta r",
            conditional_path,
            "One-to-one ortholog-concordance gain after Gejman removal",
        ),
        (
            "A",
            "study_audit",
            "geometry_tail_proportion",
            conditional.loc["pc1_stability_pearson", "empirical_tail_proportion"],
            "geometry tail proportion",
            conditional_path,
            "Conditional lower-tail proportion from 500 same-size deletions",
        ),
        (
            "B",
            "human_only",
            "pearson_mean",
            global_summary.loc["none", "pearson_mean"],
            "human-only r",
            global_path,
            "Mean Pearson across three repeated 10-fold runs",
        ),
        (
            "B",
            "fixed_target",
            "genes",
            global_summary.loc["none", "n"],
            "genes",
            global_path,
            "Global fixed-target prediction universe",
        ),
        (
            "B",
            "prior_enhanced",
            "pearson_mean",
            global_summary.loc["real", "pearson_mean"],
            "real r",
            global_path,
            "Mean Pearson across three repeated 10-fold runs",
        ),
        (
            "B",
            "prior_shuffled",
            "pearson_mean",
            global_summary.loc["shuffled", "pearson_mean"],
            "shuffled r",
            global_path,
            "Mean Pearson after fold-aware gene-prior permutation",
        ),
        (
            "C",
            "target_geometry",
            "human_pearson",
            distance.loc[
                ("orthoreg_lambda0p10_vs_human_no_gejman_pc1", "all_one2one"),
                "pearson",
            ],
            "human r",
            distance_path,
            "Shrinkage target versus target-construction human PC1",
        ),
        (
            "C",
            "target_geometry",
            "shift_rmse",
            distance.loc[
                ("orthoreg_lambda0p10_vs_human_no_gejman_pc1", "all_one2one"),
                "rmse",
            ],
            "shift RMSE",
            distance_path,
            "Z-scale target displacement on mapped pairs",
        ),
        (
            "C",
            "target_geometry",
            "ortholog_pairs",
            distance.loc[
                ("orthoreg_lambda0p10_vs_human_no_gejman_pc1", "all_one2one"),
                "n",
            ],
            "ortholog pairs",
            distance_path,
            "Mapped one-to-one pairs used for geometry checks",
        ),
        (
            "C",
            "target_prediction",
            "genes",
            target.loc["orthoreg_reconstructed_mouse_pc1_0.1", "n"],
            "genes",
            target_path,
            "Shrinkage prediction universe",
        ),
    ]
    rows = []
    for section, item, metric, value, display, source, note in values:
        rows.append(
            {
                "section": section,
                "item": item,
                "metric": metric,
                "value": float(value),
                "display": display,
                "source_file": source.relative_to(ROOT).as_posix(),
                "note": note,
            }
        )
    panel_data = pd.DataFrame(rows)
    MAIN_DATA.parent.mkdir(parents=True, exist_ok=True)
    panel_data.to_csv(MAIN_DATA, sep="\t", index=False)
    return panel_data


def _metric(panel_data: pd.DataFrame, section: str, item: str, metric: str) -> float:
    selected = panel_data.loc[
        (panel_data["section"] == section)
        & (panel_data["item"] == item)
        & (panel_data["metric"] == metric),
        "value",
    ]
    if selected.shape[0] != 1:
        raise ValueError(f"Expected one Fig. 1 value for {section}/{item}/{metric}, found {selected.shape[0]}")
    return float(selected.iloc[0])


def draw_main_workflow(panel_data: pd.DataFrame) -> None:
    saluki_delta = _metric(panel_data, "A", "study_audit", "saluki_delta_r")
    ortholog_delta = _metric(panel_data, "A", "study_audit", "ortholog_delta_r")
    geometry_tail = _metric(panel_data, "A", "study_audit", "geometry_tail_proportion")
    human_only_r = _metric(panel_data, "B", "human_only", "pearson_mean")
    fixed_target_n = round(_metric(panel_data, "B", "fixed_target", "genes"))
    real_prior_r = _metric(panel_data, "B", "prior_enhanced", "pearson_mean")
    shuffled_prior_r = _metric(panel_data, "B", "prior_shuffled", "pearson_mean")
    target_human_r = _metric(panel_data, "C", "target_geometry", "human_pearson")
    target_shift_rmse = _metric(panel_data, "C", "target_geometry", "shift_rmse")
    target_pair_n = round(_metric(panel_data, "C", "target_geometry", "ortholog_pairs"))
    target_prediction_n = round(_metric(panel_data, "C", "target_prediction", "genes"))

    fig, ax = plt.subplots(figsize=(7.25, 5.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _band(ax, 0.63, 0.35, "a", "Study-aware label reconstruction and audit", edge=BLUE, face=BLUE_LIGHT)
    source_x = 0.0325
    source_width = 0.145
    _box(
        ax,
        source_x,
        0.705,
        source_width,
        0.16,
        "Source data",
        "MOESM2\nSaluki reference\nEnsembl orthologs",
        edge=BLUE,
        body_size=6.2,
    )
    steps = [
        (0.2125, 0.105, "Coverage\nfilter", "min observed\nper gene = 3"),
        (0.3525, 0.145, "Preprocess", "sample z-score\niterative PCA\nquantile normalization"),
        (0.5325, 0.115, "Consensus\nlabel", "PC1 and\nsign alignment"),
        (0.6775, 0.110, "Study audit", "leave one\nstudy out"),
    ]
    for x, width, title, body in steps:
        _box(ax, x, 0.705, width, 0.16, title, body, edge=BLUE, title_size=7.2, body_size=5.9)
    _arrow(ax, (source_x + source_width, 0.785), (steps[0][0], 0.785))
    for (left_x, left_w, _, _), (right_x, _, _, _) in zip(steps[:-1], steps[1:]):
        _arrow(ax, (left_x + left_w, 0.785), (right_x, 0.785))
    _box(
        ax,
        0.8225,
        0.695,
        0.145,
        0.18,
        "Directional evidence",
        (
            f"Saluki delta r: {saluki_delta:+.3f}\n"
            f"ortholog delta r: {ortholog_delta:+.3f}\n"
            f"geometry lower-tail: {geometry_tail:.3f}"
        ),
        edge=ORANGE,
        face=ORANGE_LIGHT,
        title_size=6.4,
        body_size=4.9,
    )
    last_step_x, last_step_width, _, _ = steps[-1]
    _arrow(ax, (last_step_x + last_step_width, 0.785), (0.8225, 0.785), color=ORANGE)

    _band(ax, 0.325, 0.275, "b", "Fixed-target prediction: Saluki human PC1", edge=TEAL, face=TEAL_LIGHT)
    _box(ax, 0.065, 0.455, 0.25, 0.075, "Human compact features", "sequence + regulatory", edge=TEAL, title_size=7.2, body_size=6.0)
    _box(ax, 0.39, 0.455, 0.22, 0.075, "Human-only model", f"N = {fixed_target_n:,} genes", edge=TEAL, title_size=7.2, body_size=6.0)
    _box(ax, 0.70, 0.455, 0.23, 0.075, "Repeated 10-fold", f"human-only r = {human_only_r:.3f}", edge=TEAL, title_size=7.2, body_size=6.0)
    _arrow(ax, (0.315, 0.493), (0.39, 0.493), color=TEAL)
    _arrow(ax, (0.61, 0.493), (0.70, 0.493), color=TEAL)

    _box(ax, 0.065, 0.35, 0.25, 0.075, "Human features\n+ mouse priors", "same fixed human target", edge=ORANGE, face=ORANGE_LIGHT, title_size=6.8, body_size=5.8)
    _box(ax, 0.39, 0.35, 0.22, 0.075, "Cross-species transfer", "external ortholog covariates", edge=ORANGE, face=ORANGE_LIGHT, title_size=7.2, body_size=6.0)
    _box(
        ax,
        0.70,
        0.35,
        0.23,
        0.075,
        "Repeated 10-fold",
        f"real {real_prior_r:.3f}; shuffled {shuffled_prior_r:.3f}",
        edge=ORANGE,
        face=ORANGE_LIGHT,
        title_size=7.2,
        body_size=6.0,
    )
    _arrow(ax, (0.315, 0.388), (0.39, 0.388), color=ORANGE)
    _arrow(ax, (0.61, 0.388), (0.70, 0.388), color=ORANGE)

    _band(ax, 0.035, 0.255, "c", "Weak ortholog-informed target shrinkage", edge=PURPLE, face=PURPLE_LIGHT)
    _box(ax, 0.055, 0.140, 0.19, 0.065, "Target human PC1", "90% no-Gejman label", edge=PURPLE, title_size=6.8, body_size=5.7)
    _box(ax, 0.055, 0.050, 0.19, 0.065, "Target mouse PC1", "10% mapped label", edge=PURPLE, title_size=6.2, body_size=5.7)
    _box(ax, 0.335, 0.0575, 0.205, 0.14, "Ortholog-informed\nshrinkage target", "human-dominant mammalian\nstability score", edge=PURPLE, title_size=6.7, body_size=5.7)
    _arrow(ax, (0.245, 0.1725), (0.335, 0.1725), color=PURPLE)
    _arrow(ax, (0.245, 0.0825), (0.335, 0.0825), color=PURPLE)
    _box(
        ax,
        0.59,
        0.0675,
        0.17,
        0.12,
        "Geometry check",
        f"N = {target_pair_n:,} pairs\nhuman r = {target_human_r:.3f}\nshift RMSE = {target_shift_rmse:.3f}",
        edge=PURPLE,
        title_size=7.0,
        body_size=5.9,
    )
    _arrow(ax, (0.54, 0.1275), (0.59, 0.1275), color=PURPLE)
    _box(
        ax,
        0.81,
        0.0675,
        0.145,
        0.12,
        "Prediction check",
        f"N = {target_prediction_n:,} genes\nhuman-only and\nprior controls",
        edge=PURPLE,
        title_size=7.0,
        body_size=5.9,
    )
    _arrow(ax, (0.76, 0.1275), (0.81, 0.1275), color=PURPLE)
    _save(fig, MAIN_OUT)
    plt.close(fig)


def _panel_frame(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    label: str,
    title: str,
    color: str,
    face: str,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.008,rounding_size=0.012",
            linewidth=1.15,
            edgecolor=color,
            facecolor=face,
        )
    )
    ax.text(x + 0.025, y + height - 0.038, label, fontsize=10, fontweight="bold", color=color, va="center")
    ax.text(x + 0.062, y + height - 0.038, title, fontsize=8.8, fontweight="bold", color=color, va="center")


def draw_prediction_boundary() -> None:
    fig, ax = plt.subplots(figsize=(7.25, 5.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    _panel_frame(ax, 0.02, 0.52, 0.465, 0.455, "a", "Base sequence", BLUE, BLUE_LIGHT)
    _box(ax, 0.06, 0.68, 0.17, 0.14, "Available input", "5'UTR / CDS / 3'UTR\npartitioned sequence", edge=BLUE)
    _box(ax, 0.30, 0.68, 0.14, 0.14, "Base model", "526 sequence\nfeatures", edge=BLUE)
    _arrow(ax, (0.23, 0.75), (0.30, 0.75), color=BLUE)
    ax.text(
        0.252,
        0.60,
        "Applicable to new or designed sequences\nhuman-only output",
        ha="center",
        va="center",
        fontsize=6.8,
        color=MUTED,
        linespacing=1.15,
    )

    _panel_frame(ax, 0.515, 0.52, 0.465, 0.455, "b", "Compact human model", TEAL, TEAL_LIGHT)
    _box(ax, 0.555, 0.68, 0.17, 0.14, "Available input", "Ensembl human gene\n+ regulatory blocks", edge=TEAL)
    _box(ax, 0.795, 0.68, 0.14, 0.14, "Compact model", "1,802 human\nfeatures", edge=TEAL)
    _arrow(ax, (0.725, 0.75), (0.795, 0.75), color=TEAL)
    ax.text(
        0.748,
        0.60,
        "Requires precomputed CWCS, SeqWeaver\nand DeepRiPe blocks",
        ha="center",
        va="center",
        fontsize=6.8,
        color=MUTED,
        linespacing=1.15,
    )

    _panel_frame(ax, 0.02, 0.03, 0.465, 0.455, "c", "Cross-species transfer", PURPLE, PURPLE_LIGHT)
    _box(ax, 0.055, 0.16, 0.17, 0.10, "Human input", "compact features", edge=PURPLE)
    _box(ax, 0.055, 0.30, 0.17, 0.09, "Mouse ortholog priors", "values + flags", edge=PURPLE, title_size=7.0, body_size=5.9)
    _box(ax, 0.30, 0.19, 0.14, 0.175, "Transfer model", "fixed human target", edge=PURPLE)
    _arrow(ax, (0.225, 0.21), (0.30, 0.21), color=PURPLE)
    _arrow(ax, (0.225, 0.345), (0.30, 0.345), color=PURPLE)
    ax.text(
        0.252,
        0.095,
        "External mouse covariates required\nnot sequence-only prediction",
        ha="center",
        va="center",
        fontsize=6.8,
        color=MUTED,
        linespacing=1.15,
    )

    _panel_frame(ax, 0.515, 0.03, 0.465, 0.455, "d", "Coverage-aware deployment", ORANGE, ORANGE_LIGHT)
    _box(ax, 0.55, 0.315, 0.17, 0.105, "Mouse prior\navailable?", edge=ORANGE, title_size=7.2)
    _box(ax, 0.805, 0.315, 0.14, 0.105, "Use transfer\nmodel", "human + mouse", edge=ORANGE, title_size=6.7, body_size=5.7)
    _arrow(ax, (0.72, 0.368), (0.805, 0.368), color=ORANGE)
    ax.text(0.7625, 0.383, "yes", ha="center", fontsize=6.5, color=ORANGE)
    _box(ax, 0.55, 0.065, 0.17, 0.16, "Regulatory\nblocks?", edge=TEAL, title_size=7.2)
    _box(ax, 0.805, 0.16, 0.14, 0.09, "Use compact\nmodel", "human-only", edge=TEAL, title_size=6.7, body_size=5.7)
    _arrow(ax, (0.635, 0.315), (0.635, 0.225), color=ORANGE)
    ax.text(0.65, 0.267, "no", ha="left", fontsize=6.5, color=ORANGE)
    _arrow(ax, (0.72, 0.205), (0.805, 0.205), color=TEAL)
    ax.text(0.7625, 0.220, "yes", ha="center", fontsize=6.5, color=TEAL)
    _box(ax, 0.805, 0.0575, 0.14, 0.07, "Use base model", edge=BLUE, title_size=7.0)
    _arrow(ax, (0.72, 0.0925), (0.805, 0.0925), color=BLUE)
    ax.text(0.7625, 0.1075, "no", ha="center", fontsize=6.5, color=BLUE)

    _save(fig, SUPP_OUT)
    plt.close(fig)


def main() -> None:
    panel_data = _build_main_panel_data()
    draw_main_workflow(panel_data)
    draw_prediction_boundary()
    suffixes = "svg|png|tiff" if os.environ.get("BMB_SKIP_PDF") == "1" else "svg|pdf|png|tiff"
    print(f"Wrote {MAIN_OUT}.[{suffixes}]")
    print(f"Wrote {SUPP_OUT}.[{suffixes}]")


if __name__ == "__main__":
    main()
