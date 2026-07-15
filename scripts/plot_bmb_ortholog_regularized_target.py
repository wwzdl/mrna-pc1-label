#!/usr/bin/env python3
"""Plot the BMB main-text figure for ortholog-informed target shrinkage.

The four panels establish target geometry, the scale of the perturbation, and
cross-target preservation under human-only model inputs.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT_DIR = ROOT / "manuscript" / "bmb_submission" / "figures" / "main"
OUT_BASENAME = "Fig08_ortholog_regularized_target_cn"

PALETTE = {
    "human": "#0F4D92",
    "mouse": "#B64342",
    "teal": "#42949E",
    "rmse": "#3775BA",
    "mae": "#D98C3E",
    "neutral_dark": "#4D4D4D",
    "neutral_mid": "#767676",
    "neutral_light": "#D8D8D8",
}


def apply_publication_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.size"] = 8
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.linewidth"] = 0.8
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["figure.dpi"] = 600
    plt.rcParams["savefig.dpi"] = 600


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.07,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=PALETTE["neutral_dark"],
    )


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    values = pd.read_csv(RESULTS / "ortholog_label_distance" / "ortholog_label_values.tsv", sep="\t")
    distance = pd.read_csv(RESULTS / "ortholog_label_distance" / "distance_summary.tsv", sep="\t")
    noise = pd.read_csv(RESULTS / "study_noise_vs_orthoreg_shift" / "noise_summary.tsv", sep="\t")
    cross_target = pd.read_csv(
        RESULTS / "ortholog_regularized_cross_target" / "summary_aggregate.tsv",
        sep="\t",
    )
    cross_target_bootstrap = pd.read_csv(
        RESULTS / "ortholog_regularized_cross_target" / "paired_bootstrap_ci.tsv",
        sep="\t",
    )
    return values, distance, noise, cross_target, cross_target_bootstrap


def pick_distance(distance: pd.DataFrame, comparison: str, subset: str = "all_one2one") -> pd.Series:
    return distance.loc[
        (distance["comparison"] == comparison) & (distance["subset"] == subset)
    ].iloc[0]


def pick_noise(
    noise: pd.DataFrame,
    source: str,
    comparison_type: str,
    cohort: str,
    metric: str,
) -> pd.Series:
    return noise.loc[
        (noise["source"] == source)
        & (noise["comparison_type"] == comparison_type)
        & (noise["cohort"] == cohort)
        & (noise["metric"] == metric)
    ].iloc[0]


def scatter_panel(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    *,
        color: str,
        title: str,
    xlabel: str,
    ylabel: str,
    annotation_lines: list[str],
    lim: float,
) -> None:
    ax.hexbin(
        x.to_numpy(),
        y.to_numpy(),
        gridsize=55,
        cmap="Blues" if color == PALETTE["human"] else "OrRd",
        mincnt=1,
        linewidths=0,
        bins="log",
    )
    ax.plot([-lim, lim], [-lim, lim], color=PALETTE["neutral_mid"], linestyle="--", linewidth=1.0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, pad=5, fontsize=8.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, color="#E8E8E8", linewidth=0.6)
    ax.text(
        0.03,
        0.97,
        "\n".join(annotation_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#CFCFCF", "alpha": 0.95},
    )


def build_figure(
    values: pd.DataFrame,
    distance: pd.DataFrame,
    noise: pd.DataFrame,
    cross_target: pd.DataFrame,
    cross_target_bootstrap: pd.DataFrame,
) -> plt.Figure:
    human_row = pick_distance(distance, "orthoreg_lambda0p10_vs_human_no_gejman_pc1")
    mouse_row = pick_distance(distance, "orthoreg_lambda0p10_vs_reconstructed_mouse_pc1")
    shift_pearson = human_row["pearson"]
    shift_rmse = human_row["rmse"]
    shift_mae = human_row["mae"]

    pair_pearson = pick_noise(noise, "study_mean_labels", "study_mean_label_pair", "no_gejman", "pearson")
    pair_rmse = pick_noise(noise, "study_mean_labels", "study_mean_label_pair", "no_gejman", "rmse")
    pair_mae = pick_noise(noise, "study_mean_labels", "study_mean_label_pair", "no_gejman", "mae")

    ref_pearson = pick_noise(noise, "study_residuals", "study_label_vs_human_no_gejman_pc1", "no_gejman", "pearson")
    ref_rmse = pick_noise(noise, "study_residuals", "study_label_vs_human_no_gejman_pc1", "no_gejman", "rmse")
    ref_mae = pick_noise(noise, "study_residuals", "study_label_vs_human_no_gejman_pc1", "no_gejman", "mae")

    abs_lim = np.nanquantile(
        np.abs(
            values[
                [
                    "orthoreg_lambda0p10_z",
                    "human_no_gejman_pc1_z",
                    "reconstructed_mouse_pc1_z",
                ]
            ].to_numpy()
        ),
        0.995,
    )
    lim = max(2.6, float(abs_lim) * 1.06)

    fig = plt.figure(figsize=(7.15, 5.85))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.86], hspace=0.34, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    scatter_panel(
        ax_a,
        values["human_no_gejman_pc1_z"],
        values["orthoreg_lambda0p10_z"],
        color=PALETTE["human"],
        title="Ortholog-informed shrinkage target vs\nhuman no-Gejman PC1",
        xlabel="Human no-Gejman PC1 (z-score)",
        ylabel="0.10 shrinkage target (z-score)",
        annotation_lines=[
            f"N = {int(human_row['n']):,}",
            f"Pearson = {human_row['pearson']:.3f}",
            f"RMSE = {human_row['rmse']:.3f}",
            f"MAE = {human_row['mae']:.3f}",
        ],
        lim=lim,
    )
    add_panel_label(ax_a, "a")

    scatter_panel(
        ax_b,
        values["reconstructed_mouse_pc1_z"],
        values["orthoreg_lambda0p10_z"],
        color=PALETTE["mouse"],
        title="Ortholog-informed shrinkage target vs\nmouse PC1",
        xlabel="Mouse PC1 (z-score)",
        ylabel="0.10 shrinkage target (z-score)",
        annotation_lines=[
            f"N = {int(mouse_row['n']):,}",
            f"Pearson = {mouse_row['pearson']:.3f}",
            f"RMSE = {mouse_row['rmse']:.3f}",
            f"MAE = {mouse_row['mae']:.3f}",
        ],
        lim=lim,
    )
    add_panel_label(ax_b, "b")

    categories = [
        "Shrinkage-human\nlabel shift",
        "Between-study\npairs",
        "Study vs human\nno-Gejman PC1",
    ]
    rmse_values = [
        float(shift_rmse),
        float(pair_rmse["median"]),
        float(ref_rmse["median"]),
    ]
    mae_values = [
        float(shift_mae),
        float(pair_mae["median"]),
        float(ref_mae["median"]),
    ]
    pearson_labels = [
        f"r = {shift_pearson:.3f}",
        f"median r = {float(pair_pearson['median']):.3f}",
        f"median r = {float(ref_pearson['median']):.3f}",
    ]

    x = np.arange(len(categories))
    width = 0.34
    bars_rmse = ax_c.bar(x - width / 2, rmse_values, width=width, color=PALETTE["rmse"], label="RMSE")
    bars_mae = ax_c.bar(x + width / 2, mae_values, width=width, color=PALETTE["mae"], label="MAE")
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(categories, fontsize=6.4)
    ax_c.set_ylabel("Error in z-scored label units")
    ax_c.set_title("Label shift is smaller than study-level noise", pad=6, fontsize=8.2)
    ax_c.grid(True, axis="y", color="#E8E8E8", linewidth=0.6)
    ax_c.legend(loc="upper left", ncol=2)

    ymax = max(max(rmse_values), max(mae_values))
    ax_c.set_ylim(0, ymax * 1.40)
    for idx, label in enumerate(pearson_labels):
        top = max(rmse_values[idx], mae_values[idx])
        corr_y = max(top + ymax * 0.09, ymax * 0.14)
        ax_c.text(
            idx,
            corr_y,
            label,
            ha="center",
            va="bottom",
            fontsize=6.8,
            color=PALETTE["neutral_dark"],
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "none", "alpha": 0.92},
        )

    add_panel_label(ax_c, "c")

    for bar_group in (bars_rmse, bars_mae):
        for bar in bar_group:
            height = bar.get_height()
            ax_c.text(
                bar.get_x() + bar.get_width() / 2,
                height + ymax * 0.012,
                f"{height:.3f}",
                ha="center",
                va="bottom",
                fontsize=6.1,
                color=PALETTE["neutral_dark"],
            )

    human_row = cross_target.loc[
        (cross_target["training_target"] == "human_no_gejman_pc1")
        & (cross_target["evaluation_target"] == "human_no_gejman_pc1")
    ].iloc[0]
    regularized_row = cross_target.loc[
        (cross_target["training_target"] == "orthoreg_reconstructed_mouse_pc1_0.10")
        & (cross_target["evaluation_target"] == "human_no_gejman_pc1")
    ].iloc[0]
    difference_row = cross_target_bootstrap.loc[
        (cross_target_bootstrap["analysis"] == "paired_prediction_comparison")
        & (cross_target_bootstrap["evaluation_target"] == "human_no_gejman_pc1")
    ].iloc[0]
    difference = float(difference_row["difference_a_minus_b"])
    lower = float(difference_row["difference_ci95_low"])
    upper = float(difference_row["difference_ci95_high"])
    ax_d.axvline(0.0, color=PALETTE["neutral_mid"], linestyle="--", linewidth=1.0)
    ax_d.errorbar(
        difference,
        0.0,
        xerr=[[difference - lower], [upper - difference]],
        fmt="o",
        color=PALETTE["teal"],
        ecolor=PALETTE["teal"],
        markersize=6,
        capsize=4,
        linewidth=1.6,
    )
    margin = max(abs(lower), abs(upper)) * 1.75
    ax_d.set_xlim(-margin, margin)
    ax_d.set_ylim(-0.55, 0.65)
    ax_d.set_yticks([])
    ax_d.set_xlabel("Paired Pearson difference on human target")
    ax_d.set_title("No detectable change on human target", pad=6, fontsize=8.2)
    ax_d.grid(True, axis="x", color="#E8E8E8", linewidth=0.6)
    ax_d.text(
        0.03,
        0.92,
        "Human-only inputs; N = 12,307; 10-fold x 3 seeds",
        transform=ax_d.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        color=PALETTE["neutral_dark"],
        bbox={"boxstyle": "square,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    ax_d.text(
        0.03,
        0.78,
        f"Train human target: r = {human_row['pearson_mean']:.4f} ± {human_row['pearson_sd']:.4f}\n"
        f"Train shrinkage target: r = {regularized_row['pearson_mean']:.4f} ± {regularized_row['pearson_sd']:.4f}",
        transform=ax_d.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        color=PALETTE["neutral_dark"],
        linespacing=1.18,
        bbox={"boxstyle": "square,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    ax_d.text(
        difference,
        -0.20,
        f"Delta r = {difference:+.4f}\n95% CI [{lower:+.4f}, {upper:+.4f}]",
        ha="center",
        va="top",
        fontsize=7.0,
        color=PALETTE["neutral_dark"],
        bbox={"boxstyle": "square,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.9},
    )
    add_panel_label(ax_d, "d")

    return fig


def save_outputs(fig: plt.Figure, out_base: Path) -> None:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(
        out_base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    if os.environ.get("BMB_SKIP_PDF") != "1":
        fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    apply_publication_style()
    values, distance, noise, cross_target, cross_target_bootstrap = load_inputs()
    fig = build_figure(values, distance, noise, cross_target, cross_target_bootstrap)
    save_outputs(fig, OUT_DIR / OUT_BASENAME)
    print((OUT_DIR / OUT_BASENAME).with_suffix(".png"))


if __name__ == "__main__":
    main()
