#!/usr/bin/env python3
"""Reproduce sample-level PCA imputation diagnostics.

This script rebuilds the raw and processed sample-level PCA coordinates used in
the manuscript. It is intentionally limited to the visualization diagnostic: no
Saluki PC1 labels, mouse priors, model targets, or train/test folds are used.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.preprocess import (
    PreprocessConfig,
    filter_genes_by_coverage,
    iterative_pca_impute,
    preprocess_matrix,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce raw and processed sample-level PCA diagnostics."
    )
    parser.add_argument(
        "--species",
        choices=["human", "mouse", "both"],
        default="both",
        help="Species to process.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "results" / "sample_pca_imputation_reproducibility",
        help="Output directory.",
    )
    parser.add_argument("--min-observed-per-gene", type=int, default=3)
    parser.add_argument("--raw-components", type=int, default=5)
    parser.add_argument("--raw-max-iter", type=int, default=10)
    parser.add_argument("--raw-tol", type=float, default=1e-5)
    parser.add_argument("--processed-components", type=int, default=None)
    parser.add_argument("--processed-max-iter", type=int, default=30)
    parser.add_argument("--processed-tol", type=float, default=1e-6)
    return parser.parse_args()


def load_species_inputs(
    species: str,
    min_observed_per_gene: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    matrix_path = ROOT / "data" / "processed" / "real_data" / species / "moesm2_raw_sparse_matrix.tsv"
    metadata_path = ROOT / "metadata" / f"real_{species}_samples.tsv"

    matrix = load_matrix(matrix_path)
    metadata = load_metadata(metadata_path)
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)

    n_genes_input, n_samples_input = matrix.shape
    filtered = filter_genes_by_coverage(matrix, min_observed_per_gene=min_observed_per_gene)
    shape_info = {
        "n_genes_input": int(n_genes_input),
        "n_samples_input": int(n_samples_input),
        "n_genes_after_coverage_filter": int(filtered.shape[0]),
        "n_samples_after_alignment": int(filtered.shape[1]),
    }
    return filtered, metadata, shape_info


def sample_pca(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(matrix.T.to_numpy(dtype=float))
    out = pd.DataFrame(
        {
            "sample_id": matrix.columns,
            "PC1": coords[:, 0],
            "PC2": coords[:, 1],
            "explained_variance_ratio_PC1": pca.explained_variance_ratio_[0],
            "explained_variance_ratio_PC2": pca.explained_variance_ratio_[1],
        }
    )
    return out.merge(metadata, on="sample_id", how="left")


def plot_panel(raw: pd.DataFrame, processed: pd.DataFrame, species: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), sharex=False, sharey=False)
    for ax, df, title in [
        (axes[0], raw, "Raw transformed matrix"),
        (axes[1], processed, "Processed matrix"),
    ]:
        sns.scatterplot(
            data=df,
            x="PC1",
            y="PC2",
            hue="method",
            style="study",
            s=42,
            edgecolor="none",
            ax=ax,
            legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
    fig.suptitle(f"{species.capitalize()} sample-level PCA diagnostic", y=1.02)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def reproduce_species(species: str, args: argparse.Namespace) -> dict[str, object]:
    matrix, metadata, shape_info = load_species_inputs(species, args.min_observed_per_gene)
    outdir = args.outdir / species
    outdir.mkdir(parents=True, exist_ok=True)

    raw_imputed = iterative_pca_impute(
        matrix,
        n_components=args.raw_components,
        max_iter=args.raw_max_iter,
        tol=args.raw_tol,
    )
    processed = preprocess_matrix(
        matrix,
        PreprocessConfig(
            apply_log=False,
            min_observed_per_gene=args.min_observed_per_gene,
            n_components=args.processed_components,
            max_iter=args.processed_max_iter,
            tol=args.processed_tol,
        ),
    )

    raw_coords = sample_pca(raw_imputed, metadata)
    processed_coords = sample_pca(processed, metadata)

    raw_imputed.to_csv(outdir / "raw_low_rank_imputed_matrix.tsv", sep="\t")
    processed.to_csv(outdir / "processed_matrix.tsv", sep="\t")
    raw_coords.to_csv(outdir / "sample_pca_raw.tsv", sep="\t", index=False)
    processed_coords.to_csv(outdir / "sample_pca_processed.tsv", sep="\t", index=False)
    plot_panel(raw_coords, processed_coords, species, outdir / "sample_pca_raw_vs_processed.png")

    params = {
        "species": species,
        **shape_info,
        "input_matrix": str(
            Path("paper_pca") / "data" / "processed" / "real_data" / species / "moesm2_raw_sparse_matrix.tsv"
        ),
        "metadata": str(Path("paper_pca") / "metadata" / f"real_{species}_samples.tsv"),
        "coverage_filter": {
            "min_observed_per_gene": int(args.min_observed_per_gene),
            "applied_before_raw_and_processed_views": True,
        },
        "raw_view": {
            "purpose": "visualization only",
            "sample_zscore": False,
            "quantile_normalization": False,
            "imputation": {
                "method": "iterative low-rank PCA imputation",
                "n_components": int(args.raw_components),
                "max_iter": int(args.raw_max_iter),
                "tol": float(args.raw_tol),
                "observed_entries_are_kept_fixed": True,
                "initial_missing_values": "sample/column means; all-NaN columns use 0",
            },
        },
        "processed_view": {
            "purpose": "visualization of the label-reconstruction preprocessing space",
            "apply_log": False,
            "sample_zscore": True,
            "quantile_normalization": True,
            "imputation": {
                "method": "iterative low-rank PCA imputation",
                "n_components": args.processed_components,
                "effective_default": "min(5, min(n_genes, n_samples)-1) when n_components is null",
                "max_iter": int(args.processed_max_iter),
                "tol": float(args.processed_tol),
                "observed_entries_are_kept_fixed": True,
                "initial_missing_values": "sample/column means; all-NaN columns use 0",
            },
        },
        "sample_pca": {
            "matrix_orientation": "sample x gene after transposing gene x sample matrix",
            "sklearn_pca_n_components": 2,
            "random_state": 0,
        },
        "no_label_or_model_information_used": True,
    }
    with open(outdir / "imputation_and_pca_parameters.json", "w", encoding="utf-8") as handle:
        json.dump(params, handle, indent=2, ensure_ascii=False)
    return params


def main() -> None:
    args = parse_args()
    species_list = ["human", "mouse"] if args.species == "both" else [args.species]
    args.outdir.mkdir(parents=True, exist_ok=True)
    all_params = [reproduce_species(species, args) for species in species_list]
    with open(args.outdir / "run_summary.json", "w", encoding="utf-8") as handle:
        json.dump(all_params, handle, indent=2, ensure_ascii=False)
    print(args.outdir)


if __name__ == "__main__":
    main()
