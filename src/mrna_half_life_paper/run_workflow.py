from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mrna_half_life_paper.benchmark import leave_one_group_out
from mrna_half_life_paper.bias_aware import BiasAwareConfig, fit_bias_aware_consensus
from mrna_half_life_paper.io import (
    align_matrix_and_metadata,
    load_matrix,
    load_metadata,
    save_table,
)
from mrna_half_life_paper.config import ROOT
from mrna_half_life_paper.pc1_consensus import compute_mean_consensus, compute_pc1_consensus
from mrna_half_life_paper.preprocess import PreprocessConfig, preprocess_matrix


def _build_estimators():
    def mean_estimator(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.Series:
        del metadata
        return compute_mean_consensus(matrix)

    def pc1_estimator(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.Series:
        del metadata
        return compute_pc1_consensus(matrix)

    def bias_estimator(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.Series:
        return fit_bias_aware_consensus(matrix, metadata, BiasAwareConfig()).gene_consensus

    return {
        "mean": mean_estimator,
        "pc1": pc1_estimator,
        "bias_aware": bias_estimator,
    }


def _portable_path(path: Path) -> str:
    """Prefer repository-relative paths in generated metadata."""
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def run_workflow(
    matrix_path: Path,
    metadata_path: Path,
    outdir: Path,
    apply_log: bool = True,
    min_observed_per_gene: int = 3,
    n_components: int | None = None,
    max_impute_iter: int = 30,
    impute_tol: float = 1e-6,
    bias_max_iter: int = 200,
    bias_tol: float = 1e-8,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    matrix = load_matrix(matrix_path)
    metadata = load_metadata(metadata_path)
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)

    preprocess_config = PreprocessConfig(
        apply_log=apply_log,
        min_observed_per_gene=min_observed_per_gene,
        n_components=n_components,
        max_iter=max_impute_iter,
        tol=impute_tol,
    )
    bias_config = BiasAwareConfig(max_iter=bias_max_iter, tol=bias_tol)

    processed = preprocess_matrix(matrix, preprocess_config)
    save_table(processed, outdir / "processed_matrix.tsv")
    metadata.to_csv(outdir / "aligned_metadata.tsv", sep="\t", index=False)

    mean_label = compute_mean_consensus(processed)
    pc1_label = compute_pc1_consensus(processed)
    bias_result = fit_bias_aware_consensus(processed, metadata, bias_config)

    labels = pd.concat(
        [
            mean_label.rename("mean"),
            pc1_label.rename("pc1"),
            bias_result.gene_consensus.rename("bias_aware"),
            bias_result.gene_uncertainty.rename("gene_uncertainty"),
        ],
        axis=1,
    )
    save_table(labels, outdir / "consensus_labels.tsv")
    save_table(bias_result.fitted_matrix, outdir / "bias_aware_fitted_matrix.tsv")
    save_table(bias_result.residual_matrix, outdir / "bias_aware_residual_matrix.tsv")
    save_table(bias_result.study_effects.to_frame("study_effect"), outdir / "study_effects.tsv")
    save_table(bias_result.method_effects.to_frame("method_effect"), outdir / "method_effects.tsv")

    estimators = _build_estimators()
    benchmark_frames = []
    for name, estimator in estimators.items():
        benchmark_frames.append(leave_one_group_out(processed, metadata, "study", name, estimator))
        benchmark_frames.append(leave_one_group_out(processed, metadata, "method", name, estimator))
    benchmark_df = pd.concat(benchmark_frames, ignore_index=True)
    benchmark_df.to_csv(outdir / "leave_one_group_out_benchmark.tsv", sep="\t", index=False)

    summary = {
        "n_genes_input": int(matrix.shape[0]),
        "n_samples_input": int(matrix.shape[1]),
        "n_genes_processed": int(processed.shape[0]),
        "n_samples_processed": int(processed.shape[1]),
        "apply_log": bool(apply_log),
        "min_observed_per_gene": int(min_observed_per_gene),
        "n_components": None if n_components is None else int(n_components),
        "max_impute_iter": int(max_impute_iter),
        "impute_tol": float(impute_tol),
        "bias_aware_iterations": int(bias_result.iterations),
        "bias_max_iter": int(bias_max_iter),
        "bias_tol": float(bias_tol),
        "output_dir": _portable_path(outdir),
    }
    with open(outdir / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the paper 1 half-life consensus workflow.")
    parser.add_argument("--matrix", type=Path, required=True, help="Path to the gene x sample matrix.")
    parser.add_argument("--metadata", type=Path, required=True, help="Path to the sample metadata table.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--skip-log",
        action="store_true",
        help="Skip the log transform step for matrices that are already log-transformed.",
    )
    parser.add_argument(
        "--min-observed-per-gene",
        type=int,
        default=3,
        help="Minimum number of observed samples required to retain a gene.",
    )
    parser.add_argument(
        "--n-components",
        type=int,
        default=None,
        help="Number of PCA components used during iterative imputation.",
    )
    parser.add_argument(
        "--max-impute-iter",
        type=int,
        default=30,
        help="Maximum number of iterations for PCA-based imputation.",
    )
    parser.add_argument(
        "--impute-tol",
        type=float,
        default=1e-6,
        help="Convergence tolerance for PCA-based imputation.",
    )
    parser.add_argument(
        "--bias-max-iter",
        type=int,
        default=200,
        help="Maximum number of iterations for the bias-aware estimator.",
    )
    parser.add_argument(
        "--bias-tol",
        type=float,
        default=1e-8,
        help="Convergence tolerance for the bias-aware estimator.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_workflow(
        matrix_path=args.matrix,
        metadata_path=args.metadata,
        outdir=args.outdir,
        apply_log=not args.skip_log,
        min_observed_per_gene=args.min_observed_per_gene,
        n_components=args.n_components,
        max_impute_iter=args.max_impute_iter,
        impute_tol=args.impute_tol,
        bias_max_iter=args.bias_max_iter,
        bias_tol=args.bias_tol,
    )


if __name__ == "__main__":
    main()
