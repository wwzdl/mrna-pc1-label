from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np
import pandas as pd

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import RESULTS_DIR
from mrna_half_life_paper.ortholog_analysis import ONE2ONE_ORTHOLOG_PATH, extract_mouse_human_one2one
from mrna_half_life_paper.pc1_consensus import (
    compute_group_collapsed_pc1,
    compute_pc1_consensus,
    compute_study_weighted_pc1,
)
from mrna_half_life_paper.preprocess import (
    PreprocessConfig,
    filter_genes_by_coverage,
    preprocess_matrix,
    quantile_normalize,
    zscore_by_sample,
)
from mrna_half_life_paper.study_influence import _load_species_inputs


DEFAULT_RANDOM_STATE = 20260715
DEFAULT_N_REPLICATES = 500

_COMPARATOR_STATE: dict[str, object] = {}


def _build_pc1(
    matrix: pd.DataFrame,
    *,
    min_observed_per_gene: int = 3,
    n_components: int = 5,
    imputation_method: str = "iterative_pca",
) -> pd.Series:
    if imputation_method == "iterative_pca":
        processed = preprocess_matrix(
            matrix,
            PreprocessConfig(
                apply_log=False,
                min_observed_per_gene=min_observed_per_gene,
                n_components=n_components,
            ),
        )
    elif imputation_method == "sample_median":
        processed = filter_genes_by_coverage(matrix, min_observed_per_gene)
        processed = zscore_by_sample(processed)
        sample_medians = processed.median(axis=0, skipna=True).fillna(0.0)
        processed = processed.fillna(sample_medians)
        processed = quantile_normalize(processed)
    else:
        raise ValueError(f"Unknown imputation method: {imputation_method}")
    return compute_pc1_consensus(processed)


def _paired_reference_gain(
    full: pd.Series,
    candidate: pd.Series,
    reference: pd.Series,
) -> tuple[int, float, float, float]:
    common = full.index.intersection(candidate.index).intersection(reference.index)
    full_r = float(full.loc[common].corr(reference.loc[common]))
    candidate_r = float(candidate.loc[common].corr(reference.loc[common]))
    return len(common), full_r, candidate_r, candidate_r - full_r


def _paired_ortholog_gain(
    full_human: pd.Series,
    candidate_human: pd.Series,
    mouse_pc1: pd.Series,
    mapping: pd.DataFrame,
) -> tuple[int, float, float, float]:
    pairs = mapping.loc[
        mapping["human_gene_id"].isin(full_human.index)
        & mapping["human_gene_id"].isin(candidate_human.index)
        & mapping["mouse_gene_id"].isin(mouse_pc1.index),
        ["human_gene_id", "mouse_gene_id"],
    ].drop_duplicates()
    full_values = pairs["human_gene_id"].map(full_human)
    candidate_values = pairs["human_gene_id"].map(candidate_human)
    mouse_values = pairs["mouse_gene_id"].map(mouse_pc1)
    full_r = float(full_values.corr(mouse_values))
    candidate_r = float(candidate_values.corr(mouse_values))
    return len(pairs), full_r, candidate_r, candidate_r - full_r


def _comparison_row(
    *,
    study: str,
    n_samples_removed: int,
    full_pc1: pd.Series,
    candidate_pc1: pd.Series,
    saluki_human_pc1: pd.Series,
    mouse_pc1: pd.Series,
    mapping: pd.DataFrame,
) -> dict[str, object]:
    stability = label_correlation(full_pc1, candidate_pc1)
    saluki_n, saluki_full, saluki_candidate, saluki_gain = _paired_reference_gain(
        full_pc1,
        candidate_pc1,
        saluki_human_pc1,
    )
    ortholog_n, ortholog_full, ortholog_candidate, ortholog_gain = _paired_ortholog_gain(
        full_pc1,
        candidate_pc1,
        mouse_pc1,
        mapping,
    )
    return {
        "study": study,
        "n_samples_removed": n_samples_removed,
        "stability_n": int(stability["n"]),
        "pc1_stability_pearson": float(stability["pearson"]),
        "pc1_stability_spearman": float(stability["spearman"]),
        "saluki_common_n": saluki_n,
        "saluki_full_pearson": saluki_full,
        "saluki_candidate_pearson": saluki_candidate,
        "delta_saluki_pearson": saluki_gain,
        "ortholog_common_n": ortholog_n,
        "ortholog_full_pearson": ortholog_full,
        "ortholog_candidate_pearson": ortholog_candidate,
        "delta_ortholog_pearson": ortholog_gain,
    }


def _fixed_gene_universe(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    min_observed_per_gene: int,
) -> pd.Index:
    keep = matrix.notna().sum(axis=1) >= min_observed_per_gene
    for study in sorted(metadata["study"].unique()):
        samples = metadata.loc[metadata["study"] != study, "sample_id"].tolist()
        keep &= matrix.loc[:, samples].notna().sum(axis=1) >= min_observed_per_gene
    return matrix.index[keep]


def _leave_one_study_out(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    saluki_human_pc1: pd.Series,
    mouse_pc1: pd.Series,
    mapping: pd.DataFrame,
    *,
    min_observed_per_gene: int,
    n_components: int,
    universe_mode: str,
    imputation_method: str = "iterative_pca",
) -> tuple[pd.DataFrame, pd.Series]:
    if universe_mode not in {"dynamic", "fixed"}:
        raise ValueError(f"Unknown universe mode: {universe_mode}")

    working = matrix
    if universe_mode == "fixed":
        working = matrix.loc[
            _fixed_gene_universe(matrix, metadata, min_observed_per_gene)
        ].copy()

    full_pc1 = _build_pc1(
        working,
        min_observed_per_gene=min_observed_per_gene,
        n_components=n_components,
        imputation_method=imputation_method,
    )
    rows: list[dict[str, object]] = []
    for study in sorted(metadata["study"].unique()):
        keep_meta = metadata.loc[metadata["study"] != study]
        candidate = _build_pc1(
            working.loc[:, keep_meta["sample_id"]],
            min_observed_per_gene=min_observed_per_gene,
            n_components=n_components,
            imputation_method=imputation_method,
        )
        rows.append(
            _comparison_row(
                study=str(study),
                n_samples_removed=int((metadata["study"] == study).sum()),
                full_pc1=full_pc1,
                candidate_pc1=candidate,
                saluki_human_pc1=saluki_human_pc1,
                mouse_pc1=mouse_pc1,
                mapping=mapping,
            )
        )
    return pd.DataFrame(rows).sort_values("pc1_stability_pearson"), full_pc1


def _initialize_comparator_worker(state: dict[str, object]) -> None:
    global _COMPARATOR_STATE
    _COMPARATOR_STATE = state


def _random_removal_worker(task: tuple[int, tuple[str, ...]]) -> dict[str, object]:
    replicate, removed = task
    matrix = _COMPARATOR_STATE["matrix"]
    full_pc1 = _COMPARATOR_STATE["full_pc1"]
    saluki_human_pc1 = _COMPARATOR_STATE["saluki_human_pc1"]
    mouse_pc1 = _COMPARATOR_STATE["mouse_pc1"]
    mapping = _COMPARATOR_STATE["mapping"]
    min_observed = int(_COMPARATOR_STATE["min_observed_per_gene"])
    n_components = int(_COMPARATOR_STATE["n_components"])
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("Comparator worker matrix was not initialized")
    if not all(isinstance(value, pd.Series) for value in (full_pc1, saluki_human_pc1, mouse_pc1)):
        raise TypeError("Comparator worker labels were not initialized")
    if not isinstance(mapping, pd.DataFrame):
        raise TypeError("Comparator worker ortholog mapping was not initialized")

    candidate = _build_pc1(
        matrix.drop(columns=list(removed)),
        min_observed_per_gene=min_observed,
        n_components=n_components,
    )
    row = _comparison_row(
        study=f"random_{replicate:04d}",
        n_samples_removed=len(removed),
        full_pc1=full_pc1,
        candidate_pc1=candidate,
        saluki_human_pc1=saluki_human_pc1,
        mouse_pc1=mouse_pc1,
        mapping=mapping,
    )
    row["replicate"] = replicate
    row["removed_samples"] = ";".join(removed)
    return row


def _conditional_deletion_summary(
    observed: dict[str, object],
    comparators: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    definitions = (
        ("pc1_stability_pearson", "lower"),
        ("delta_saluki_pearson", "higher"),
        ("delta_ortholog_pearson", "higher"),
    )
    for metric, direction in definitions:
        values = comparators[metric].to_numpy(dtype=float)
        value = float(observed[metric])
        if direction == "lower":
            extreme = int(np.count_nonzero(values <= value))
        else:
            extreme = int(np.count_nonzero(values >= value))
        rows.append(
            {
                "metric": metric,
                "extreme_direction": direction,
                "observed_gejman": value,
                "comparator_n": len(values),
                "comparator_median": float(np.median(values)),
                "comparator_q025": float(np.quantile(values, 0.025)),
                "comparator_q975": float(np.quantile(values, 0.975)),
                "observed_percentile": float(np.mean(values <= value)),
                "empirical_tail_proportion": float((extreme + 1) / (len(values) + 1)),
            }
        )
    return pd.DataFrame(rows)


def _study_weighted_analysis(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    saluki_human_pc1: pd.Series,
    *,
    min_observed_per_gene: int,
    n_components: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed_full = preprocess_matrix(
        matrix,
        PreprocessConfig(
            apply_log=False,
            min_observed_per_gene=min_observed_per_gene,
            n_components=n_components,
        ),
    )
    labels = {
        "sample_weighted": compute_pc1_consensus(processed_full),
        "equal_study_weighted": compute_study_weighted_pc1(processed_full, metadata),
        "study_mean_collapsed": compute_group_collapsed_pc1(processed_full, metadata),
    }
    rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for estimator, full_label in labels.items():
        saluki = label_correlation(full_label, saluki_human_pc1)
        estimator_rows: list[dict[str, object]] = []
        for study in sorted(metadata["study"].unique()):
            keep_meta = metadata.loc[metadata["study"] != study]
            processed = preprocess_matrix(
                matrix.loc[:, keep_meta["sample_id"]],
                PreprocessConfig(
                    apply_log=False,
                    min_observed_per_gene=min_observed_per_gene,
                    n_components=n_components,
                ),
            )
            if estimator == "sample_weighted":
                candidate = compute_pc1_consensus(processed)
            elif estimator == "equal_study_weighted":
                candidate = compute_study_weighted_pc1(processed, keep_meta)
            else:
                candidate = compute_group_collapsed_pc1(processed, keep_meta)
            stability = label_correlation(full_label, candidate)
            row = {
                "estimator": estimator,
                "study": str(study),
                "n_samples_removed": int((metadata["study"] == study).sum()),
                "n_genes_overlap": int(stability["n"]),
                "pearson": float(stability["pearson"]),
                "spearman": float(stability["spearman"]),
            }
            rows.append(row)
            estimator_rows.append(row)
        ranked = sorted(estimator_rows, key=lambda row: float(row["pearson"]))
        gejman = next(row for row in estimator_rows if row["study"] == "Gejman")
        summary_rows.append(
            {
                "estimator": estimator,
                "min_observed_per_gene": min_observed_per_gene,
                "n_components": n_components,
                "n_genes": len(full_label),
                "saluki_pearson": float(saluki["pearson"]),
                "worst_left_out": ranked[0]["study"],
                "worst_stability_pearson": ranked[0]["pearson"],
                "gejman_rank": 1 + next(i for i, row in enumerate(ranked) if row["study"] == "Gejman"),
                "gejman_stability_pearson": gejman["pearson"],
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def run_study_influence_sensitivity(
    *,
    results_root: Path | None = None,
    n_replicates: int = DEFAULT_N_REPLICATES,
    random_state: int = DEFAULT_RANDOM_STATE,
    n_jobs: int = 8,
) -> dict[str, Path]:
    results_root = results_root or RESULTS_DIR / "study_influence_sensitivity"
    results_root.mkdir(parents=True, exist_ok=True)

    matrix, metadata, saluki_human_pc1 = _load_species_inputs("human")
    mouse_matrix, _mouse_metadata, _saluki_mouse = _load_species_inputs("mouse")
    mouse_pc1 = _build_pc1(mouse_matrix, min_observed_per_gene=3, n_components=5)
    extract_mouse_human_one2one()
    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")

    primary_dynamic, full_pc1 = _leave_one_study_out(
        matrix,
        metadata,
        saluki_human_pc1,
        mouse_pc1,
        mapping,
        min_observed_per_gene=3,
        n_components=5,
        universe_mode="dynamic",
    )
    primary_fixed, _ = _leave_one_study_out(
        matrix,
        metadata,
        saluki_human_pc1,
        mouse_pc1,
        mapping,
        min_observed_per_gene=3,
        n_components=5,
        universe_mode="fixed",
    )
    primary_dynamic.assign(universe_mode="dynamic").to_csv(
        results_root / "primary_dynamic_loo.tsv", sep="\t", index=False
    )
    primary_fixed.assign(universe_mode="fixed").to_csv(
        results_root / "primary_fixed_loo.tsv", sep="\t", index=False
    )

    rng = np.random.default_rng(random_state)
    non_gejman_samples = metadata.loc[metadata["study"] != "Gejman", "sample_id"].to_numpy()
    tasks = [
        (
            replicate,
            tuple(sorted(rng.choice(non_gejman_samples, size=15, replace=False).tolist())),
        )
        for replicate in range(1, n_replicates + 1)
    ]
    comparator_state: dict[str, object] = {
        "matrix": matrix,
        "full_pc1": full_pc1,
        "saluki_human_pc1": saluki_human_pc1,
        "mouse_pc1": mouse_pc1,
        "mapping": mapping,
        "min_observed_per_gene": 3,
        "n_components": 5,
    }
    context = mp.get_context("fork")
    with context.Pool(
        processes=max(1, n_jobs),
        initializer=_initialize_comparator_worker,
        initargs=(comparator_state,),
    ) as pool:
        comparator_rows = list(pool.imap_unordered(_random_removal_worker, tasks, chunksize=1))
    comparators = pd.DataFrame(comparator_rows).sort_values("replicate")
    comparators.to_csv(
        results_root / "conditional_same_size_deletions.tsv",
        sep="\t",
        index=False,
    )
    observed = primary_dynamic.loc[primary_dynamic["study"] == "Gejman"].iloc[0].to_dict()
    comparator_summary = _conditional_deletion_summary(observed, comparators)
    comparator_summary.to_csv(
        results_root / "conditional_same_size_summary.tsv",
        sep="\t",
        index=False,
    )

    sensitivity_rows: list[dict[str, object]] = []
    parameter_sets = [
        (3, 5, "iterative_pca"),
        (5, 5, "iterative_pca"),
        (10, 5, "iterative_pca"),
        (3, 3, "iterative_pca"),
        (3, 10, "iterative_pca"),
        (3, 5, "sample_median"),
    ]
    for min_observed, rank, imputation_method in parameter_sets:
        for universe_mode in ("dynamic", "fixed"):
            table, label = _leave_one_study_out(
                matrix,
                metadata,
                saluki_human_pc1,
                mouse_pc1,
                mapping,
                min_observed_per_gene=min_observed,
                n_components=rank,
                universe_mode=universe_mode,
                imputation_method=imputation_method,
            )
            ranked = table.sort_values("pc1_stability_pearson").reset_index(drop=True)
            gejman = ranked.loc[ranked["study"] == "Gejman"].iloc[0]
            sensitivity_rows.append(
                {
                    "min_observed_per_gene": min_observed,
                    "pca_imputation_rank": rank,
                    "imputation_method": imputation_method,
                    "universe_mode": universe_mode,
                    "n_full_genes": len(label),
                    "worst_left_out": ranked.iloc[0]["study"],
                    "gejman_rank": int(ranked.index[ranked["study"] == "Gejman"][0] + 1),
                    "gejman_stability_pearson": float(gejman["pc1_stability_pearson"]),
                    "gejman_delta_saluki_pearson": float(gejman["delta_saluki_pearson"]),
                    "gejman_delta_ortholog_pearson": float(gejman["delta_ortholog_pearson"]),
                }
            )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(results_root / "pipeline_sensitivity_summary.tsv", sep="\t", index=False)

    weighted_loo, weighted_summary = _study_weighted_analysis(
        matrix,
        metadata,
        saluki_human_pc1,
        min_observed_per_gene=3,
        n_components=5,
    )
    weighted_loo.to_csv(results_root / "study_weighted_loo.tsv", sep="\t", index=False)
    weighted_summary.to_csv(results_root / "study_weighted_summary.tsv", sep="\t", index=False)

    manifest = {
        "n_replicates": n_replicates,
        "random_state": random_state,
        "n_jobs": n_jobs,
        "conditional_comparator_sampling_pool": "39 non-Gejman human samples; Gejman retained",
        "samples_removed_per_replicate": 15,
        "primary_min_observed_per_gene": 3,
        "primary_pca_imputation_rank": 5,
        "parameter_sensitivity": [
            {
                "min_observed_per_gene": threshold,
                "pca_imputation_rank": rank,
                "imputation_method": imputation_method,
            }
            for threshold, rank, imputation_method in parameter_sets
        ],
    }
    with (results_root / "analysis_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    return {
        "primary_dynamic": results_root / "primary_dynamic_loo.tsv",
        "primary_fixed": results_root / "primary_fixed_loo.tsv",
        "conditional_deletions": results_root / "conditional_same_size_deletions.tsv",
        "conditional_summary": results_root / "conditional_same_size_summary.tsv",
        "pipeline_sensitivity": results_root / "pipeline_sensitivity_summary.tsv",
        "study_weighted_loo": results_root / "study_weighted_loo.tsv",
        "study_weighted_summary": results_root / "study_weighted_summary.tsv",
        "manifest": results_root / "analysis_manifest.json",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Study-influence sample-size and gene-universe sensitivity analysis.")
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "study_influence_sensitivity")
    parser.add_argument("--n-replicates", type=int, default=DEFAULT_N_REPLICATES)
    parser.add_argument("--random-state", type=int, default=DEFAULT_RANDOM_STATE)
    parser.add_argument("--n-jobs", type=int, default=8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_study_influence_sensitivity(
        results_root=args.results_root,
        n_replicates=args.n_replicates,
        random_state=args.random_state,
        n_jobs=args.n_jobs,
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
