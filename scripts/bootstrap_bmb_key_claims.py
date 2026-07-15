#!/usr/bin/env python3
"""Paired gene-bootstrap uncertainty for the central BMB manuscript claims."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
PROCESSED = ROOT / "data" / "processed" / "real_data"


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def _paired_bootstrap(
    reference: np.ndarray,
    candidate_a: np.ndarray,
    candidate_b: np.ndarray,
    *,
    n_bootstrap: int,
    random_state: int,
) -> dict[str, float | int]:
    rng = np.random.default_rng(random_state)
    n = reference.shape[0]
    observed_a = _corr(reference, candidate_a)
    observed_b = _corr(reference, candidate_b)
    correlations_a = np.empty(n_bootstrap, dtype=float)
    correlations_b = np.empty(n_bootstrap, dtype=float)
    differences = np.empty(n_bootstrap, dtype=float)
    for idx in range(n_bootstrap):
        sample = rng.integers(0, n, size=n)
        correlations_a[idx] = _corr(reference[sample], candidate_a[sample])
        correlations_b[idx] = _corr(reference[sample], candidate_b[sample])
        differences[idx] = correlations_a[idx] - correlations_b[idx]
    a_low, a_high = np.quantile(correlations_a, [0.025, 0.975])
    b_low, b_high = np.quantile(correlations_b, [0.025, 0.975])
    d_low, d_high = np.quantile(differences, [0.025, 0.975])
    return {
        "n": int(n),
        "n_bootstrap": int(n_bootstrap),
        "pearson_a": observed_a,
        "pearson_a_ci95_low": float(a_low),
        "pearson_a_ci95_high": float(a_high),
        "pearson_b": observed_b,
        "pearson_b_ci95_low": float(b_low),
        "pearson_b_ci95_high": float(b_high),
        "difference_a_minus_b": observed_a - observed_b,
        "difference_ci95_low": float(d_low),
        "difference_ci95_high": float(d_high),
        "bootstrap_fraction_le_zero": float(np.mean(differences <= 0.0)),
    }


def _load_saluki_human() -> pd.Series:
    path = PROCESSED / "human" / "moesm3_saluki_human_pc1.tsv"
    frame = pd.read_csv(path, sep="\t")
    return frame.set_index("gene_id")["saluki_human_pc1"].astype(float)


def _load_consensus(path: Path, name: str) -> pd.Series:
    frame = pd.read_csv(path, sep="\t")
    return frame.set_index("gene_id")["pc1"].astype(float).rename(name)


def _three_seed_ensemble(setting: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for seed in (13, 42, 101):
        path = RESULTS / "tenfold_global_prior" / "oof_predictions" / f"{setting}__folds10__seed{seed}.tsv"
        frame = pd.read_csv(path, sep="\t").sort_values("gene_id").reset_index(drop=True)
        frames.append(frame)
    out = frames[0][["gene_id", "saluki_human_pc1"]].copy()
    out[setting] = np.mean([frame["prediction"].to_numpy(dtype=float) for frame in frames], axis=0)
    return out


def run(*, output_dir: Path, n_bootstrap: int) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | int | str]] = []

    full = _load_consensus(RESULTS / "real_human_moesm2" / "consensus_labels.tsv", "human_full_pc1")
    no_gejman = _load_consensus(
        RESULTS / "real_human_moesm2_no_gejman" / "consensus_labels.tsv",
        "human_no_gejman_pc1",
    )
    saluki = _load_saluki_human().rename("saluki_human_pc1")
    agreement = pd.concat([saluki, no_gejman, full], axis=1).dropna()
    rows.append(
        {
            "claim": "Saluki agreement after Gejman exclusion",
            "candidate_a": "human_no_gejman_pc1",
            "candidate_b": "human_full_pc1",
            **_paired_bootstrap(
                agreement["saluki_human_pc1"].to_numpy(dtype=float),
                agreement["human_no_gejman_pc1"].to_numpy(dtype=float),
                agreement["human_full_pc1"].to_numpy(dtype=float),
                n_bootstrap=n_bootstrap,
                random_state=20260713,
            ),
        }
    )

    ortholog = pd.read_csv(RESULTS / "ortholog_analysis" / "common_ortholog_pairs.tsv", sep="\t").dropna(
        subset=["mouse_pc1", "human_pruned_pc1", "human_full_pc1"]
    )
    for subset_name, subset in (
        ("all one-to-one orthologs", ortholog),
        ("Ensembl high-confidence orthologs", ortholog.loc[ortholog["is_high_confidence"] == 1]),
    ):
        rows.append(
            {
                "claim": f"Ortholog concordance after Gejman exclusion: {subset_name}",
                "candidate_a": "human_no_gejman_pc1",
                "candidate_b": "human_full_pc1",
                **_paired_bootstrap(
                    subset["mouse_pc1"].to_numpy(dtype=float),
                    subset["human_pruned_pc1"].to_numpy(dtype=float),
                    subset["human_full_pc1"].to_numpy(dtype=float),
                    n_bootstrap=n_bootstrap,
                    random_state=20260714,
                ),
            }
        )

    settings = {
        "human_only": "compact_all_global",
        "real_prior": "compact_all_plus_both_mouse_priors_global",
        "shuffled_prior": "compact_all_plus_both_mouse_priors_shuffled_global",
    }
    ensembles = {name: _three_seed_ensemble(setting) for name, setting in settings.items()}
    prediction = ensembles["human_only"]
    for name in ("real_prior", "shuffled_prior"):
        prediction = prediction.merge(
            ensembles[name][["gene_id", settings[name]]],
            on="gene_id",
            how="inner",
        )
    for comparison_name, baseline_name in (
        ("Real prior versus human-only", "human_only"),
        ("Real prior versus shuffled prior", "shuffled_prior"),
    ):
        rows.append(
            {
                "claim": comparison_name,
                "candidate_a": "real_prior",
                "candidate_b": baseline_name,
                **_paired_bootstrap(
                    prediction["saluki_human_pc1"].to_numpy(dtype=float),
                    prediction[settings["real_prior"]].to_numpy(dtype=float),
                    prediction[settings[baseline_name]].to_numpy(dtype=float),
                    n_bootstrap=n_bootstrap,
                    random_state=20260715,
                ),
            }
        )

    output_path = output_dir / "paired_bootstrap_summary.tsv"
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    metadata = {
        "bootstrap_unit": "gene or one-to-one ortholog pair",
        "bootstrap_iterations": n_bootstrap,
        "interval": "percentile 95% confidence interval",
        "prior_prediction_summary": "mean of three aligned 10-fold OOF prediction vectors (seeds 13, 42, 101)",
        "difference_definition": "Pearson(candidate_a, reference) - Pearson(candidate_b, reference)",
    }
    (output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bootstrap key BMB manuscript correlation differences.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS / "key_claim_bootstrap")
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(run(output_dir=args.output_dir, n_bootstrap=args.bootstrap_iterations))


if __name__ == "__main__":
    main()
