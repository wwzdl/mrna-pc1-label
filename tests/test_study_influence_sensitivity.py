from __future__ import annotations

import pandas as pd

from mrna_half_life_paper.study_influence_sensitivity import (
    _conditional_deletion_summary,
    _fixed_gene_universe,
)


def test_conditional_deletion_summary_uses_prespecified_tails() -> None:
    observed = {
        "pc1_stability_pearson": 0.05,
        "delta_saluki_pearson": 1.0,
        "delta_ortholog_pearson": 1.0,
    }
    comparators = pd.DataFrame(
        {
            "pc1_stability_pearson": [0.10, 0.20, 0.30],
            "delta_saluki_pearson": [0.10, 0.20, 0.30],
            "delta_ortholog_pearson": [0.10, 0.20, 0.30],
        }
    )

    summary = _conditional_deletion_summary(observed, comparators).set_index("metric")

    assert summary.loc["pc1_stability_pearson", "extreme_direction"] == "lower"
    assert summary.loc["delta_saluki_pearson", "extreme_direction"] == "higher"
    assert summary.loc["delta_ortholog_pearson", "extreme_direction"] == "higher"
    assert (summary["empirical_tail_proportion"] == 0.25).all()


def test_fixed_gene_universe_requires_coverage_after_every_study_removal() -> None:
    matrix = pd.DataFrame(
        {
            "a1": [1.0, 1.0, 1.0],
            "a2": [2.0, 2.0, None],
            "b1": [3.0, None, 3.0],
            "b2": [4.0, None, None],
        },
        index=["always_covered", "study_a_only", "one_per_study"],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["a1", "a2", "b1", "b2"],
            "study": ["A", "A", "B", "B"],
        }
    )

    fixed = _fixed_gene_universe(matrix, metadata, min_observed_per_gene=2)

    assert fixed.tolist() == ["always_covered"]
