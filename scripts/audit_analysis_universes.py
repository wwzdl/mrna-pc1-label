#!/usr/bin/env python3
"""Build and validate the analysis-universe ledger used by the BMB manuscript."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DATA = ROOT / "data" / "processed" / "real_data"
SUBMISSION = ROOT / "manuscript" / "bmb_submission"
LEDGER_PATH = RESULTS / "bmb_analysis_universe_ledger.tsv"


def _ids(path: Path, column: str = "gene_id") -> set[str]:
    frame = pd.read_csv(path, sep="\t", usecols=[column])
    values = frame[column].astype(str)
    if values.duplicated().any():
        raise AssertionError(f"Duplicate {column} values in {path}")
    return set(values)


def _pair_ids(path: Path) -> tuple[pd.DataFrame, set[str]]:
    frame = pd.read_csv(path, sep="\t")
    pair_cols = ["human_gene_id", "mouse_gene_id"]
    if frame.duplicated(pair_cols).any():
        raise AssertionError(f"Duplicate ortholog pairs in {path}")
    return frame, set(frame["human_gene_id"].astype(str))


def _expect(name: str, observed: int, expected: int) -> None:
    if observed != expected:
        raise AssertionError(f"{name}: expected {expected:,}, observed {observed:,}")


def build_ledger() -> pd.DataFrame:
    human_full = _ids(RESULTS / "real_human_moesm2" / "consensus_labels.tsv")
    mouse_full = _ids(RESULTS / "real_mouse_moesm2" / "consensus_labels.tsv")
    human_no_gejman = _ids(RESULTS / "real_human_moesm2_no_gejman" / "consensus_labels.tsv")
    saluki_human = _ids(DATA / "human" / "moesm3_saluki_human_pc1.tsv")
    features = _ids(RESULTS / "human_sequence_features_regions3mer4mer.tsv.gz")

    global_prediction = _ids(
        RESULTS
        / "tenfold_global_prior"
        / "oof_predictions"
        / "compact_all_global__folds10__seed13.tsv"
    )
    both_priors = _ids(RESULTS / "prior_residual_analysis" / "oof_predictions.tsv")
    shrinkage_construction = _ids(
        RESULTS
        / "ortholog_regularized_label_10fold_main"
        / "targets"
        / "orthoreg_reconstructed_mouse_pc1_0.1.tsv"
    )
    shrinkage_prediction = _ids(
        RESULTS
        / "ortholog_regularized_cross_target"
        / "oof_predictions"
        / "cross_target__folds10__seed13.tsv.gz"
    )

    ortholog_audit_frame, ortholog_audit = _pair_ids(
        RESULTS / "ortholog_analysis" / "common_ortholog_pairs.tsv"
    )
    ortholog_high_frame, ortholog_high = _pair_ids(
        RESULTS / "ortholog_analysis" / "common_ortholog_pairs_high_confidence.tsv"
    )
    shrinkage_pair_frame, shrinkage_pairs = _pair_ids(
        RESULTS / "ortholog_label_distance" / "ortholog_label_values.tsv"
    )
    shrinkage_high_frame = shrinkage_pair_frame.loc[
        shrinkage_pair_frame["is_high_confidence"].astype(int) == 1
    ]
    shrinkage_high_pairs = set(shrinkage_high_frame["human_gene_id"].astype(str))

    fixed_loo = pd.read_csv(
        RESULTS / "study_influence_sensitivity" / "primary_fixed_loo.tsv",
        sep="\t",
    )
    fixed_loo_counts = set(fixed_loo["stability_n"].astype(int))
    if len(fixed_loo_counts) != 1:
        raise AssertionError(f"Fixed LOO universe is not constant: {sorted(fixed_loo_counts)}")
    fixed_loo_n = fixed_loo_counts.pop()

    feature_target = features & saluki_human
    fixed_human_saluki = human_full & human_no_gejman & saluki_human
    expected_global = features & human_full & human_no_gejman & saluki_human
    expected_shrinkage_prediction = shrinkage_construction & features & saluki_human
    shrinkage_model_pair_overlap = shrinkage_prediction & shrinkage_pairs
    fixed_target_ortholog = global_prediction & ortholog_audit
    fixed_target_high_conf = global_prediction & ortholog_high

    subset_scan = pd.read_csv(RESULTS / "human_compact_all_subset_scan.tsv", sep="\t").set_index("subset")
    reg_top33_n = int(subset_scan.loc["reg_top33", "n"])

    checks = {
        "human full reconstruction": (len(human_full), 16444),
        "mouse full reconstruction": (len(mouse_full), 16951),
        "human no-Gejman reconstruction": (len(human_no_gejman), 15788),
        "fixed LOO sensitivity": (fixed_loo_n, 14244),
        "fixed human-Saluki comparison": (len(fixed_human_saluki), 13265),
        "ortholog concordance": (len(ortholog_audit_frame), 12592),
        "high-confidence ortholog concordance": (len(ortholog_high_frame), 12376),
        "broad feature-target universe": (len(feature_target), 13532),
        "global prediction universe": (len(global_prediction), 12916),
        "both-priors subset": (len(both_priors), 11107),
        "shrinkage target-construction universe": (len(shrinkage_construction), 12644),
        "shrinkage prediction universe": (len(shrinkage_prediction), 12307),
        "shrinkage mapped-pair set": (len(shrinkage_pair_frame), 10768),
        "high-confidence shrinkage pairs": (len(shrinkage_high_frame), 10626),
        "shrinkage prediction/mapped overlap": (len(shrinkage_model_pair_overlap), 10682),
        "fixed-target reconstructed-prior subset": (len(fixed_target_ortholog), 11569),
        "fixed-target high-confidence subset": (len(fixed_target_high_conf), 11389),
        "fixed-target top-regulatory-33 subset": (reg_top33_n, 4262),
    }
    for name, (observed, expected) in checks.items():
        _expect(name, observed, expected)

    if global_prediction != expected_global:
        raise AssertionError("Global prediction IDs do not match features/full/no-Gejman/Saluki intersection")
    if shrinkage_prediction != expected_shrinkage_prediction:
        raise AssertionError("Shrinkage prediction IDs do not match target/features/Saluki intersection")
    if not both_priors <= global_prediction:
        raise AssertionError("Both-priors set is not a subset of the global prediction universe")
    if not shrinkage_pairs <= shrinkage_construction:
        raise AssertionError("Shrinkage mapped-pair set is not a subset of target construction")
    if not ortholog_high <= ortholog_audit:
        raise AssertionError("High-confidence ortholog set is not a subset of all ortholog pairs")
    if not shrinkage_high_pairs <= shrinkage_pairs:
        raise AssertionError("High-confidence shrinkage set is not a subset of mapped shrinkage pairs")

    only_reconstructed = len(fixed_target_ortholog - both_priors)
    missing_both = len(global_prediction - fixed_target_ortholog)
    _expect("only reconstructed mouse prior", only_reconstructed, 462)
    _expect("missing both mouse priors", missing_both, 1347)

    rows = [
        (
            "L1",
            "label reconstruction",
            "core",
            "human full reconstruction",
            16444,
            "genes",
            "MOESM2 human genes observed in at least 3 of 54 samples",
            "Full-human PC1, sample PCA, and the starting point for human LOO auditing",
        ),
        (
            "L2",
            "label reconstruction",
            "core",
            "mouse full reconstruction",
            16951,
            "genes",
            "MOESM2 mouse genes observed in at least 3 of 27 samples",
            "Mouse PC1, mouse PCA/LOO, and the reconstructed mouse prior source",
        ),
        (
            "L3",
            "label reconstruction",
            "core",
            "human no-Gejman reconstruction",
            15788,
            "genes",
            "MOESM2 human genes observed in at least 3 of the 39 non-Gejman samples",
            "No-Gejman PC1 and the primary dynamic-coverage Gejman comparison",
        ),
        (
            "A1",
            "label audit",
            "core sensitivity",
            "fixed LOO sensitivity universe",
            14244,
            "genes",
            "Genes meeting the 3-sample coverage rule in the full data and after every human study removal",
            "Fixed-universe LOO sensitivity analysis only",
        ),
        (
            "A2",
            "label audit",
            "core comparison",
            "fixed human-Saluki comparison",
            13265,
            "genes",
            "Intersection of human full PC1, human no-Gejman PC1, and Saluki human PC1",
            "Fixed-gene Saluki-agreement estimate and paired bootstrap",
        ),
        (
            "A3",
            "label audit",
            "core comparison",
            "human-mouse ortholog concordance",
            12592,
            "ortholog pairs",
            "Ensembl one-to-one pairs with human full/no-Gejman PC1 and reconstructed mouse PC1",
            "Cross-species audit of the Gejman-related label change",
        ),
        (
            "P1",
            "fixed-target prediction",
            "core secondary",
            "broad feature-target universe",
            13532,
            "genes",
            "Intersection of human sequence features and Saluki human PC1",
            "Regulatory-block ablation and MOESM3-derived-label controls",
        ),
        (
            "P2",
            "fixed-target prediction",
            "core primary",
            "global prediction universe",
            12916,
            "genes",
            "P1 genes also present in human full and no-Gejman reconstructions",
            "Primary human-only, cross-species transfer, permutation, and repeated-CV benchmark",
        ),
        (
            "P3",
            "fixed-target prediction",
            "core subset",
            "both-priors subset",
            11107,
            "genes",
            "P2 genes with both reconstructed and Saluki mouse prior values",
            "Prior-complete coverage analysis and two-stage residual decomposition",
        ),
        (
            "S1",
            "target shrinkage",
            "core construction",
            "shrinkage target-construction universe",
            12644,
            "genes",
            "Human no-Gejman genes retained under the Saluki-like 10-sample coverage rule",
            "Construction of the 0.10 ortholog-informed target",
        ),
        (
            "S2",
            "target shrinkage",
            "core primary",
            "shrinkage prediction universe",
            12307,
            "genes",
            "S1 genes with human compact features and Saluki human PC1",
            "Human-only/prior target prediction, ablation, and cross-target evaluation",
        ),
        (
            "S3",
            "target shrinkage",
            "core comparison",
            "shrinkage mapped-pair set",
            10768,
            "ortholog pairs",
            "S1 human genes with an Ensembl one-to-one reconstructed-mouse-PC1 match",
            "Human-mouse label distance, target geometry, and study-noise comparison",
        ),
        (
            "D1",
            "derived subsets",
            "secondary",
            "high-confidence ortholog audit",
            12376,
            "ortholog pairs",
            "A3 pairs with Ensembl is_high_confidence = 1",
            "High-confidence ortholog sensitivity",
        ),
        (
            "D2",
            "derived subsets",
            "secondary",
            "fixed-target reconstructed-prior subset",
            11569,
            "genes",
            "P2 genes present in A3",
            "Prior coverage and ortholog-subset sensitivity",
        ),
        (
            "D3",
            "derived subsets",
            "secondary",
            "fixed-target high-confidence subset",
            11389,
            "genes",
            "P2 genes present in D1",
            "High-confidence gene-subset sensitivity",
        ),
        (
            "D4",
            "derived subsets",
            "secondary",
            "fixed-target top-regulatory-33 subset",
            4262,
            "genes",
            "Top third of P2 by regulatory score",
            "Gene-subset sensitivity in Fig. 6d",
        ),
        (
            "D5",
            "derived subsets",
            "secondary",
            "only reconstructed mouse prior",
            462,
            "genes",
            "P2 genes in D2 but not P3",
            "Prior-availability accounting",
        ),
        (
            "D6",
            "derived subsets",
            "secondary",
            "missing both mouse priors",
            1347,
            "genes",
            "P2 genes absent from D2 and P3",
            "Coverage-aware fallback and deployment boundary",
        ),
        (
            "D7",
            "derived subsets",
            "secondary overlap",
            "shrinkage prediction/mapped overlap",
            10682,
            "genes",
            "Intersection of S2 prediction genes and S3 mapped human genes",
            "Set accounting; S2 and S3 are not nested",
        ),
        (
            "D8",
            "derived subsets",
            "secondary",
            "high-confidence shrinkage pairs",
            10626,
            "ortholog pairs",
            "S3 pairs with Ensembl is_high_confidence = 1",
            "Label-distance sensitivity in Supplementary Table S12",
        ),
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "set_id",
            "track",
            "level",
            "set_name",
            "n",
            "unit",
            "eligibility_rule",
            "main_use",
        ],
    )


def _check_manuscripts() -> None:
    main_required = [
        "16,444",
        "16,951",
        "15,788",
        "14,244",
        "13,265",
        "12,592",
        "13,532",
        "12,916",
        "11,107",
        "12,644",
        "12,307",
        "10,768",
        "10,682",
    ]
    supplement_required = main_required + [
        "12,376",
        "11,569",
        "11,389",
        "4,262",
        "462",
        "1,347",
        "10,626",
    ]
    for language in ("en", "cn"):
        main = (SUBMISSION / f"bmb_main_manuscript_{language}.md").read_text(encoding="utf-8")
        supplement = (SUBMISSION / f"bmb_supplementary_material_{language}.md").read_text(encoding="utf-8")
        missing_main = [value for value in main_required if value not in main]
        missing_supplement = [value for value in supplement_required if value not in supplement]
        if missing_main:
            raise AssertionError(f"{language} main manuscript misses universe counts: {missing_main}")
        if missing_supplement:
            raise AssertionError(f"{language} supplement misses universe counts: {missing_supplement}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Compare the computed ledger with the versioned TSV and verify manuscript counts.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ledger = build_ledger()
    if args.check:
        if not LEDGER_PATH.exists():
            raise FileNotFoundError(f"Missing versioned universe ledger: {LEDGER_PATH}")
        stored = pd.read_csv(LEDGER_PATH, sep="\t")
        pd.testing.assert_frame_equal(stored, ledger, check_dtype=False)
        _check_manuscripts()
        print(f"[ok] analysis-universe ledger: {len(ledger)} fixed/derived rows")
        return 0

    LEDGER_PATH.parent.mkdir(parents=True, exist_ok=True)
    ledger.to_csv(LEDGER_PATH, sep="\t", index=False)
    print(f"Wrote {LEDGER_PATH} ({len(ledger)} rows)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, FileNotFoundError) as exc:
        print(f"[fail] {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
