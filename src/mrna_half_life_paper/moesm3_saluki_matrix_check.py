from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import METADATA_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.pc1_consensus import compute_pc1_consensus
from mrna_half_life_paper.preprocess import PreprocessConfig, preprocess_matrix, zscore_by_sample
from mrna_half_life_paper.saluki_naming import saluki_label_column, saluki_label_filename


def _load_species_inputs(species: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    matrix_path = PROCESSED_DIR / "real_data" / species / "moesm3_saluki_processed_matrix.tsv"
    saluki_path = PROCESSED_DIR / "real_data" / species / saluki_label_filename(species)
    metadata_path = METADATA_DIR / f"real_{species}_samples.tsv"

    matrix = load_matrix(matrix_path)
    metadata = load_metadata(metadata_path)
    matrix, metadata = align_matrix_and_metadata(matrix, metadata)

    saluki_frame = pd.read_csv(saluki_path, sep="\t")
    saluki_col = saluki_label_column(species)
    saluki_label = saluki_frame.set_index("gene_id")[saluki_col].astype(float).rename(saluki_col)
    return matrix, metadata, saluki_label


def _subset_matrix_and_meta(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    subset_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if subset_name == "all_samples":
        return matrix.copy(), metadata.copy()
    if subset_name == "no_gejman":
        keep_meta = metadata.loc[metadata["study"] != "Gejman"].copy()
        keep_samples = keep_meta["sample_id"].tolist()
        return matrix.loc[:, keep_samples].copy(), keep_meta
    raise ValueError(f"Unknown subset: {subset_name}")


def _compute_variant_labels(matrix: pd.DataFrame) -> dict[str, pd.Series]:
    variants: dict[str, pd.Series] = {}
    variants["as_is_pc1"] = compute_pc1_consensus(matrix)
    variants["zscore_only_pc1"] = compute_pc1_consensus(zscore_by_sample(matrix))
    rerun = preprocess_matrix(matrix, PreprocessConfig(apply_log=False))
    variants["rerun_pipeline_pc1"] = compute_pc1_consensus(rerun)
    return variants


def run_moesm3_saluki_matrix_check(results_root: Path | None = None) -> pd.DataFrame:
    results_root = (
        results_root if results_root is not None else RESULTS_DIR / "moesm3_saluki_matrix_check"
    )
    results_root.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []

    for species in ("human", "mouse"):
        species_dir = results_root / species
        species_dir.mkdir(parents=True, exist_ok=True)
        matrix, metadata, saluki_label = _load_species_inputs(species)

        subsets = ["all_samples"]
        if species == "human":
            subsets.append("no_gejman")

        for subset_name in subsets:
            subset_matrix, subset_meta = _subset_matrix_and_meta(matrix, metadata, subset_name)
            variant_labels = _compute_variant_labels(subset_matrix)
            subset_dir = species_dir / subset_name
            subset_dir.mkdir(parents=True, exist_ok=True)
            subset_meta.to_csv(subset_dir / "metadata.tsv", sep="\t", index=False)

            for variant_name, labels in variant_labels.items():
                labels.rename("pc1").to_frame().reset_index().rename(
                    columns={"index": "gene_id"}
                ).to_csv(subset_dir / f"{variant_name}.tsv", sep="\t", index=False)
                corr = label_correlation(labels, saluki_label)
                summary_rows.append(
                    {
                        "species": species,
                        "subset": subset_name,
                        "variant": variant_name,
                        "n_samples": int(subset_matrix.shape[1]),
                        "n_genes_labels": int(labels.shape[0]),
                        "n_common_genes": int(corr["n"]),
                        "pearson_vs_saluki_label": corr["pearson"],
                        "spearman_vs_saluki_label": corr["spearman"],
                    }
                )

    summary = pd.DataFrame(summary_rows).sort_values(
        ["species", "subset", "pearson_vs_saluki_label"], ascending=[True, True, False]
    )
    summary.to_csv(results_root / "summary.tsv", sep="\t", index=False)
    with open(results_root / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary_rows, handle, indent=2)
    return summary


def main() -> None:
    summary = run_moesm3_saluki_matrix_check()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
