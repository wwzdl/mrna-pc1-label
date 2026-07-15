from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import METADATA_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.real_data import extract_real_data
from mrna_half_life_paper.run_workflow import run_workflow
from mrna_half_life_paper.saluki_naming import saluki_label_column, saluki_label_filename


def _read_saluki_label(path: Path, species: str) -> pd.Series:
    df = pd.read_csv(path, sep="\t")
    column = saluki_label_column(species)
    return df.set_index("gene_id")[column].astype(float)


def _run_one_species(species: str, extracted_root: Path, metadata_root: Path, results_root: Path) -> None:
    matrix_path = extracted_root / species / "moesm2_raw_sparse_matrix.tsv"
    metadata_path = metadata_root / f"real_{species}_samples.tsv"
    saluki_label_path = extracted_root / species / saluki_label_filename(species)
    outdir = results_root / f"real_{species}_moesm2"

    run_workflow(
        matrix_path=matrix_path,
        metadata_path=metadata_path,
        outdir=outdir,
        apply_log=False,
    )

    labels = pd.read_csv(outdir / "consensus_labels.tsv", sep="\t").set_index("gene_id")
    saluki_label = _read_saluki_label(saluki_label_path, species)

    rows = []
    for label_name in ("mean", "pc1", "bias_aware"):
        corr = label_correlation(labels[label_name], saluki_label)
        rows.append(
            {
                "species": species,
                "label": label_name,
                "n_genes_overlap": int(corr["n"]),
                "pearson_vs_saluki_label": corr["pearson"],
                "spearman_vs_saluki_label": corr["spearman"],
            }
        )

    pd.DataFrame(rows).to_csv(outdir / "saluki_label_comparison.tsv", sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the paper workflow on the public Saluki half-life supplements.")
    parser.add_argument(
        "--species",
        nargs="+",
        default=["human", "mouse"],
        choices=["human", "mouse"],
        help="Species to analyze.",
    )
    parser.add_argument(
        "--extracted-root",
        type=Path,
        default=PROCESSED_DIR / "real_data",
        help="Directory containing extracted supplement TSVs.",
    )
    parser.add_argument(
        "--metadata-root",
        type=Path,
        default=METADATA_DIR,
        help="Directory containing sample metadata TSVs.",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR,
        help="Directory for workflow outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extract_real_data(args.extracted_root, args.metadata_root)
    for species in args.species:
        _run_one_species(species, args.extracted_root, args.metadata_root, args.results_root)
        print(args.results_root / f"real_{species}_moesm2")


if __name__ == "__main__":
    main()
