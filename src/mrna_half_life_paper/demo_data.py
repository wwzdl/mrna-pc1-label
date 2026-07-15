from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mrna_half_life_paper.config import METADATA_DIR, RAW_DIR


def generate_demo_data(
    out_matrix: str | Path | None = None,
    out_metadata: str | Path | None = None,
    random_state: int = 42,
) -> tuple[Path, Path]:
    rng = np.random.default_rng(random_state)

    genes = [f"gene_{i:04d}" for i in range(1, 201)]
    samples = [
        "sample_01",
        "sample_02",
        "sample_03",
        "sample_04",
        "sample_05",
        "sample_06",
        "sample_07",
        "sample_08",
    ]
    study = ["study_a", "study_a", "study_b", "study_b", "study_c", "study_c", "study_d", "study_d"]
    method = ["4sU", "4sU", "ActD", "ActD", "4sU", "4sU", "ActD", "ActD"]
    cell_type = ["HEK293", "HEK293", "HeLa", "HeLa", "K562", "K562", "HepG2", "HepG2"]

    theta = rng.normal(loc=0.0, scale=1.0, size=len(genes))
    study_effects = {"study_a": 0.15, "study_b": -0.25, "study_c": 0.20, "study_d": -0.10}
    method_effects = {"4sU": 0.12, "ActD": -0.12}

    matrix = pd.DataFrame(index=genes, columns=samples, dtype=float)
    for sample, s, m in zip(samples, study, method):
        noise = rng.normal(loc=0.0, scale=0.25, size=len(genes))
        values = theta + study_effects[s] + method_effects[m] + noise
        matrix[sample] = values

    missing_mask = rng.random(matrix.shape) < 0.08
    matrix = matrix.mask(missing_mask)
    matrix = np.exp(matrix) + 0.2

    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "study": study,
            "method": method,
            "cell_type": cell_type,
            "species": ["human"] * len(samples),
        }
    )

    out_matrix = Path(out_matrix) if out_matrix is not None else RAW_DIR / "demo_half_life_matrix.tsv"
    out_metadata = Path(out_metadata) if out_metadata is not None else METADATA_DIR / "demo_sample_metadata.tsv"
    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    out_metadata.parent.mkdir(parents=True, exist_ok=True)

    matrix.reset_index().rename(columns={"index": "gene_id"}).to_csv(out_matrix, sep="\t", index=False)
    metadata.to_csv(out_metadata, sep="\t", index=False)
    return out_matrix, out_metadata


def main() -> None:
    out_matrix, out_metadata = generate_demo_data()
    print(out_matrix)
    print(out_metadata)


if __name__ == "__main__":
    main()
