from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
METADATA_DIR = ROOT / "metadata"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
MANUSCRIPT_DIR = ROOT / "manuscript"


def ensure_project_dirs() -> None:
    for path in (
        DATA_DIR,
        RAW_DIR,
        INTERIM_DIR,
        PROCESSED_DIR,
        METADATA_DIR,
        RESULTS_DIR,
        FIGURES_DIR,
        MANUSCRIPT_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)
