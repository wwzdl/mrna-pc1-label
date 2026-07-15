from __future__ import annotations

import importlib
import platform
import sys

from mrna_half_life_paper.config import ROOT, ensure_project_dirs


def main() -> None:
    ensure_project_dirs()

    modules = [
        "numpy",
        "pandas",
        "scipy",
        "sklearn",
        "statsmodels",
        "matplotlib",
        "seaborn",
        "xgboost",
        "openpyxl",
        "docx",
        "PIL",
        "markdown",
        "weasyprint",
        "bibtexparser",
    ]

    print("Smoke test for mRNA half-life paper environment")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    print(f"Project root: {ROOT}")

    for module in modules:
        importlib.import_module(module)
        print(f"[OK] import {module}")

    print("[OK] project directories ready")


if __name__ == "__main__":
    main()
