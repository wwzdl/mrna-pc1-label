#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source scripts/activate_env.sh

echo "[preflight] Python syntax"
python -m compileall -q src/mrna_half_life_paper
mapfile -t PYTHON_SCRIPTS < <(
  find scripts -maxdepth 1 -type f -name '*.py' -print | sort
)
python -m py_compile "${PYTHON_SCRIPTS[@]}"

echo "[preflight] shell syntax"
mapfile -t SHELL_SCRIPTS < <(
  find scripts -maxdepth 1 -type f -name '*.sh' -print | sort
)
bash -n "${SHELL_SCRIPTS[@]}"

echo "[preflight] focused unit tests"
python -m pytest -q \
  tests/test_prior_ablation_inputs.py \
  tests/test_study_influence_sensitivity.py

if [[ -f results/ortholog_label_distance/ortholog_label_values.tsv ]]; then
  echo "[preflight] analysis-universe ledger"
  python scripts/audit_analysis_universes.py --check
else
  echo "[skip] full analysis-universe audit requires regenerated per-gene label values"
fi

if [[ -d results/tenfold_global_prior/oof_predictions ]]; then
  echo "[preflight] OOF integrity"
  python scripts/audit_oof_integrity.py
else
  echo "[skip] OOF integrity audit requires regenerated fold-level predictions"
fi

echo "[done] public reproducibility preflight passed."
