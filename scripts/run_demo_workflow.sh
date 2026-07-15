#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${ROOT_DIR}/scripts/activate_env.sh"
python -m mrna_half_life_paper.demo_data
python -m mrna_half_life_paper.run_workflow \
  --matrix "${ROOT_DIR}/data/raw/demo_half_life_matrix.tsv" \
  --metadata "${ROOT_DIR}/metadata/demo_sample_metadata.tsv" \
  --outdir "${ROOT_DIR}/results/demo_workflow"

echo "[run_demo_workflow] outputs written to results/demo_workflow"
