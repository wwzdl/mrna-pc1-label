#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: bash scripts/run_paper1_workflow.sh <matrix.tsv> <metadata.tsv> <outdir>" >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MATRIX_PATH="$1"
METADATA_PATH="$2"
OUTDIR="$3"

source "${ROOT_DIR}/scripts/activate_env.sh"

python -m mrna_half_life_paper.run_workflow \
  --matrix "${MATRIX_PATH}" \
  --metadata "${METADATA_PATH}" \
  --outdir "${OUTDIR}"

echo "[run_paper1_workflow] outputs written to ${OUTDIR}"
