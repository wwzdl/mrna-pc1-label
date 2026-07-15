#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT_DIR/../.venv/bin/activate"
PYTHONPATH="$ROOT_DIR/src" python -m mrna_half_life_paper.moesm3_saluki_matrix_check
