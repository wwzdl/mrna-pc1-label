#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${ROOT_DIR}/scripts/activate_env.sh"
python "${ROOT_DIR}/scripts/reproduce_sample_pca_imputation.py" "$@"
