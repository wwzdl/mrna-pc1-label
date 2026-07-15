#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${ROOT_DIR}/scripts/activate_env.sh"

if [ -z "${VIRTUAL_ENV:-}" ]; then
  echo "[install_optional_extras] no active virtual environment after activation" >&2
  exit 1
fi

"${VIRTUAL_ENV}/bin/pip" install \
  "pymc>=5" \
  arviz \
  xarray \
  r_pca

echo "[install_optional_extras] installed Bayesian extras"
