#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_PREFIX="${MRNA_HALF_LIFE_ENV_PREFIX:-${ROOT_DIR}/.venv}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

echo "[create_env] root: ${ROOT_DIR}"
echo "[create_env] env prefix: ${ENV_PREFIX}"
echo "[create_env] python: ${PYTHON_BIN}"

mkdir -p "$(dirname "${ENV_PREFIX}")"
"${PYTHON_BIN}" -m venv "${ENV_PREFIX}"
"${ENV_PREFIX}/bin/python" -m pip install --upgrade pip
"${ENV_PREFIX}/bin/python" -m pip install -r "${ROOT_DIR}/requirements.txt"

echo "[create_env] done"
echo "[create_env] activate with: source scripts/activate_env.sh"
