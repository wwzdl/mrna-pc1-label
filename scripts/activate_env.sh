#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PARENT_DIR="$(cd "${ROOT_DIR}/.." && pwd)"

resolve_env_prefix() {
  if [ -d "${ROOT_DIR}/.venv" ]; then
    printf '%s\n' "${ROOT_DIR}/.venv"
    return 0
  fi
  if [ -d "${PARENT_DIR}/.venv" ]; then
    printf '%s\n' "${PARENT_DIR}/.venv"
    return 0
  fi
  printf '%s\n' "${ROOT_DIR}/.venv"
}

ENV_PREFIX="$(resolve_env_prefix)"

if [ ! -f "${ENV_PREFIX}/bin/activate" ]; then
  echo "[activate_env] missing virtual environment at ${ENV_PREFIX}" >&2
  echo "[activate_env] run: bash scripts/create_env.sh" >&2
  return 1 2>/dev/null || exit 1
fi
source "${ENV_PREFIX}/bin/activate"

export PYTHONPATH="${ROOT_DIR}/src:${PYTHONPATH:-}"
export MRNA_HALF_LIFE_PROJECT_ROOT="${ROOT_DIR}"
if [ -d "/usr/local/cuda-12.8" ]; then
  export CUDA_HOME="/usr/local/cuda-12.8"
  export PATH="${CUDA_HOME}/bin:${PATH}"
  export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH:-}"
fi

echo "[activate_env] active env: ${VIRTUAL_ENV}"
echo "[activate_env] project root: ${ROOT_DIR}"
echo "[activate_env] PYTHONPATH=${PYTHONPATH}"
if [ -n "${CUDA_HOME:-}" ]; then
  echo "[activate_env] CUDA_HOME=${CUDA_HOME}"
fi
