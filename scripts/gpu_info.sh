#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "${ROOT_DIR}/scripts/activate_env.sh" >/dev/null

echo "[gpu_info] GPUs"
nvidia-smi -L
echo
echo "[gpu_info] Driver and memory"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
echo
echo "[gpu_info] nvcc"
"${CUDA_HOME}/bin/nvcc" --version
