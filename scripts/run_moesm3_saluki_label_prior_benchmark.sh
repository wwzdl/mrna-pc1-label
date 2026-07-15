#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

source scripts/activate_env.sh
python -m mrna_half_life_paper.moesm3_saluki_label_prior_benchmark "$@"
