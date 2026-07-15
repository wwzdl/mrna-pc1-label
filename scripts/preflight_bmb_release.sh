#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source scripts/activate_env.sh

echo "[preflight] py_compile"
python3 -m py_compile \
  scripts/audit_bmb_submission_package.py \
  scripts/check_bmb_cn_consistency.py \
  scripts/render_markdown_to_docx.py \
  scripts/render_bmb_markdown_with_references.py \
  scripts/render_markdown_to_pdf.py

bash -n scripts/assemble_bmb_upload_candidate.sh

echo "[preflight] manuscript consistency"
python3 scripts/check_bmb_cn_consistency.py --lang both

echo "[preflight] refresh BMB DOCX bundle"
bash scripts/render_bmb_submission_bundle.sh

echo "[preflight] BMB submission package audit"
python3 scripts/audit_bmb_submission_package.py

if [[ "${BMB_SKIP_PDF:-0}" == "1" ]]; then
  echo "[skip] upload-candidate assembly requires the frozen supplementary PDF"
else
  echo "[preflight] assemble BMB upload candidate"
  bash scripts/assemble_bmb_upload_candidate.sh
fi

echo "[done] BMB release preflight passed."
