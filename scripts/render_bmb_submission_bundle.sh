#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

source scripts/activate_env.sh

render_one() {
  local input_md="$1"
  local output_docx="$2"
  echo "[render] ${input_md} -> ${output_docx}"
  python3 scripts/render_markdown_to_docx.py "${input_md}" "${output_docx}"
}

render_one \
  manuscript/bmb_submission/bmb_main_manuscript_cn.md \
  manuscript/bmb_submission/bmb_main_manuscript_cn.docx

render_one \
  manuscript/bmb_submission/bmb_supplementary_material_cn.md \
  manuscript/bmb_submission/bmb_supplementary_material_cn.docx

render_one \
  manuscript/bmb_submission/bmb_main_manuscript_en.md \
  manuscript/bmb_submission/bmb_main_manuscript_en.docx

render_one \
  manuscript/bmb_submission/bmb_supplementary_material_en.md \
  manuscript/bmb_submission/bmb_supplementary_material_en.docx

render_one \
  manuscript/bmb_submission/bmb_title_page_cn.md \
  manuscript/bmb_submission/bmb_title_page_cn.docx

render_one \
  manuscript/bmb_submission/bmb_title_page_en.md \
  manuscript/bmb_submission/bmb_title_page_en.docx

render_one \
  manuscript/bmb_submission/bmb_cover_letter_cn.md \
  manuscript/bmb_submission/bmb_cover_letter_cn.docx

render_one \
  manuscript/bmb_submission/bmb_cover_letter_en.md \
  manuscript/bmb_submission/bmb_cover_letter_en.docx

render_one \
  manuscript/bmb_submission/bmb_statements_and_declarations_cn.md \
  manuscript/bmb_submission/bmb_statements_and_declarations_cn.docx

render_one \
  manuscript/bmb_submission/bmb_statements_and_declarations_en.md \
  manuscript/bmb_submission/bmb_statements_and_declarations_en.docx

render_one \
  manuscript/bmb_submission/bmb_author_contributions_cn.md \
  manuscript/bmb_submission/bmb_author_contributions_cn.docx

render_one \
  manuscript/bmb_submission/bmb_editor_highlights_cn.md \
  manuscript/bmb_submission/bmb_editor_highlights_cn.docx

if [[ "${BMB_SKIP_PDF:-0}" == "1" ]]; then
  echo "[skip] supplementary PDF export (BMB_SKIP_PDF=1)"
else
  echo "[render] supplementary Markdown -> paginated PDF"
  python3 scripts/render_markdown_to_pdf.py \
    manuscript/bmb_submission/bmb_supplementary_material_cn.md \
    manuscript/bmb_submission/bmb_supplementary_material_cn.pdf
  python3 scripts/render_markdown_to_pdf.py \
    manuscript/bmb_submission/bmb_supplementary_material_en.md \
    manuscript/bmb_submission/bmb_supplementary_material_en.pdf
fi

echo "[done] BMB manuscript bundle refreshed."
