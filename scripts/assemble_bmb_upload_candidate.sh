#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUB_DIR="${ROOT_DIR}/manuscript/bmb_submission"
OUT_DIR="${SUB_DIR}/upload_candidate"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}/figures/main" "${OUT_DIR}/figures/supplement" "${OUT_DIR}/references"

copy_one() {
  local src="$1"
  local dst="$2"
  if [[ ! -f "${src}" ]]; then
    echo "[missing] ${src}" >&2
    exit 1
  fi
  cp "${src}" "${dst}"
}

copy_one "${SUB_DIR}/bmb_main_manuscript_en.docx" "${OUT_DIR}/main_manuscript.docx"
copy_one "${SUB_DIR}/bmb_title_page_en.docx" "${OUT_DIR}/title_page.docx"
copy_one "${SUB_DIR}/bmb_cover_letter_en.docx" "${OUT_DIR}/cover_letter.docx"
copy_one "${SUB_DIR}/bmb_statements_and_declarations_en.docx" "${OUT_DIR}/statements_and_declarations.docx"
copy_one "${SUB_DIR}/bmb_supplementary_material_en.pdf" "${OUT_DIR}/Online_Resource_1_supplementary_material.pdf"

copy_one "${SUB_DIR}/figures/main/Fig01_workflow.tiff" "${OUT_DIR}/figures/main/Figure_1.tiff"
copy_one "${SUB_DIR}/figures/main/Fig02_dataset_summary_cn.tiff" "${OUT_DIR}/figures/main/Figure_2.tiff"
copy_one "${SUB_DIR}/figures/main/Fig03_human_pca_panel.tiff" "${OUT_DIR}/figures/main/Figure_3.tiff"
copy_one "${SUB_DIR}/figures/main/Fig04_human_study_influence.tiff" "${OUT_DIR}/figures/main/Figure_4.tiff"
copy_one "${SUB_DIR}/figures/main/Fig05_ortholog_summary_cn.tiff" "${OUT_DIR}/figures/main/Figure_5.tiff"
copy_one "${SUB_DIR}/figures/main/Fig06_sequence_progression_cn.tiff" "${OUT_DIR}/figures/main/Figure_6.tiff"
copy_one "${SUB_DIR}/figures/main/Fig07_prior_summary_cn.tiff" "${OUT_DIR}/figures/main/Figure_7.tiff"
copy_one "${SUB_DIR}/figures/main/Fig08_ortholog_regularized_target_cn.tiff" "${OUT_DIR}/figures/main/Figure_8.tiff"

copy_one "${SUB_DIR}/figures/supplement/FigS01_prediction_boundary_flow.tiff" "${OUT_DIR}/figures/supplement/Online_Resource_1_Figure_S1.tiff"
copy_one "${SUB_DIR}/figures/supplement/FigS_mouse_pca_panel.tiff" "${OUT_DIR}/figures/supplement/Online_Resource_1_Figure_S2.tiff"
copy_one "${SUB_DIR}/figures/supplement/FigS_mouse_study_influence.tiff" "${OUT_DIR}/figures/supplement/Online_Resource_1_Figure_S3.tiff"
copy_one "${SUB_DIR}/figures/supplement/FigS_study_influence_sensitivity.tiff" "${OUT_DIR}/figures/supplement/Online_Resource_1_Figure_S4.tiff"

copy_one "${SUB_DIR}/references/bmb_references.bib" "${OUT_DIR}/references/bmb_references.bib"
copy_one "${SUB_DIR}/references/bmb_references_for_endnote.ris" "${OUT_DIR}/references/bmb_references_for_endnote.ris"

copy_one "${SUB_DIR}/upload_candidate_manifest.md" "${OUT_DIR}/MANIFEST.md"

echo "[done] upload candidate: ${OUT_DIR}"
