#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

usage() {
  cat <<'EOF'
Usage: bash scripts/reproduce_bmb_key_results.sh <labels|models|figures|all>

  labels   Download MOESM2/MOESM3 and reproduce PC1, study-influence,
           sample-count sensitivity, ortholog concordance, label distance,
           and study-noise results.
  models   Reproduce fixed-target and target-shrinkage OOF benchmarks.
           This stage requires the Saluki datapack and a CUDA-capable XGBoost setup.
  figures  Rebuild all active BMB figures, DOCX/PDF files, and run preflight checks.
  all      Run labels, models, and figures in that order.

Set MRNA_BMB_DEVICE=cpu only for modules that expose a CPU device option. The
published benchmark used CUDA; global-prior modules currently use CUDA directly.
EOF
}

if [ "$#" -ne 1 ]; then
  usage >&2
  exit 2
fi

STAGE="$1"
DEVICE="${MRNA_BMB_DEVICE:-cuda}"

source scripts/activate_env.sh

require_file() {
  if [ ! -s "$1" ]; then
    echo "[reproduce] missing required file: $1" >&2
    exit 1
  fi
}

run_labels() {
  echo "[reproduce] public matrices and consensus labels"
  bash scripts/prepare_real_data.sh
  bash scripts/run_real_data_workflow.sh
  bash scripts/run_study_influence_analysis.sh --run-pruned-human
  python -m mrna_half_life_paper.study_influence_sensitivity \
    --n-replicates "${MRNA_STUDY_NULL_REPLICATES:-500}" \
    --n-jobs "${MRNA_STUDY_NULL_JOBS:-1}"
  bash scripts/run_ortholog_analysis.sh
  python -m mrna_half_life_paper.ortholog_label_distance
  python -m mrna_half_life_paper.study_noise_vs_orthoreg_shift
  bash scripts/run_sample_pca_imputation_reproducibility.sh
}

run_models() {
  require_file data/raw/saluki_datapack/human/CWCS.txt.gz
  require_file data/raw/saluki_datapack/human/seqFeatWithKmerFreqs.txt.gz
  require_file data/raw/saluki_datapack/human/SeqWeaver_predictions/5pUTR_avg.txt.gz
  require_file data/raw/saluki_datapack/human/SeqWeaver_predictions/ORF_avg.txt.gz
  require_file data/raw/saluki_datapack/human/SeqWeaver_predictions/3pUTR_avg.txt.gz
  require_file data/raw/saluki_datapack/human/DeepRiPe_predictions_v83/5pUTR_avg.txt.gz
  require_file data/raw/saluki_datapack/human/DeepRiPe_predictions_v83/ORF_avg.txt.gz
  require_file data/raw/saluki_datapack/human/DeepRiPe_predictions_v83/3pUTR_avg.txt.gz

  if [ ! -s results/human_sequence_features_regions3mer4mer.tsv.gz ]; then
    echo "[reproduce] build human sequence features from Ensembl release 115"
    python -m mrna_half_life_paper.sequence_baseline \
      --species human \
      --skip-benchmark \
      --skip-note
  fi
  require_file results/human_sequence_features_regions3mer4mer.tsv.gz

  echo "[reproduce] preliminary XGBoost parameter scan and provenance"
  python scripts/reproduce_xgb_parameter_scan.py --device "${DEVICE}"

  echo "[reproduce] fixed-target global prior benchmark"
  python -m mrna_half_life_paper.global_mouse_prior --skip-note
  python -m mrna_half_life_paper.fold_robustness \
    --random-states 13 42 101 \
    --cv-splits 10 \
    --results-root results/tenfold_global_prior
  python -m mrna_half_life_paper.prior_missingness_analysis \
    --random-states 13 42 101 \
    --cv-splits 10 \
    --results-root results/prior_missingness_analysis
  python -m mrna_half_life_paper.prior_residual_analysis \
    --random-state 42 \
    --cv-splits 5 \
    --nuisance-cv-splits 5 \
    --device "${DEVICE}" \
    --results-root results/prior_residual_analysis

  echo "[reproduce] 0.10 ortholog-informed shrinkage target"
  python -m mrna_half_life_paper.ortholog_regularized_label \
    --sources reconstructed_mouse_pc1 \
    --lambdas 0.1 \
    --include-saluki-baseline \
    --include-controls \
    --random-states 13 42 101 \
    --cv-splits 10 \
    --device "${DEVICE}" \
    --results-root results/ortholog_regularized_label_10fold_main
  python -m mrna_half_life_paper.ortholog_regularized_label \
    --sources reconstructed_mouse_pc1 \
    --lambdas 0.05 0.1 0.3 \
    --include-controls \
    --random-states 13 42 101 \
    --device "${DEVICE}" \
    --results-root results/ortholog_regularized_label_multiseed
  python -m mrna_half_life_paper.ortholog_prior_ablation \
    --source reconstructed_mouse_pc1 \
    --lambdas 0.1 \
    --random-states 13 42 101 \
    --cv-splits 10 \
    --saluki-prior-baseline 0.839426 \
    --saluki-pure-baseline 0.754235 \
    --device "${DEVICE}" \
    --results-root results/ortholog_prior_ablation_lambda0p1_10fold
  python -m mrna_half_life_paper.ortholog_regularized_cross_target \
    --random-states 13 42 101 \
    --cv-splits 10 \
    --device "${DEVICE}" \
    --bootstrap-iterations 2000 \
    --results-root results/ortholog_regularized_cross_target
  python scripts/bootstrap_bmb_key_claims.py
  python scripts/audit_oof_integrity.py
}

run_figures() {
  echo "[reproduce] active BMB figures and documents"
  python scripts/redraw_bmb_flow_figures.py
  python scripts/redraw_bmb_dataset_summary_cn.py
  python scripts/redraw_bmb_human_pca_panel.py
  python scripts/redraw_bmb_study_influence_cn.py
  python scripts/redraw_bmb_ortholog_combined_cn.py
  python scripts/plot_bmb_ortholog_regularized_target.py
  python scripts/redraw_bmb_sequence_progression_cn.py
  python scripts/redraw_bmb_prior_summary_cn.py
  python scripts/redraw_bmb_supplementary_diagnostics.py
  bash scripts/preflight_bmb_release.sh
}

case "${STAGE}" in
  labels)
    run_labels
    ;;
  models)
    run_models
    ;;
  figures)
    run_figures
    ;;
  all)
    run_labels
    run_models
    run_figures
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac

echo "[reproduce] completed stage: ${STAGE}"
