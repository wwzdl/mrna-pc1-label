# mrna-pc1-label

Reproducible code and manuscript sources for a study-aware mammalian mRNA half-life benchmark. The project separates three estimands that should not be placed on one leaderboard:

1. consensus-label auditing across studies;
2. human-only prediction versus cross-species transfer under a fixed Saluki human PC1 target;
3. a human-dominant target with weak ortholog regularization at `lambda = 0.10`.

The active manuscript targets *Bulletin of Mathematical Biology*. English files are submission sources; Chinese files are retained for author-side checking.

## Main Results

- A reference-free leave-one-study-out stability screen identifies `Gejman` as the dominant influence on the human consensus-label geometry. Saluki agreement and human-mouse ortholog concordance are subsequent validation axes, not selection criteria.
- On 12,916 human genes, the repeated 10-fold by 3-seed human-only model reaches Pearson `0.748 +/- 0.001`; cross-species transfer with two mouse ortholog priors reaches `0.830 +/- 0.001`; fold-wise prior permutation returns to `0.748 +/- 0.001`.
- The `0.10 ortholog-regularized target` remains nearly identical to human no-Gejman PC1 (`r = 0.9982`, `RMSE = 0.065`) and shows no detectable loss of original-human-label predictability in cross-target evaluation.

See `docs/status_report.md` for the current claim boundaries and complete key values.

## Analysis Universes

| Analysis | Size | Role |
|:--|--:|:--|
| Global fixed-target training/evaluation | 12,916 human genes | Main human-only and cross-species-transfer benchmark |
| Both-priors subset | 11,107 human genes | Prior coverage and residual decomposition only |
| Ortholog-regularized target | 12,307 human genes | Regularized-target training/evaluation |
| One-to-one label geometry | 10,768 ortholog pairs | Human/mouse target-distance analysis |
| Ortholog validation | 12,592 ortholog pairs | Cross-species label validation |

The 11,107-gene both-priors subset is not the global training universe.

## Repository Layout

- `src/mrna_half_life_paper/`: analysis modules
- `scripts/`: data preparation, result reproduction, figure generation, document rendering, and submission audit
- `results/`: compact result tables and summaries
- `metadata/`: sample metadata and templates
- `manuscript/bmb_submission/`: active Chinese/English manuscripts, figure sources, and references
- `docs/`: current status, method notes, and chronological experiment log

Raw data, intermediate matrices, OOF prediction arrays, local submission exports, and 600-dpi TIFF bundles are intentionally excluded from the lightweight source repository.

## Environment

Create and activate a local virtual environment from the repository root:

```bash
bash scripts/create_env.sh
source scripts/activate_env.sh
bash scripts/smoke_test.sh
```

`environment.yml` provides a Python 3.11 Conda alternative. `requirements-validated.txt` records the direct package versions that passed the 2026-07-14 pre-submission audit; `requirements.txt` remains the portable installation entry point. The published tree-model benchmarks use XGBoost/CUDA; label reconstruction, auditing, document generation, and most figure generation do not require a GPU.

## Public Inputs

MOESM2 and MOESM3 are downloaded from the public supplementary files of Agarwal and Kelley (2022):

```bash
bash scripts/prepare_real_data.sh
```

The download script verifies SHA-256 checksums before extraction: MOESM2 `ddf9441d1fe85d091d4b7f0a444b8dcf3658aa1d4158f66948bff00d5406de98`; MOESM3 `00c49a9869e982527b07ac27896340fa34a3188199551f386027a8e8eb4f2217`.

The compact regulatory feature blocks use the public Saluki datapack (`datasets.zip`, DOI [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)). Because that archive is approximately 18.8 GB, it is not downloaded automatically. Stage the required human files as:

```text
data/raw/saluki_datapack/human/
  seqFeatWithKmerFreqs.txt.gz
  CWCS.txt.gz
  SeqWeaver_predictions/*_avg.txt.gz
  DeepRiPe_predictions_v83/*_avg.txt.gz
```

Ensembl sequence and ortholog resources are downloaded by the corresponding analysis modules and cached under ignored `data/raw/` paths.
The URLs are pinned to Ensembl release 115. The mouse homology source is verified against SHA-256 `ed82f17dda96df6b47d0885e07bb3adc7d274c760e10437fa3a58e7081206dd0`.

## Reproduce Key Results

The staged entry point keeps CPU label analysis, CUDA modeling, and figure/document generation separate:

```bash
bash scripts/reproduce_bmb_key_results.sh labels
bash scripts/reproduce_bmb_key_results.sh models
bash scripts/reproduce_bmb_key_results.sh figures
```

Run all stages only after the large Saluki datapack files and a CUDA-capable XGBoost environment are available:

```bash
bash scripts/reproduce_bmb_key_results.sh all
```

The `models` stage builds the human sequence feature table from Ensembl release 115 when it is absent, then reproduces the global prior benchmark, 10-fold robustness, missing-prior analysis, two-stage residual decomposition, `lambda = 0.10` target benchmark, lambda sensitivity, prior ablation, cross-target evaluation, and paired-bootstrap summaries. It writes to the result paths used by the manuscripts.

After model reproduction, `scripts/audit_oof_integrity.py` verifies the expected evaluation-universe size, unique gene coverage, shared gene order across settings/seeds, fold coverage, and absence of missing predictions.

## Rebuild Submission Files

Rebuild all active figures, DOCX/PDF files, and the upload-candidate directory:

```bash
bash scripts/reproduce_bmb_key_results.sh figures
```

Run only the release preflight when result tables and figures are already current:

```bash
bash scripts/preflight_bmb_release.sh
```

The preflight checks manuscript/SI consistency, quantitative claims against result tables, figure numbering and resolution, references, BibTeX/RIS identity, DOCX line/page fields, SI metadata, and upload-candidate assembly.

During editing, rebuild and audit DOCX files without regenerating PDF or assembling a mixed-version upload candidate:

```bash
BMB_SKIP_PDF=1 bash scripts/preflight_bmb_release.sh
```

## Interpretation Boundaries

- Mouse priors are external mouse-side covariates fixed before human cross-validation; they are not part of the Saluki human PC1 target.
- Prior-enhanced results are cross-species transfer, not performance available from human sequence alone.
- Ortholog regularization changes the target definition; its scores are not ranked against fixed-target results as if the estimands were identical.
- The base-sequence model can accept a newly partitioned 5'UTR/CDS/3'UTR sequence. `compact_all` requires precomputed regulatory blocks, and the transfer model additionally requires mappable mouse ortholog priors.

## Manuscript Sources

- `manuscript/bmb_submission/bmb_main_manuscript_en.md`
- `manuscript/bmb_submission/bmb_supplementary_material_en.md`
- `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
- `manuscript/bmb_submission/bmb_supplementary_material_cn.md`

All active main figures (Fig. 1-8) and supplementary figures (S1-S4) are generated from scripts and fixed result tables. Current PNG/SVG assets are kept in the source release; 600-dpi TIFF files and final PDFs are generated for the frozen submission package. The submission chain does not use cropped or AI-generated artwork as final figures.

## Release

Before public release, create the GitHub repository `wwzdl/mrna-pc1-label`, run the preflight, and tag the reviewed state as `mRNA-PC1-label-v1.0`. Update the manuscript Code Availability statement from future to present tense only after the repository and release are public.

中文说明：本 README 采用英文以便审稿人与读者直接复现；中文定稿状态和术语说明见 `docs/status_report.md`。
