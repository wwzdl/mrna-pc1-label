# mrna-pc1-label

Reproducible code and manuscript sources for a study-aware mammalian mRNA half-life benchmark. The project separates three estimands that should not be placed on one leaderboard:

1. consensus-label auditing across studies;
2. human-only prediction versus cross-species transfer under a fixed Saluki human PC1 target;
3. a human-dominant target defined by weak ortholog-informed shrinkage at `lambda = 0.10`.

The active manuscript targets *Bulletin of Mathematical Biology*. English files are submission sources; Chinese files are retained for author-side checking.

## Main Results

- A Saluki-label-independent leave-one-study-out screen ranks `Gejman` first in the primary sample-weighted analysis. Its geometric displacement is compatible with 500 conditional deletions of 15 non-Gejman samples and its rank falls to third under study balancing, whereas those comparator deletions do not reproduce its positive Saluki-agreement and ortholog-concordance changes.
- On 12,916 human genes, the repeated 10-fold by 3-seed human-only model reaches Pearson `0.748 +/- 0.001`; cross-species transfer with two mouse ortholog priors reaches `0.830 +/- 0.001`; partition-wise prior permutation within each outer fold returns to `0.748 +/- 0.001`.
- Across the 12,307 model-eligible genes, the `0.10 ortholog-informed shrinkage target` remains nearly identical to human no-Gejman PC1 (`r = 0.9982`); in a separately defined set of 10,768 mapped one-to-one ortholog pairs, its shift RMSE is `0.065`. Cross-target evaluation detects no reduction in original-human-label predictability but is not a formal equivalence test.
- In a strict target-vector ablation, exact target-construction mouse PC1 plus released Saluki mouse PC1 reaches `0.8537 +/- 0.0004`; removing the exact vector while retaining the Saluki prior and all indicators gives `0.8524 +/- 0.0005`. Direct reuse of the mouse vector contributing 10% of the target therefore accounts for an increment of `0.0013`.

See `docs/status_report.md` for the current claim boundaries and complete key values.

## Analysis Universes

The manuscript uses 12 fixed-definition core sets in three branches, plus eight derived coverage/sensitivity sets. They are task-specific eligibility sets, not successive attrition from one cohort. The two primary model training/evaluation sets are P2 (fixed-target, 12,916 genes) and S2 (shrinkage target, 12,307 genes).

| Analysis | Size | Role |
|:--|--:|:--|
| Global fixed-target training/evaluation | 12,916 human genes | Main human-only and cross-species-transfer benchmark |
| Both-priors subset | 11,107 human genes | Prior coverage and residual decomposition only |
| Shrinkage target construction | 12,644 human genes | Defines the 0.10 ortholog-informed target before model eligibility filters |
| Ortholog-informed shrinkage prediction | 12,307 human genes | Shrinkage-target training/evaluation |
| One-to-one label geometry | 10,768 ortholog pairs | Human/mouse target-distance analysis |
| Ortholog concordance | 12,592 ortholog pairs | Cross-species label comparison |

The 11,107-gene both-priors subset is not the global training universe. The 12,307-gene prediction set and 10,768-pair geometry set overlap at 10,682 human genes but are not nested. The complete 12-core/8-derived definition is in Supplementary Table S6; `results/bmb_analysis_universe_ledger.tsv` is the machine-readable ledger, and `python3 scripts/audit_analysis_universes.py --check` validates it from the source/result tables.

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

`environment.yml` provides a Python 3.11 Conda alternative. `requirements-validated.txt` records the direct package versions that passed the 2026-07-15 pre-submission audit; `requirements.txt` remains the portable installation entry point. The published tree-model benchmarks use XGBoost/CUDA; label reconstruction, auditing, document generation, and most figure generation do not require a GPU.

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

The `labels` stage also reproduces the dynamic/fixed-universe leave-one-study-out analyses, the 500-replicate conditional same-size deletion analysis, preprocessing sensitivity, and two study-balanced estimators. The `models` stage builds the human sequence feature table from Ensembl release 115 when it is absent, then reproduces the preliminary XGBoost parameter scan, global prior benchmark, 10-fold robustness, missing-prior analysis, cross-fitted two-stage residual decomposition, `lambda = 0.10` target benchmark, lambda sensitivity, prior ablation, cross-target evaluation, and paired-bootstrap summaries. Both stages write to the result paths used by the manuscripts.

After model reproduction, `scripts/audit_oof_integrity.py` verifies the expected evaluation-universe size, unique gene coverage, shared gene order across settings/seeds, fold coverage, and absence of missing predictions. `scripts/audit_analysis_universes.py --check` independently verifies the named core sets, derived coverage subsets, and the 10,682-gene S2/S3 overlap.

## Rebuild Submission Files

Rebuild all active figures, DOCX/PDF files, and the upload-candidate directory:

```bash
bash scripts/reproduce_bmb_key_results.sh figures
```

Run only the release preflight when result tables and figures are already current:

```bash
bash scripts/preflight_bmb_release.sh
```

The preflight runs unit tests and OOF-integrity checks, then checks manuscript/SI consistency, quantitative claims against result tables, figure numbering and resolution, references, BibTeX/RIS identity, DOCX line/page fields, SI metadata, and upload-candidate assembly.

During editing, rebuild and audit DOCX files without regenerating PDF or assembling a mixed-version upload candidate:

```bash
BMB_SKIP_PDF=1 bash scripts/preflight_bmb_release.sh
```

## Interpretation Boundaries

- Mouse priors are external mouse-side covariates fixed before human cross-validation; they are not part of the Saluki human PC1 target.
- Prior-enhanced results are cross-species transfer, not performance available from human sequence alone.
- Ortholog-informed shrinkage changes the target definition; its scores are not ranked against fixed-target results as if the estimands were identical.
- The base-sequence model can accept a newly partitioned 5'UTR/CDS/3'UTR sequence. `compact_all` requires precomputed regulatory blocks, and the transfer model additionally requires mappable mouse ortholog priors.

## Manuscript Sources

- `manuscript/bmb_submission/bmb_main_manuscript_en.md`
- `manuscript/bmb_submission/bmb_supplementary_material_en.md`
- `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
- `manuscript/bmb_submission/bmb_supplementary_material_cn.md`

All active main figures (Fig. 1-8) and supplementary figures (S1-S4) are generated from scripts and fixed result tables. Current PNG/SVG assets are kept in the source release; 600-dpi TIFF files and final PDFs are generated for the frozen submission package. The submission chain does not use cropped or AI-generated artwork as final figures.

## Release

The public source repository is [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label). The audited two-author manuscript state is tagged as `mRNA-PC1-label-v1.4.3`; run the release preflight before deriving any later submission package from a newer commit.

中文说明：本 README 采用英文以便审稿人与读者直接复现；中文定稿状态和术语说明见 `docs/status_report.md`。
