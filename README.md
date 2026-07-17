# mRNA-PC1 Label

Reproducible analysis code and compact result tables for study-aware mammalian
mRNA half-life label auditing, fixed-target prediction, cross-species transfer,
and weak ortholog-informed target shrinkage.

This public repository is intentionally limited to the material needed to
understand and reproduce the computational results. Manuscript drafts,
submission documents, internal review notes, generated figures, large
intermediate matrices, and intellectual-property application files are kept
outside the public source tree.

## Main Results

- A Saluki-label-independent leave-one-study-out screen ranks `Gejman` first in
  the primary sample-weighted analysis. Study-balanced estimators place it
  third, showing that sample-count leverage contributes to its geometric
  influence.
- On 12,916 human genes, repeated 10-fold by 3-seed evaluation gives Pearson
  `0.748 +/- 0.001` for the human-only compact model and `0.830 +/- 0.001` for
  cross-species transfer with two mouse ortholog priors. Partition-wise prior
  permutation returns to `0.748 +/- 0.001`.
- Across 12,307 model-eligible genes, the `0.10 ortholog-informed shrinkage
  target` remains nearly identical to human no-Gejman PC1 (`r = 0.9982`).
  Within 10,768 mapped one-to-one ortholog pairs, its shift RMSE is `0.065`.
- A strict target-vector ablation gives `0.8537 +/- 0.0004` with the exact
  target-construction mouse PC1 and `0.8524 +/- 0.0005` after removing that
  vector while retaining the released Saluki mouse PC1 and all indicators.

## Analysis Universes

The analysis uses task-specific eligibility sets rather than a single
successive attrition cohort.

| Analysis | Size | Role |
|:--|--:|:--|
| Global fixed-target training/evaluation | 12,916 human genes | Main human-only and cross-species-transfer benchmark |
| Both-priors subset | 11,107 human genes | Prior coverage and residual decomposition only |
| Shrinkage target construction | 12,644 human genes | Target definition before model-eligibility filters |
| Shrinkage-target prediction | 12,307 human genes | Target-shrinkage training/evaluation |
| One-to-one label geometry | 10,768 ortholog pairs | Human/mouse target-distance analysis |
| Ortholog concordance | 12,592 ortholog pairs | Cross-species label comparison |

The complete machine-readable ledger is
`results/bmb_analysis_universe_ledger.tsv`. The 11,107-gene both-priors subset
is not the global training universe.

## Repository Layout

- `src/mrna_half_life_paper/`: analysis modules
- `scripts/`: data preparation, result reproduction, scientific plotting, and
  result-integrity checks
- `results/`: compact result tables and summaries supporting the reported
  values
- `data/processed/`: small processed annotations and pinned ortholog mappings
- `metadata/`: sample metadata and reusable templates
- `tests/`: focused tests for high-risk analysis inputs and study-influence
  controls

Generated figures are written under `figures/paper/` and are ignored by Git.

## Environment

Create and activate a local virtual environment from the repository root:

```bash
bash scripts/create_env.sh
source scripts/activate_env.sh
bash scripts/smoke_test.sh
```

`environment.yml` provides a Python 3.11 Conda alternative.
`requirements-validated.txt` records the package versions used for the audited
analysis; `requirements.txt` is the portable installation entry point. The
published tree-model benchmarks use XGBoost/CUDA, while label reconstruction,
study auditing, and most scientific plotting can run without a GPU.

## Public Inputs

MOESM2 and MOESM3 are downloaded from the public supplementary files of
Agarwal and Kelley (2022):

```bash
bash scripts/prepare_real_data.sh
```

The downloader verifies SHA-256 checksums before extraction:

- MOESM2: `ddf9441d1fe85d091d4b7f0a444b8dcf3658aa1d4158f66948bff00d5406de98`
- MOESM3: `00c49a9869e982527b07ac27896340fa34a3188199551f386027a8e8eb4f2217`

The compact regulatory feature blocks use the public Saluki datapack
(`datasets.zip`, DOI
[10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)). Because the
archive is approximately 18.8 GB, it is not downloaded automatically. Stage
the required human files as:

```text
data/raw/saluki_datapack/human/
  seqFeatWithKmerFreqs.txt.gz
  CWCS.txt.gz
  SeqWeaver_predictions/*_avg.txt.gz
  DeepRiPe_predictions_v83/*_avg.txt.gz
```

Ensembl sequence and one-to-one ortholog resources are pinned to release 115
and cached under ignored `data/raw/` paths.

## Reproduce Results

The staged entry point separates CPU label analysis, CUDA modeling, and figure
generation:

```bash
bash scripts/reproduce_bmb_key_results.sh labels
bash scripts/reproduce_bmb_key_results.sh models
bash scripts/reproduce_bmb_key_results.sh figures
```

Run every stage only after the Saluki datapack and a CUDA-capable XGBoost
environment are available:

```bash
bash scripts/reproduce_bmb_key_results.sh all
```

The `labels` stage rebuilds consensus labels, study-influence analyses,
same-size deletion controls, study-balanced estimators, ortholog concordance,
label-distance results, and sample-level PCA diagnostics.

The `models` stage rebuilds the human feature table when absent and reproduces
the fixed-target benchmark, prior permutation, repeated 10-fold robustness,
missing-prior analysis, cross-fitted residual decomposition, target-shrinkage
benchmark, lambda sensitivity, prior ablation, cross-target evaluation, and
paired-bootstrap summaries.

The `figures` stage regenerates publication figures from fixed result tables
under the ignored `figures/paper/` directory.

## Validate a Checkout

Run the lightweight public preflight:

```bash
bash scripts/preflight_public_release.sh
```

It checks Python and shell syntax and runs the focused unit tests. Analysis-
universe and OOF-integrity checks run automatically when their regenerated
per-gene or fold-level inputs are present.

## Interpretation Boundaries

- Mouse priors are external mouse-side covariates fixed before human
  cross-validation; they are not part of the fixed Saluki human PC1 target.
- Prior-enhanced results represent cross-species transfer, not performance
  available from human sequence alone.
- Ortholog-informed shrinkage changes the target definition, so its prediction
  score is not ranked against fixed-target results as if the estimands were
  identical.
- The base-sequence model can accept newly partitioned 5'UTR/CDS/3'UTR
  sequence. The compact model additionally requires precomputed regulatory
  blocks, and the transfer model requires mappable mouse ortholog priors.

## Public Release Scope

Included:

- analysis source code and focused tests;
- environment specifications and citation metadata;
- compact summaries and result tables needed to audit reported values;
- small processed annotations and metadata needed to connect the workflow.

Excluded:

- manuscript and supplementary-material source files;
- Word/PDF rendering and journal-upload scripts;
- internal review logs and working notes;
- generated PNG/SVG/TIFF/PDF figures;
- raw downloads, feature packages, intermediate matrices, and per-gene OOF
  arrays;
- backups, submission archives, and patent/software-copyright materials.

The public repository is
[wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label). Citation
metadata are provided in `CITATION.cff`.
