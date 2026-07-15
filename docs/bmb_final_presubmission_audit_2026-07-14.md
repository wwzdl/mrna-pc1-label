# BMB Final Pre-submission Audit

Date: 2026-07-14; release status updated 2026-07-15  
Target: *Bulletin of Mathematical Biology*  
Scope: Chinese/English main manuscripts, Supplementary Information, figures, tables, references, submission materials, code, public inputs, and reproducibility

## Executive verdict

The scientific manuscript is ready for author approval and journal submission. No blocking inconsistency was found in the three-estimand logic, canonical analysis universes, quantitative claims, figures, tables, citations, or OOF implementation. The package is **conditionally upload-ready**, because three actions remain external to the scientific source: author metadata/contribution confirmation, software-license selection, and final generation plus visual inspection of the English Supplementary PDF and upload ZIP. The public GitHub release gate was completed on 2026-07-15.

The 0.10 ortholog-regularized target remains a core contribution. It is not demoted to an incidental supplement, but it is consistently separated from the fixed Saluki-human-PC1 ranking. This framing preserves the result's novelty while avoiding the claim that changing a target improves prediction of an unchanged ground truth.

## Reviewer 1: scientific logic and positioning

### Strengths

- The paper has one coherent question: how study heterogeneity, external ortholog information, and weak target regularization change the interpretation of mammalian mRNA half-life benchmarks.
- The three claims are separable and independently testable: reference-free label influence screening; fixed-target human-only versus cross-species transfer; and a human-dominant ortholog-regularized target.
- `Gejman` is ranked by PC1 stability without Saluki labels or downstream model scores. Saluki agreement and ortholog concordance are validation axes rather than selection criteria.
- The 0.10 target is located quantitatively relative to human, mouse, and ordinary study-level variability, then checked against the original human target by cross-target evaluation.
- The Discussion and Limitations avoid claiming a new experimental ground truth, sequence-only superiority over Saluki, or mechanistic causality from predictive priors.

### Residual scientific risks already disclosed

- Label auditing uses one public compendium and lacks an independent second half-life resource.
- Leave-one-study-out influence combines sample count and profile direction; without a size-matched null, `Gejman` is a dominant influence rather than a sample-count-adjusted statistical outlier.
- `lambda = 0.10` is a conservative empirical choice whose portability requires external validation.
- The benchmark uses gene-level random OOF splits and does not recreate Saluki's fixed test split; cross-paper score comparisons therefore remain descriptive.

## Reviewer 2: statistics, leakage, and result integrity

### Findings

- No target or `gene_id` column enters the model feature matrices.
- XGBoost early stopping uses a validation subset drawn only from each outer training fold; the held-out fold is used only for prediction.
- Mouse priors are constructed from mouse measurements and Ensembl mapping before human cross-validation. They are correctly framed as external covariates in a cross-species transfer setting.
- Prior permutation is fold-aware and shuffles values only among genes with an available prior while preserving availability/high-confidence indicators. The control therefore disrupts gene-prior identity without changing column count or missingness structure.
- Residual decomposition is outer-fold safe: held-out human labels do not enter either stage. Its 12.8% quantity is appropriately presented as descriptive remaining-variance decomposition, not a new sequence-only benchmark.
- OOF integrity checks passed for 9 global fixed-target files at 12,916 genes, 3 cross-target files at 12,307 genes, 18 ortholog-regularized files at 12,307 genes, and the 11,107-gene residual subset. IDs are unique, universes are shared within each comparison, folds are complete, and predictions contain no missing values.
- Main-text bootstrap values, cross-target values, target-distance metrics, and residual metrics match their result tables.

### Interpretation boundary

The prior-enhanced result is not information leakage under the declared estimand, because mouse-side gene covariates are part of the allowed prediction input. It would be leakage only if reported as human sequence-only prediction or if held-out human target values entered prior construction; neither occurs here.

## Reviewer 3: journal format, figures, references, and package

### Passed checks

- English abstract: 238 words; keywords: 6.
- Author-year citations and unnumbered alphabetized reference lists match BMB style.
- Main manuscript: 8 contiguous figures and 3 captioned tables.
- Supplement: 4 contiguous figures and 15 captioned tables, with journal/title/authors/affiliations/corresponding-author metadata.
- All 12 referenced PNGs report 600 dpi; all corresponding TIFFs report 600 dpi and lossless LZW compression. No duplicate main/SI image is referenced.
- Current visual review found no cropping, omitted panel label, obvious overlap, unreadable annotation, or caption-image mismatch.
- The English main DOCX contains continuous line numbering and page fields; the English Supplementary DOCX contains page fields. All Chinese/English DOCX files are newer than their Markdown and figure sources.
- Main reference list: 44 entries; Supplementary reference list: 10 entries; unique union: 45. The BibTeX and RIS libraries contain exactly these 45 records. All 43 DOI records resolved and matched basic bibliographic metadata; the two no-DOI records were checked against official proceedings/journal pages.
- Automated package audit: 136 checks, 0 failures, 0 warnings in `BMB_SKIP_PDF=1` mode.

## Reproducibility and GitHub audit

- `requirements.txt` now includes the document-rendering dependencies actually imported by the code; `requirements-validated.txt` records the tested direct versions.
- Python compile, shell syntax, dependency consistency (`pip check`), and project smoke tests passed.
- MOESM2 and MOESM3 downloads are SHA-256 verified. The Ensembl homology URL is pinned to release 115 and its 106-MB source file is SHA-256 verified.
- The `models` stage now rebuilds the ignored human sequence-feature table when absent and checks all required Saluki regulatory blocks before training.
- Generated manifests and summaries use repository-relative paths; no local absolute path, private key, or GitHub/API token pattern was found in the release candidate.
- Legacy Saluki result-schema tokens (`*_vs_author` and `mouse_author`) were removed from retained result snapshots and normalized to the same `*_vs_saluki` and `saluki_mouse_prior` vocabulary used by the current code and manuscripts; numerical values were unchanged.
- The public standalone repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label) contains 277 source/result assets totaling 24.60 MiB in its initial commit, with a largest file of 2.29 MB. Raw data, large matrices, per-gene OOF arrays, feature packages, DOCX, PDF, ZIP, and TIFF submission artifacts are excluded.
- `CITATION.cff` is present. A legal software license is intentionally not chosen on the authors' behalf.

## Required release gates

1. Confirm Wenzhuo Wang and Yuebin Zhang's ORCID status and approve the final Author Contributions, author order, funding, and affiliations.
2. Select and add a software license.
3. Generate the frozen English Online Resource 1 PDF, inspect every page, assemble the upload candidate, and create one final ZIP with integrity and SHA-256 checks.

## Verification limitation

This final audit verified source code, existing OOF predictions, summaries, inputs, documents, and figures. It did not repeat the full GPU training suite from raw inputs, because that is a long-running reproduction job; the staged command and integrity checks are prepared for that independent rerun.
