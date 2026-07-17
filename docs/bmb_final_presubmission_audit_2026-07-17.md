# BMB final pre-submission audit

Date: 2026-07-17

## Overall assessment

The active English manuscript, supplementary information, title page, cover letter, declarations, reference library, figure set, and reproduction paths were reviewed as one submission chain. No unresolved internal contradiction was found in the reported estimands, sample universes, fold structure, equation numbering, figure/table numbering, author metadata, or active reference lists. The package is technically ready for submission after the remaining author confirmations listed below.

## Corrections made in this audit

- Separated the 12,307 model-eligible shrinkage-prediction universe from the independently filtered 10,768 mapped one-to-one-ortholog pairs. The former supports `r = 0.9982`; the latter supports `r = 0.9979`, `RMSE = 0.065`, and `MAE = 0.050`.
- Replaced the ambiguous single-cohort wording with a twelve-set analysis-universe ledger. The main text now identifies P2 (12,916 genes) and S2 (12,307 genes) as the two primary modeling universes and records the non-nested S2/S3 overlap of 10,682 genes.
- Distinguished the dynamic-coverage Fig. 4 agreement gain (`+0.0349`) from the fixed 13,265-gene paired estimate used for bootstrap inference (`+0.0314`).
- Rewrote the two-stage residual procedure to state explicitly that both RidgeCV and XGBoost are trained within each outer fold and that no held-out human target enters residual construction.
- Corrected Fig. 5d captions to describe both Pearson and Spearman gains.
- Added the exact sample universes and fold counts to Fig. 6 and its captions.
- Replaced unsupported visual claims about PCA spread with qualitative, separately scaled interpretation.
- Recast cross-target results as “no reduction detected” and explicitly stated that no formal equivalence claim is made because no equivalence margin was prespecified.
- Updated feature-method citations for CWCS, SeqWeaver, and DeepRiPe, and updated Ensembl release 115 to the Ensembl 2025 reference.
- Aligned the internal CRediT table with the formal author-contribution statements.
- Synchronized repository URL, release-candidate version, abstract count, and active reference counts across the manuscript-support files.
- Replaced machine identifiers in reader-facing supplementary tables with concise display labels and repaired landscape-section pagination so each wide table remains with its caption.
- Clarified the roles of MOESM2 and MOESM3: MOESM2 is the reconstruction/audit input, Saluki human PC1 is the fixed supervision target, released Saluki mouse PC1 is an external prior, and MOESM3 processed sample columns are secondary controls.
- Distinguished the 13,921-gene full-reconstruction comparison (`r = 0.948`) from the fixed 13,265-gene paired Gejman-inference universe (`r = 0.9515`).
- Normalized manuscript typography so study names, dataset names, statistics, and explanatory expressions use body or mathematical fonts; monospace is reserved for code fields, file paths, and exact parameter identifiers.
- Expanded the Fig. 5c correlation axis to avoid visually exaggerating small between-setting differences, and corrected reference-author metadata in the synchronized BibTeX/RIS libraries.
- Replaced the remaining ambiguous “fold-wise shuffling” wording with the implemented partition-wise permutation within each outer fold.
- Added the gene-level split boundary: paralog families were not grouped across folds, so gene-family-held-out generalization is not established.
- Synchronized the fixed DOCX table grid with the equation-cell widths so WPS allocates the intended 5.89-inch center cell; added vertical padding, row-split protection, a unit test, and submission-package checks for this layout.

## Quantitative and reproducibility checks

- Full preflight: `169 checks`, `0 failures`, `0 warnings`.
- Unit tests: `4 passed`.
- OOF integrity: nine global fixed-target files cover 12,916 genes; three cross-target files and 18 shrinkage-target files cover 12,307 genes; residual outputs cover 11,107 genes across folds 1-5.
- English abstract: 228 words; six keywords.
- Main figures: Fig. 1-8; supplementary figures: S1-S4. All referenced PNG files report 600 dpi, and matching lossless 600-dpi TIFF files are present.
- All 12 referenced figures were inspected at readable scale after the final redraw; no clipped labels, overlapping text, missing panel labels, or ambiguous numeric annotations remain. Fig. 1 now labels the conditional result explicitly as the geometry lower-tail proportion.
- Main equations: (1)-(9); supplementary equations: (S1)-(S2). All are rendered as editable OMML in DOCX and pass the synchronized full-width WPS-layout audit.
- Main tables: 1-4; supplementary tables: S1-S16.
- Active references: 46 in the main manuscript, 11 in the supplement, 47 unique records in the BibTeX/RIS union. All 45 DOI records resolve through Crossref; the XGBoost full title was checked against the ACM record because Crossref returns an abbreviated title. The two no-DOI records were checked against the official NeurIPS and JMLR pages.
- Quantitative claims for label audit, ortholog concordance, prior permutation, cross-target evaluation, shrinkage geometry, residual decomposition, and study-influence sensitivity match the versioned result tables.
- The intended release tree contains 299 files totaling 26.17 MiB; the largest tracked file is 2.29 MB. Raw inputs, large matrices, fold-level prediction arrays, DOCX, PDF, TIFF, ZIP, and local caches remain excluded.
- A redacted high-confidence credential scan found no private keys or GitHub, OpenAI, or AWS tokens in the 270 release text files or the 21.95 MB Git patch history examined before release.

## Remaining author-controlled gates

1. Both authors must confirm author order, contribution wording, shared affiliation, no-funding statement, and approval of the submitted version.
2. Wenzhuo Wang's ORCID should be added if available; it is optional if no ORCID is held.
3. A software license should be selected for downstream reuse. This does not block manuscript submission but affects repository reuse terms.
4. Generate and inspect the final supplementary PDF and upload ZIP only after this source state is frozen.

## Release state

The intended immutable repository tag for this audited state is `mRNA-PC1-label-v1.4.4`. PDF and ZIP generation remain intentionally deferred until the final author-controlled submission freeze.
