# BMB final pre-submission audit

Date: 2026-07-17

## Overall assessment

The active English manuscript, supplementary information, title page, cover letter, declarations, reference library, figure set, and reproduction paths were reviewed as one submission chain. No unresolved internal contradiction was found in the reported estimands, sample universes, fold structure, equation numbering, figure/table numbering, author metadata, or active reference lists. The package is technically ready for submission after the remaining author confirmations listed below.

## Corrections made in this audit

- Separated the complete 12,307-gene shrinkage-target universe from the 10,768 mapped one-to-one-ortholog universe. The former supports `r = 0.9982`; the latter supports `r = 0.9979`, `RMSE = 0.065`, and `MAE = 0.050`.
- Distinguished the dynamic-coverage Fig. 4 agreement gain (`+0.0349`) from the fixed 13,265-gene paired estimate used for bootstrap inference (`+0.0314`).
- Rewrote the two-stage residual procedure to state explicitly that both RidgeCV and XGBoost are trained within each outer fold and that no held-out human target enters residual construction.
- Corrected Fig. 5d captions to describe both Pearson and Spearman gains.
- Added the exact sample universes and fold counts to Fig. 6 and its captions.
- Replaced unsupported visual claims about PCA spread with qualitative, separately scaled interpretation.
- Recast cross-target results as “no reduction detected” and explicitly stated that no formal equivalence claim is made because no equivalence margin was prespecified.
- Updated feature-method citations for CWCS, SeqWeaver, and DeepRiPe, and updated Ensembl release 115 to the Ensembl 2025 reference.
- Aligned the internal CRediT table with the formal author-contribution statements.
- Synchronized repository URL, release-candidate version, abstract count, and active reference counts across the manuscript-support files.

## Quantitative and reproducibility checks

- Full preflight: `158 checks`, `0 failures`, `1 warning`.
- Unit tests: `3 passed`.
- OOF integrity: nine global fixed-target files cover 12,916 genes; three cross-target files and 18 shrinkage-target files cover 12,307 genes; residual outputs cover 11,107 genes across folds 1-5.
- English abstract: 232 words; six keywords.
- Main figures: Fig. 1-8; supplementary figures: S1-S4. All referenced PNG files report 600 dpi, and matching lossless 600-dpi TIFF files are present.
- Main equations: (1)-(9); supplementary equations: (S1)-(S2). All are rendered as editable OMML in DOCX.
- Main tables: 1-3; supplementary tables: S1-S16.
- Active references: 46 in the main manuscript, 11 in the supplement, 47 unique records in the BibTeX/RIS union. All 45 DOI records resolve through Crossref; the XGBoost full title was checked against the ACM record because Crossref returns an abbreviated title. The two no-DOI records were checked against the official NeurIPS and JMLR pages.
- Quantitative claims for label audit, ortholog concordance, prior permutation, cross-target evaluation, shrinkage geometry, residual decomposition, and study-influence sensitivity match the versioned result tables.

## Remaining author-controlled gates

1. Both authors must confirm author order, contribution wording, shared affiliation, no-funding statement, and approval of the submitted version.
2. Wenzhuo Wang's ORCID should be added if available; it is optional if no ORCID is held.
3. A software license should be selected for downstream reuse. This does not block manuscript submission but affects repository reuse terms.
4. After the manuscript is closed in WPS, remove the temporary Office lock files. Generate and inspect the final supplementary PDF and upload ZIP only after this source state is frozen.

## Release state

The intended immutable repository tag for this audited state is `mRNA-PC1-label-v1.3`. PDF and ZIP generation remain intentionally deferred until the final author-controlled submission freeze.
