# Supplementary Information: Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction

**Journal:** *Bulletin of Mathematical Biology*  
**Article title:** Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction  
**Authors:** Wenzhuo Wang, Ying Shao  
**Affiliation:** School of Science, Dalian Maritime University, Dalian 116026, China  
**Corresponding author:** Ying Shao; Email: yshao@dlmu.edu.cn
**ORCID:** Ying Shao: https://orcid.org/0000-0002-4056-5757

## S1. Supplementary overview

This Supplementary Information supports the three claims in the main text. First, study-aware label auditing uses Saluki-label-independent PC1 stability to screen influential studies, separates sample-count leverage from directional comparison gains, and then evaluates the direction of change through Saluki agreement and ortholog concordance. Second, with Saluki human PC1 fixed, the human-only baseline and cross-species transfer setting are evaluated separately and supported by permutation, paired bootstrap, residual decomposition, and multi-seed checks. Third, weak ortholog-informed target shrinkage at 0.10 defines a human-dominant mammalian stability target and is bounded by label distance, study-noise scale, cross-target evaluation, and prior ablation. This target is a core label-construction result, but it does not enter the fixed Saluki human PC1 ranking.

Compared with the main text, the Supplementary Information preserves the methodological details required to support the conclusions, including sample-level PCA definitions, complete study-influence paths, conditional same-size deletion and preprocessing sensitivity, prior-permutation and residual controls, ortholog label distances, training-universe boundaries, 10-fold split-sensitivity checks, and implementation notes (Agarwal and Kelley 2022; Chen and Guestrin 2016; Pedregosa et al. 2011).

## S2. Data extraction, preprocessing, and ortholog definitions

This section expands the fields, parameters, and definitions only briefly described in Section 3.1 of the main text. The source data are taken from the public supplement of the Saluki study, specifically the MOESM2 and MOESM3 tables (Agarwal and Kelley 2022). MOESM2 contains transformed human and mouse mRNA half-life gene by sample matrices. The first two columns are the Ensembl gene ID and gene name, and the remaining columns are transformed values from different studies, experimental methods, cell types, and biological replicates. This table is the primary input for label reconstruction and leave-one-study-out analysis. MOESM3 is the more heavily processed matrix released with Saluki. Its third column is the published "half-life (PC1)" summary label. **Saluki human PC1** denotes the released human label and is the fixed supervision target in the fixed-target benchmarks. **Released Saluki mouse PC1** denotes the released mouse label after ortholog mapping for use as an external mouse-side covariate; `mouse_pc1` instead denotes our locally reconstructed mouse prior. The human worksheet indicates that Saluki human PC1 was derived after removal of Gejman, whereas the mouse worksheet uses the full mouse dataset. The released PC1 columns therefore define the fixed target and published processing reference, while the processed sample columns support secondary controls. Neither is used to reconstruct the local labels from MOESM2. Because MOESM2 already contains transformed values, we do not apply an additional log transform and instead begin with sample-wise z-score, followed by iterative PCA missing-value imputation and quantile normalization (Tipping and Bishop 1999; Troyanskaya et al. 2001; Bolstad et al. 2003).

The precomputed sequence/regulatory feature blocks used by `compact_all` are from the associated public Saluki dataset (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)).

Ortholog mappings were obtained from the mouse-human homology table in Ensembl Compara release 115 (Dyer et al. 2025). We retained only records with `homology_type = ortholog_one2one` and `homology_species = homo_sapiens`; the interpretation of orthology follows established evolutionary-genomics usage (Gabaldon and Koonin 2013). The high-confidence ortholog subset was not redefined in this work; we directly used the binary Ensembl field `is_high_confidence`. According to the Ensembl orthology QC documentation, this released flag is assigned from identity together with gene-order-conservation or whole-genome-alignment thresholds when those data are available, with a tree-compliance fallback otherwise. We did not reconstruct the upstream QC and used `is_high_confidence == 1` only to define a stricter sensitivity subset.

The main text has already reported the global agreement between our reconstructed labels and the Saluki labels. On the 13,921 genes shared by the full human reconstruction and Saluki human PC1, Pearson is 0.948; the corresponding reconstructed-mouse value is 0.958. The fixed A2 universe used for paired Gejman inference additionally requires a no-Gejman label and therefore contains 13,265 genes, on which the full-label baseline is 0.9515. We call these minimum-three-observation vectors the **audit-threshold** reconstructions. Target construction instead uses the **target-construction human no-Gejman PC1** (minimum 10 observations) and **target-construction reconstructed mouse PC1** (minimum 5 observations). Coverage-threshold sensitivity is small but nonzero: audit- and target-threshold mouse PC1 values overlap on 11,882 genes (Pearson = 0.99907; z-scale RMSE = 0.04312), and the corresponding human no-Gejman values overlap on 12,644 genes (Pearson = 0.99827; z-scale RMSE = 0.05878). The names are therefore kept distinct even though their global geometry is similar. Exact reconstruction parameters are summarized in Supplementary Table S1, and modeling-input boundaries are summarized in Supplementary Table S2.

### Supplementary Table S1. Reproducibility parameters for label reconstruction and leave-one-study-out analysis

| Step | Parameter / definition | Output | Notes |
|:--|:--|:--|:--|
| Input matrix | MOESM2 transformed human/mouse half-life matrix | species-specific gene by sample matrix | No additional log transform; starts from transformed values |
| Metadata parsing | `study`, `method`, `cell type`, `replicate`, `species` | standardized sample metadata | Sample names were parsed and manually checked for study-level leave-one-out analysis |
| Audit-threshold coverage | `min_observed_per_gene=3` | audit-threshold human/mouse PC1 vectors | Used for full-label auditing, leave-one-study-out analysis, and the reconstructed prior in the fixed-target transfer branch |
| Target-construction coverage | human `min_observed_per_gene=10`; mouse `min_observed_per_gene=5` | target-construction human no-Gejman and reconstructed mouse PC1 vectors | Used for 0.10 ortholog-informed target shrinkage |
| Sample standardization | sample-wise z-score | standardized matrix | Reduces scale differences across experiments |
| Missing-value handling | iterative PCA imputation; up to 30 iterations; tolerance `1e-6` | complete matrix | Avoids distorting gene-level covariance by simple mean filling |
| Distribution alignment | quantile normalization | processed matrix | Reduces the effect of sample-distribution differences on PC1 |
| Consensus extraction | PCA on the gene by sample matrix; PC1 sign aligned to the gene-wise mean | reconstructed PC1 label | human full vs Saluki human PC1 Pearson = 0.948; mouse = 0.958 |
| Study influence | full rerun of the above pipeline after leave-one-study-out removal | PC1-stability screen and Saluki-agreement validation | Saluki-label-independent stability screens studies; Saluki agreement validates; downstream scores are not used |

### Supplementary Table S2. Modeling inputs, applicability, and claim boundaries

| Model tier | Required inputs | Can it be used for a new sequence? | Main output | Claims that should not be made |
|:--|:--|:--|:--|:--|
| Base sequence model | partitioned 5'UTR, CDS, and 3'UTR sequences | Yes | base-sequence human-only prediction | Should not be presented as reaching `compact_all` or prior-enhanced peak performance |
| Compact human-only model | Saluki-datapack Ensembl `gene_id` + base features + precomputed CWCS / SeqWeaver / DeepRiPe blocks | Only after the regulatory blocks are regenerated | pure human sequence/regulatory benchmark | Should be described as a transparent human-only baseline, not as a sequence-only win claim |
| Global prior-enhanced model | `compact_all` + two mouse-prior values + availability / high-confidence binary indicators | Requires mappable mouse ortholog priors | cross-species transfer using mouse-side covariates external to the human target | Should not be described as sequence-only or raw-sequence end-to-end prediction |
| Coverage-aware fallback | prior coverage status + human-only prediction | Useful for tiered deployment | use prior-enhanced when priors exist; otherwise fall back to human-only | Should not extrapolate prior-enhanced conclusions to genes with missing priors |
| 0.10 ortholog-informed target shrinkage | target-construction human no-Gejman PC1 + mapped target-construction mouse PC1 + $\lambda$ | Not a standalone predictor for arbitrary new sequences | human-dominant mammalian stability target | Should not be presented as direct human half-life or ranked with the fixed target |

![](figures/supplement/FigS01_prediction_boundary_flow.png)

**Supplementary Fig. S1** Prediction-use boundaries under different input conditions. a, When only partitioned 5'UTR/CDS/3'UTR sequences are available, only the base-sequence human-only model can be used. b, Once an Ensembl human gene and precomputed regulatory blocks are available, the compact sequence/regulatory model can be used. c, If mappable mouse ortholog priors are also available, the prior-enhanced model can be used in a cross-species transfer setting whose mouse-side covariates are external to the human target. d, In deployment, the model tier should be chosen according to prior coverage; genes without mouse priors should fall back to the compact human-only or base-sequence tier

## S3. Sample-level PCA diagnosis

Sample-level PCA is the most direct visual diagnosis before leave-one-study-out analysis. The main text already retains the human raw/processed PCA as the primary evidence in Fig. 3, so we do not repeat the same human panel here. Instead, we provide the algorithmic definition and the mouse-side result. In the human raw matrix, sample PCA shows that samples cluster mainly by study and method rather than by cell type, consistent with the key observation reported in the original paper (Agarwal and Kelley 2022). Preprocessing changes the PC geometry, but study-of-origin structure remains. The mouse views also retain study-of-origin structure. Because the human and mouse panels use separately fitted PCA axes and different numerical scales, their absolute spread is not compared.

This sample-level PCA is used only for visualization and diagnosis. It is not used to construct the supervised labels and does not participate in downstream model training. The difference between the raw and processed views is as follows. The raw view preserves the original scale of the MOESM2 transformed matrix as much as possible and fills missing entries only because PCA cannot be computed directly on incomplete data. The processed view follows the same preprocessing pipeline used for label reconstruction and is used to ask whether study/method structure remains after standardization and distribution alignment. The low-rank PCA imputation is not supervised by Saluki human PC1 or any downstream prediction label; it uses only low-rank structure internal to the gene by sample matrix to reconstruct missing positions.

The imputation algorithm is as follows. Let $X$ be the input matrix, with genes in rows and samples in columns, and let $\mathcal{M}$ denote the set of originally missing positions. We first retain genes with at least three observed samples. In step 1, missing entries are initialized with the mean of the corresponding sample column; if a column were entirely missing, it would be initialized to 0, although this case does not occur in the present data. In step 2, PCA is performed on the current filled matrix using `sklearn.decomposition.PCA(svd_solver="full", random_state=0)`, and the rank-$k$ low-rank approximation is reconstructed with `inverse_transform`. In step 3, only entries in $\mathcal{M}$ are replaced with their low-rank reconstructed values; observed entries remain fixed throughout the iteration. At iteration $t$, convergence is monitored by

```equation
number: S1
\delta_t = \max_{(i,j)\in\mathcal{M}} \left|X_{ij}^{(t)} - X_{ij}^{(t-1)}\right|
```

The algorithm stops early when $\delta_t<\mathrm{tol}$; otherwise, it continues to the maximum number of iterations. The raw view uses $k=5$, `max_iter=10`, and `tol=1e-5`. The processed view uses the label-reconstruction defaults, $k=\min\{5,\min(n_{\mathrm{genes}},n_{\mathrm{samples}})-1\}$, `max_iter=30`, and `tol=1e-6`. After imputation, the gene by sample matrix is transposed to a sample by gene matrix, and the first two sample PCs are computed.

The reproducibility script is:

```bash
bash scripts/run_sample_pca_imputation_reproducibility.sh --species both
```

The default output directory is `results/sample_pca_imputation_reproducibility/`. Each species subdirectory contains `raw_low_rank_imputed_matrix.tsv`, `processed_matrix.tsv`, `sample_pca_raw.tsv`, `sample_pca_processed.tsv`, `sample_pca_raw_vs_processed.png`, and `imputation_and_pca_parameters.json`. The parameter JSON records the input matrix, coverage filter, raw/processed imputation parameters, sample-PCA orientation, and the boundary note that no label or model information was used.

### Supplementary Table S3. Definitions of raw and processed matrices in sample-level PCA diagnosis

| PCA view | Input matrix | Missing-value handling | Additional preprocessing | PCA procedure | Purpose and boundary |
|:--|:--|:--|:--|:--|:--|
| Raw transformed view | MOESM2 transformed gene by sample matrix after `min_observed_per_gene=3` coverage filtering | Missing values filled only for visualization; iterative low-rank PCA imputation with `n_components=5`, `max_iter=10`, `tol=1e-5` | No sample-wise z-score, no quantile normalization, no PC1-label extraction | PCA computed on the first two PCs after transposing the imputed matrix to sample by gene | Visualizes study/method structure in the transformed matrix; not used as a supervised label or model input |
| Processed view | the same coverage-filtered MOESM2 matrix | iterative PCA imputation from the label-reconstruction pipeline; default `max_iter=30`, `tol=1e-6` | `apply_log=False`; sample-wise z-score and quantile normalization | PCA computed on the first two PCs after transposing the processed gene by sample matrix to sample by gene | Evaluates whether study/method structure remains after standardization, imputation, and distribution alignment; used to motivate leave-one-study-out auditing |

![](figures/supplement/FigS_mouse_pca_panel.png)

**Supplementary Fig. S2** Mouse raw and processed sample-level PCA. a, Raw transformed matrix after missing-value imputation for visualization only. b, Processed matrix after sample-wise standardization, iterative PCA imputation, and quantile normalization. Each point is one mouse half-life sample; color indicates experimental method and marker shape indicates study. Study-of-origin structure remains after preprocessing; absolute PC spread is not compared with main-text Fig. 3 because the species-specific PCA axes and scales were fitted separately

## S4. Supplementary leave-one-study-out rankings and result paths

The main text presents the primary sample-weighted human leave-one-study-out result in Fig. 4, so we retain only the full mouse ranking here. Complete human and mouse values are stored in `results/study_influence/human/study_influence.tsv` and `results/study_influence/mouse/study_influence.tsv`, respectively. In the primary human estimator, PC1 stability does not use Saluki PC1 and ranks Gejman first. Saluki-PC1 agreement gain then measures recovery of the published processing choice, and ortholog concordance supplies a cross-species comparison. Sample-count, fixed-universe, and preprocessing controls are reported separately in Section S5. No downstream model score entered this assessment. In mouse, removal of Hanna changes full-label geometry but does not improve agreement with the released Saluki mouse PC1, so it is not interpreted as a mouse analogue of the human case.

This human/mouse asymmetry shows that geometric influence does not automatically justify removal. For human, the primary geometric ranking is sample-count sensitive, but the positive Saluki and ortholog gains are not reproduced by conditional same-size deletion of other samples. The conclusion is a directional study effect superimposed on sample-count leverage. The mouse compendium lacks the same convergence of evidence and is kept intact.

![](figures/supplement/FigS_mouse_study_influence.png)

**Supplementary Fig. S3** Mouse leave-one-study-out full ranking. a, Points show the correlation between each leave-one-study-out PC1 and the full mouse PC1; lower values indicate greater influence on consensus-label geometry. The dashed line marks $r = 1$. b, Bars show the change in agreement with the released Saluki mouse PC1 after removing each study. Hanna has the lowest stability but a negative agreement change, so mouse does not show the same convergence of evidence observed for human Gejman in main-text Fig. 4

## S5. Sample-count, gene-universe, and preprocessing sensitivity of human study influence

The primary leave-one-study-out estimator gives every sample equal weight. Because Gejman contributes 15 of 54 human samples, its removal changes the sample count much more than removal of most other studies. We therefore drew 15 samples without replacement from the 39 non-Gejman samples, reran the full reconstruction pipeline, and repeated this procedure 500 times. Every comparator retained Gejman. The observed Gejman stability of 0.955 was not unusually low relative to this conditional deletion distribution (median 0.944; central 95% comparator range, 0.932-0.956; empirical lower-tail proportion = 0.964). Geometric influence alone is therefore confounded by the number of samples removed.

The two directional comparison gains behaved differently. The observed Saluki-agreement gain was +0.0314, whereas the comparator median was -0.0926 (central 95% comparator range, -0.1151 to -0.0713). The observed ortholog-concordance gain was +0.0147, whereas the comparator median was -0.0616 (central 95% comparator range, -0.0800 to -0.0476). No comparator reached either observed gain, giving the minimum attainable empirical upper-tail proportion, 1/501 = 0.002, for each. These values describe deletion sensitivity conditional on examining Gejman. They are not exchangeable whole-study p values, do not account for selecting Gejman, and are not adjusted for examining multiple studies. The result does not support a sample-count-adjusted geometric outlier; deleting the same number of other samples did not produce the same direction of change.

Gene-universe and preprocessing sensitivity gave a complementary result. Gejman ranked first in all 12 combinations of dynamic/fixed coverage, coverage thresholds, PCA-imputation ranks, and iterative-PCA/sample-median imputation, and both comparison gains remained positive. By contrast, when study sample counts were balanced either by $1/\sqrt{n_s}$ PCA weights or by collapsing each study to its mean profile, Gejman ranked third rather than first. These analyses show which part of the claim is robust: the primary ranking is stable across preprocessing choices but sensitive to study weighting, whereas the directions of the Saluki-processing and ortholog-comparison gains are stable.

![](figures/supplement/FigS_study_influence_sensitivity.png)

**Supplementary Fig. S4** Sample-size and pipeline sensitivity of human study influence. a, Conditional comparator distribution of PC1 stability from 500 full-pipeline reruns after removing 15 non-Gejman samples while retaining Gejman; lower values indicate more geometric change, and the observed Gejman value is not extreme in the lower tail. b, Conditional comparator distributions of Saluki-agreement and ortholog-concordance gains; orange points show the observed Gejman gains, both beyond all 500 comparator deletions. c, Both comparison gains remain positive across 12 dynamic/fixed gene-universe, coverage, rank, and imputation settings; Gejman ranks first in each primary-style setting. d, Study balancing changes the geometric rank from first to third under both estimators. Panel source data are stored in `manuscript/bmb_submission/figures/data/FigS04_panel_data.tsv`

### Supplementary Table S4. Sample-count and pipeline sensitivity of the Gejman study-influence result

| Analysis | Gejman result | Conditional comparator | Empirical tail proportion or rank | Interpretation |
|:--|:--|:--|:--|:--|
| Primary sample-weighted geometry | stability 0.9550 | median 0.9438; central 95% range 0.9319-0.9561 | lower-tail proportion = 0.964 | Geometry is not extreme after matching the number of removed samples |
| Saluki agreement gain | +0.0314 | median -0.0926; central 95% range -0.1151 to -0.0713 | upper-tail proportion = 1/501 = 0.002 | Conditional deletions do not reproduce the directional gain |
| Ortholog concordance gain | +0.0147 | median -0.0616; central 95% range -0.0800 to -0.0476 | upper-tail proportion = 1/501 = 0.002 | Conditional deletions do not reproduce the cross-species gain |
| Fixed gene universe, primary preprocessing | stability 0.9527; Saluki gain +0.0319; ortholog gain +0.0153 | 14,244 fixed-universe genes | rank 1/19 | Dynamic coverage does not create the directional result |
| Preprocessing sensitivity | both gains positive in all 12 settings | minimum observed samples 3/5/10; PCA rank 3/5/10; iterative PCA or sample median | rank 1/19 in all settings | Primary ranking is stable across preprocessing choices |
| Equal-study PCA weighting | stability 0.9796 | worst study Dieterich, stability 0.9653 | rank 3/19 | Sample count contributes to primary geometric rank |
| Study-mean collapse | stability 0.9776 | worst study Dieterich, stability 0.9634 | rank 3/19 | Independent study balancing gives the same qualification |

## S6. Additional explanation of permutation control

The main text reports the key permutation result: the global model with real double priors achieves Pearson = 0.829, whereas performance falls back to 0.745 after the priors are shuffled. The control isolates what is being permuted. First, within each outer fold, the fit, inner-validation, and held-out partitions are permuted separately, so a held-out gene-prior assignment cannot be restored from a training partition. Second, permutation occurs only among genes that originally have prior values; the identities of genes with missing priors, the corresponding availability indicators, and the high-confidence 0/1 indicators are kept fixed. Third, the control specifically disrupts the mapping between gene identity and prior value. The real priors improve Pearson by roughly +0.085 and R2 by roughly +0.135 relative to the shuffled priors, identifying gene-specific correspondence rather than extra feature columns as the source of the gain.

Across three independent 10-fold OOF evaluations, the real global model with both priors achieves Pearson 0.830 ± 0.001, whereas both the shuffled-prior control and human-only `compact_all` achieve 0.748 ± 0.001; the coverage-aware fallback achieves 0.832 ± 0.001. The ± values are across-seed SDs. Averaging the three seed-specific OOF predictions per gene forms a seed-averaged OOF ensemble. Conditional on that fitted ensemble and gene resampling, paired bootstrap gives a real-versus-shuffled gain of 0.0794 (95% CI, 0.0735-0.0855) and a real-versus-human-only gain of 0.0799 (95% CI, 0.0740-0.0862). These intervals do not include uncertainty from retraining or drawing new folds.

## S7. Additional explanation of residual decomposition

The goal of residual decomposition is to separate prior signal from human-feature signal. A direct correlation analysis already shows that `mouse_pc1` correlates with Saluki human PC1 at nearly 0.78, so the prior itself is a strong predictor. To avoid misreading this as evidence that the prior-enhanced model merely copies the prior, we explicitly decompose prediction into a two-stage out-of-fold procedure rather than putting priors and human features into a single black box and interpreting the result afterward.

The implementation is as follows. First, we retain only genes with both `mouse_pc1` and the released Saluki mouse PC1, i.e. the `both_priors_available` subset (N = 11,107), because prior-only decomposition requires both mouse-side prior values for each gene. Second, within each 5-fold outer split, a 5-fold inner cross-fitting loop fits RidgeCV on four inner folds and predicts the fifth, yielding nuisance-prior predictions for every outer-training gene without using that gene's target in its own nuisance fit. These inner-OOF predictions define the residual target as Saluki human PC1 minus nuisance-prior prediction. A separate RidgeCV fit on all outer-training genes predicts the outer held-out genes. Third, an XGBoost/CUDA model is trained only on outer-training human `compact_all` features and the cross-fitted residuals. Fourth, the final outer-held-out prediction is the sum of the held-out prior prediction and held-out compact-residual prediction. No outer-held-out human target enters either training stage.

The reproducibility command is:

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.prior_residual_analysis
```

### Supplementary Table S5. Two-stage OOF decomposition in prior residual analysis

| Model | Genes | Prior features | Human features | Pearson r | Spearman ρ | R² | ΔR² vs prior | Remaining variance explained |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| Prior only (linear) | 11,107 | 2 | 0 | 0.792 | 0.782 | 0.627 | 0.000 | 0.00% |
| Human compact only | 11,107 | 0 | 1,802 | 0.739 | 0.731 | 0.544 | -0.084 | NA |
| Prior plus human residual | 11,107 | 2 | 1,802 | 0.826 | 0.817 | 0.675 | 0.048 | 12.85% |

The reported residual fraction is computed as

```equation
number: S2
R_{\mathrm{remaining}}^2 = \frac{R_{\mathrm{final}}^2 - R_{\mathrm{prior}}^2}{1 - R_{\mathrm{prior}}^2}
```

Supplementary Eq. (S2) shows that, after the prior-only model has explained part of the variance, the human `compact_all` features explain about 12.9% of the remaining variance. The prior-enhanced model is therefore more accurately described as a strong mouse prior plus a human correction than as a setting in which the prior alone is sufficient.

These results also clarify the method boundary. Supervised labels for held-out human genes do not enter the input, and all predictions are generated out of fold. The mouse priors nevertheless constitute external gene-level information correlated with the target; the corresponding result is therefore a cross-species-transfer estimate rather than sequence-only prediction. The complete transfer model is not a simple copy of the priors, because it continues to learn human-specific structure on top of them.

## S8. Additional interpretation and cross-target evaluation of 0.10 ortholog-informed target shrinkage

This analysis changes the **supervised target**, whereas the global prior benchmark changes **input information**. The 0.10 ortholog-informed shrinkage target combines 90% target-construction human no-Gejman PC1 (minimum 10 observations) with 10% mapped target-construction reconstructed mouse PC1 (minimum 5 observations), after z-scoring each vector; genes without an ortholog signal retain the human label. We interpret the result as a human-dominant mammalian stability score derived from a multi-study compendium.

Label geometry supports human dominance. Across the 12,307 model-eligible genes, the shrinkage target correlates with target-construction human no-Gejman PC1 at 0.9982 (95% bootstrap CI, 0.9981-0.9983). In the separately defined 10,768-pair mapped set, target-construction human no-Gejman PC1 and target-construction reconstructed mouse PC1 have Pearson = 0.789, RMSE = 0.663, and MAE = 0.512 before shrinkage. After shrinkage, the target correlates with the human label at 0.9979 (RMSE = 0.065; MAE = 0.050) and with the mouse label at 0.827. Its displacement is also much smaller than ordinary study-to-study differences (Section S9 and Supplementary Table S13).

Sensitivity analysis explains the choice of $\lambda=0.10$. In matched 5-fold evaluations across three random seeds, $\lambda=0.05$, 0.10, and 0.30 gave target-versus-Saluki correlations of 0.988, 0.987, and 0.975; target-versus-mapped-mouse correlations of 0.809, 0.827, and 0.895; and prior-enhanced target-prediction correlations of 0.842, 0.854, and 0.896, respectively. The largest weight therefore produced the highest prediction score but also shifted the target more strongly toward mouse. We use 0.10 as a conservative operating point that limits cross-species contribution, not as the weight that maximizes cross-validation performance. Full results are provided in `results/ortholog_regularized_label_multiseed/summary_by_label.tsv`.

Cross-target evaluation directly tests whether changing the target harms the original human label. Using the same 12,307 genes, 1802 human-only features, 10 folds by 3 seeds, and identical model parameters, training on human no-Gejman PC1 gives r = 0.7476 ± 0.0010 against that target. Training on the shrinkage target and still evaluating against the original human target gives r = 0.7480 ± 0.0010. The seed-averaged OOF ensemble gives a paired difference of +0.0003 (95% CI, -0.0006 to 0.0012), conditional on the fitted predictions and gene resampling. The interval includes zero, so no reduction was detected; because no equivalence margin was prespecified, the result does not establish formal equivalence.

With explicit mouse priors, prediction of the shrinkage target reaches 0.8551 ± 0.0007 in the primary setting, whereas the fixed Saluki human PC1 prior-enhanced baseline on the same universe is 0.8394 ± 0.0007. The primary shrinkage-target predictor uses the audit-threshold reconstructed mouse PC1 together with released Saluki mouse PC1. To test direct reuse of the mouse vector that contributed 10% of the target, a stricter paired analysis substitutes the exact target-construction reconstructed mouse PC1. The exact dual-prior setting gives 0.8537 ± 0.0004; removing that exact vector while retaining released Saluki mouse PC1 and all availability/confidence indicators gives 0.8524 ± 0.0005. The direct-vector increment is therefore 0.0013, whereas indicators without mouse values give 0.7538 ± 0.0007. Supplementary Table S11 reports all settings. These results distinguish a small direct-target-vector contribution from the larger signal carried by external cross-species covariates.

## S9. Scale comparison between the $\lambda=0.10$ label shift and study-level noise

The main text provides the corresponding evidence in Fig. 8c and Table 4 that the label shift induced by shrinkage is much smaller than study-level noise. This section records the exact calculation protocol. Two concepts must first be separated. human no-Gejman PC1 is a consensus label integrated across many studies, whereas a study mean label is a single-study proxy obtained by averaging the samples within one study. The latter naturally carries much stronger study-specific measurement noise, so pairwise correlations between study mean labels should not be misread as the stability of the final consensus label.

The calculation proceeds in three steps. First, within the one-to-one-ortholog human gene universe, we apply sample-wise z-score to the MOESM2 human samples so that samples with different original scales are not mixed directly. Second, samples are averaged within each study to obtain a study mean label, and then each study mean label is z-scored so that it lies on the same standardized label scale as the $\lambda=0.10$ target and human no-Gejman PC1. Third, only sign alignment is applied: if the Pearson correlation between a study label and the reference label is negative, the study label is multiplied by $-1$. No linear fit, no slope recalibration, and no downstream model predictions are used. The 18 no-Gejman studies yield $\binom{18}{2}=153$ distinct study pairs. In the main analysis, the sign-alignment reference is human no-Gejman PC1. If the reference is changed directly to the $\lambda=0.10$ target, the median Pearson, RMSE, and MAE across these pairs do not change, because the $\lambda=0.10$ target is itself almost perfectly aligned with human no-Gejman PC1 (Pearson = 0.9979).

The key results are summarized in Supplementary Table S13. The $\lambda=0.10$ target is almost identical to the original human no-Gejman PC1, with Pearson = 0.9979, RMSE = 0.0646, and MAE = 0.0501. By contrast, pairwise comparisons among the normal-study mean labels show a median Pearson = 0.5303, median RMSE = 0.9540, and median MAE = 0.7203. The low median Pearson primarily reflects independent noise attenuation in two single-study proxies rather than failed sign alignment. A more direct magnitude comparison comes from the error metrics: the median RMSE and MAE of normal-study pairs are approximately 14.8-fold and 14.4-fold larger, respectively, than the $\lambda=0.10$ label shift. This scale relationship remains even if only studies with at least two samples are retained or if full Saluki-like preprocessing is performed before study averaging. The safer conclusion is therefore not that the shrinkage shift is literally negligible, but that $\lambda=0.10$ represents a shrinkage that is small relative to the scale of experimental study noise.

## S10. Prediction inputs, regulatory-feature sources, and usage boundaries

This framework cannot simply be described as “give a raw sequence and obtain the full prediction.” The current models operate at three input tiers, each with distinct applicability and performance boundaries. The first tier is base-sequence prediction. Here, the user must provide explicitly partitioned 5'UTR, CDS, and 3'UTR sequences, from which the model computes region lengths, GC content, codon frequencies, region-level 3-mers, and 3'UTR 4-mers. This tier is most suitable for fully new synthetic sequences or newly defined transcripts without external annotation, but it supports only base-sequence or human-only prediction and should not be claimed to reach the performance of the full regulatory or prior-enhanced setting.

The second tier is compact sequence/regulatory prediction, which applies to Ensembl human genes represented in the Saluki feature datapack. In this setting, the CWCS, SeqWeaver, and DeepRiPe regulatory features are not computed on the fly from a bare sequence. They come from precomputed tables: `CWCS.txt.gz`, `SeqWeaver_predictions/*_avg.txt.gz`, and `DeepRiPe_predictions_v83/*_avg.txt.gz` (Agarwal et al. 2015; Park et al. 2021; Ghanbari and Ohler 2020). For genes already represented in these tables, the regulatory blocks can be merged with base-sequence features to form `compact_all`. An Ensembl identifier alone is insufficient when the gene is absent from the datapack, and a new synthetic sequence requires rerunning the corresponding miRNA, RBP, and regulatory-prediction pipelines.

The third tier is cross-species prior-enhanced prediction for human genes with mappable mouse ortholog priors. Inputs include the two prior values (`mouse_pc1` and released Saluki mouse PC1) and their availability/high-confidence indicators. These are mouse-side covariates external to the human target. Unlike 0.10 ortholog-informed target shrinkage, this setting does **not** change the supervised target, which remains Saluki human PC1. Genes without either mouse prior should use coverage-aware fallback to the human-only compact model.

The prior-missingness stratified analysis supports this usage boundary. Across all 12,916 genes, the real both-prior model reaches Pearson = 0.8300 ± 0.0007, and the coverage-aware fallback reaches 0.8318 ± 0.0008. Among the 11,107 genes for which both mouse priors are available, the real both-prior model reaches Pearson = 0.8433 ± 0.0010, whereas the human-only model reaches 0.7456 ± 0.0014. By contrast, among the 1,347 genes for which both mouse priors are missing, the real both-prior model reaches Pearson = 0.7555 ± 0.0013, lower than the human-only value of 0.7695 ± 0.0006. The appropriate practical description is therefore that the model supports feature-based gene-level prediction; its highest performance depends on regulatory feature blocks and mouse ortholog priors; and for new sequences with only 5'UTR / CDS / 3'UTR information, the supported claim should be limited to base-sequence human-only prediction.

### Supplementary Table S6. Complete analysis-universe ledger and set relationships

The active manuscript uses 12 fixed-definition core sets arranged in three analysis branches. Eight additional sets are derived only for coverage accounting or sensitivity analysis, giving 20 named fixed rows in the machine-readable ledger. These sets are not successive exclusions from one master cohort. Gene sets are used for label construction or modeling, whereas ortholog-pair sets are used for paired cross-species comparisons. Dynamic leave-one-study-out analysis forms an additional family of per-removal intersections rather than one fixed universe: for Gejman, the primary dynamic calculation uses 15,788 genes for PC1 stability, 13,265 for Saluki agreement, and 12,592 ortholog pairs.

| ID | Track and set | N and unit | Eligibility rule | Main use |
| --- | --- | ---: | --- | --- |
| L1 | Audit-threshold human full reconstruction | 16,444 genes | MOESM2 human genes observed in at least 3 of 54 samples | Full-human PC1, sample PCA, and starting label for human LOO auditing |
| L2 | Audit-threshold mouse full reconstruction | 16,951 genes | MOESM2 mouse genes observed in at least 3 of 27 samples | Mouse PC1, mouse PCA/LOO, and fixed-target reconstructed mouse-prior source |
| L3 | Audit-threshold human no-Gejman reconstruction | 15,788 genes | Human genes observed in at least 3 of the 39 non-Gejman samples | No-Gejman PC1 and primary dynamic-coverage Gejman comparison |
| A1 | Fixed-LOO sensitivity universe | 14,244 genes | Genes passing the 3-sample rule in the full data and after every human study removal | Fixed-universe LOO sensitivity only |
| A2 | Fixed human-Saluki comparison | 13,265 genes | Intersection of human full PC1, human no-Gejman PC1, and Saluki human PC1 | Fixed-gene Saluki agreement and paired bootstrap |
| A3 | Human-mouse ortholog concordance | 12,592 ortholog pairs | Ensembl one-to-one pairs with L1, L3, and reconstructed mouse PC1 | Cross-species audit of the Gejman-related label change |
| P1 | Broad feature-target universe | 13,532 genes | Intersection of the human sequence-feature table and Saluki human PC1 | Regulatory-block ablation and MOESM3-derived-label controls |
| P2 | Global prediction universe | 12,916 genes | P1 genes also present in both L1 and L3 | Primary fixed-target human-only, transfer, permutation, and repeated-CV benchmark |
| P3 | Both-priors subset | 11,107 genes | P2 genes with both reconstructed mouse PC1 and released Saluki mouse PC1 | Prior-complete coverage analysis and two-stage residual decomposition |
| S1 | Target-construction human no-Gejman universe | 12,644 genes | Human no-Gejman genes retained at ≥10 observations; mouse comparator reconstructed at ≥5 | Construction of the 0.10 ortholog-informed target |
| S2 | Shrinkage prediction universe | 12,307 genes | S1 genes with human compact features and Saluki human PC1 | Human-only/prior target prediction, ablation, and cross-target evaluation |
| S3 | Shrinkage mapped-pair set | 10,768 ortholog pairs | S1 human genes with an Ensembl one-to-one reconstructed-mouse-PC1 match | Human-mouse label distance, target geometry, and study-noise comparison |
| D1 | High-confidence ortholog audit | 12,376 ortholog pairs | A3 pairs with Ensembl `is_high_confidence = 1` | Ortholog-confidence sensitivity |
| D2 | Fixed-target reconstructed-prior subset | 11,569 genes | P2 genes with a reconstructed mouse ortholog prior | Prior coverage and all-ortholog subset sensitivity |
| D3 | Fixed-target high-confidence subset | 11,389 genes | P2 genes with a high-confidence reconstructed ortholog | High-confidence subset sensitivity |
| D4 | Fixed-target top-regulatory-33 subset | 4,262 genes | Top third of P2 by regulatory score | Regulatory subset sensitivity in Fig. 6d |
| D5 | Only reconstructed mouse prior | 462 genes | D2 genes outside P3 | Prior-availability accounting |
| D6 | Missing both mouse priors | 1,347 genes | P2 genes outside D2 and P3 | Coverage-aware fallback and deployment boundary |
| D7 | Shrinkage prediction-mapping overlap | 10,682 genes | Human-gene intersection of S2 and S3 | Set accounting; S2 and S3 are not nested |
| D8 | High-confidence shrinkage pairs | 10,626 ortholog pairs | S3 pairs with Ensembl `is_high_confidence = 1` | Label-distance sensitivity in Supplementary Table S12 |

The fixed-target coverage groups partition P2: 11,107 genes have both priors, 462 have only the reconstructed mouse prior, none have only the released Saluki mouse prior, and 1,347 have neither. In the shrinkage branch, 10,682 genes occur in both S2 and S3. The remaining 1,625 S2 genes lack a reconstructed mouse mapping and retain the human label, whereas 86 S3 genes lack the complete prediction inputs and are excluded from S2. This 10,682-gene S2/S3 overlap was not imposed on all analyses because doing so would condition every task on mouse mapping and remove 2,234 genes (17.3%) from P2. Comparisons are instead matched within each task, using identical gene identifiers and folds.

## S11. Control experiments with MOESM3-derived supervision labels

The MOESM3 control experiments address a very specific question: because the Saluki supplement already provides a processed matrix and a "half-life (PC1)" summary label, should we replace Saluki human PC1 with labels re-derived from MOESM3 itself as the main supervision target? Supplementary Table S7 shows that the answer is no. When all human MOESM3 samples, including Gejman, are retained, the ability of the derived labels to recover Saluki human PC1 declines overall. After removing Gejman, the `as_is_pc1` label essentially returns only to a level equivalent to Saluki human PC1. Only `no_gejman_rerun_pipeline_pc1` gives a very small human-only gain, on the order of about +0.001.

The point of this result is not to create a new main label, but to explain why the main workflow does not directly switch to MOESM3-derived supervision labels. It shows that MOESM3 is better treated as a control and diagnostic resource for the Saluki processing pipeline rather than as a better primary source of task labels.

### Supplementary Table S7. Downstream benchmark for MOESM3-derived supervision labels

| Feature set | Supervision label | Genes | Label-Saluki r | OOF r, training target | OOF r, Saluki target | OOF R², Saluki target |
|:--|:--|--:|--:|--:|--:|--:|
| Compact human | MOESM3 no-Gejman: pipeline rerun | 13,532 | 1.000 | 0.730 | 0.730 | 0.533 |
| Compact human | Saluki human PC1 | 13,532 | 1.000 | 0.729 | 0.729 | 0.530 |
| Compact human | MOESM3 no-Gejman: as-is PC1 | 13,532 | 1.000 | 0.729 | 0.729 | 0.530 |
| Compact human | MOESM3 no-Gejman: z-score only | 13,532 | 1.000 | 0.729 | 0.729 | 0.530 |
| Compact human | MOESM3 all samples: z-score only | 13,532 | 0.963 | 0.707 | 0.722 | 0.520 |
| Compact human | MOESM3 all samples: pipeline rerun | 13,532 | 0.963 | 0.707 | 0.721 | 0.519 |
| Compact human | MOESM3 all samples: as-is PC1 | 13,532 | 0.963 | 0.707 | 0.721 | 0.520 |
| Base sequence | MOESM3 no-Gejman: z-score only | 13,532 | 1.000 | 0.720 | 0.720 | 0.518 |
| Base sequence | Saluki human PC1 | 13,532 | 1.000 | 0.720 | 0.720 | 0.517 |
| Base sequence | MOESM3 no-Gejman: as-is PC1 | 13,532 | 1.000 | 0.720 | 0.720 | 0.517 |
| Base sequence | MOESM3 no-Gejman: pipeline rerun | 13,532 | 1.000 | 0.718 | 0.718 | 0.516 |
| Base sequence | MOESM3 all samples: pipeline rerun | 13,532 | 0.963 | 0.698 | 0.713 | 0.507 |
| Base sequence | MOESM3 all samples: z-score only | 13,532 | 0.963 | 0.697 | 0.713 | 0.506 |
| Base sequence | MOESM3 all samples: as-is PC1 | 13,532 | 0.963 | 0.696 | 0.712 | 0.507 |

Reader-facing labels are used in the table. Exact setting identifiers and full-precision values are retained in `results/moesm3_saluki_label_benchmark/summary.tsv`.

## S12. Focused control of MOESM3-derived labels under mouse-prior settings

A stricter question is whether these MOESM3-derived labels suddenly become better than directly training on Saluki human PC1 once mouse priors are allowed. Supplementary Table S8 again shows that the answer is no. In the focused prior benchmark, Saluki human PC1 (code field `Saluki_human_PC1`) combined with `compact_all_plus_both_mouse_priors_global` reaches a Pearson correlation of 0.815917, whereas the best alternative label, `moesm3_no_gejman_rerun_pipeline_pc1 + compact_all_plus_both_mouse_priors_global`, reaches 0.815588, which is still slightly lower.

Taken together, these two control experiments support a more stable submission-ready statement: MOESM3-derived labels can serve as interpretive controls for the Saluki labeling process, but even in the cross-species transfer setting they are not sufficient to replace Saluki human PC1 as the primary supervision target.

### Supplementary Table S8. Focused benchmark of MOESM3-derived supervision labels under prior settings

| Supervision label | Feature setting | Genes | Reconstructed-prior genes | Saluki-prior genes | Label-Saluki r | OOF r, training target | OOF r, Saluki target | OOF R², Saluki target |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|
| Saluki human PC1 | Compact human plus both mouse priors | 13,532 | 11,569 | 11,423 | 1.000 | 0.816 | 0.816 | 0.666 |
| MOESM3 no-Gejman: pipeline rerun | Compact human plus both mouse priors | 13,532 | 11,569 | 11,423 | 1.000 | 0.816 | 0.816 | 0.665 |
| MOESM3 no-Gejman: pipeline rerun | Compact human plus reconstructed mouse prior | 13,532 | 11,569 | 11,423 | 1.000 | 0.814 | 0.814 | 0.662 |
| Saluki human PC1 | Compact human plus reconstructed mouse prior | 13,532 | 11,569 | 11,423 | 1.000 | 0.814 | 0.814 | 0.662 |
| Saluki human PC1 | Compact human plus Saluki mouse prior | 13,532 | 11,569 | 11,423 | 1.000 | 0.812 | 0.812 | 0.659 |
| MOESM3 no-Gejman: pipeline rerun | Compact human plus Saluki mouse prior | 13,532 | 11,569 | 11,423 | 1.000 | 0.812 | 0.811 | 0.658 |
| MOESM3 no-Gejman: pipeline rerun | Compact human only | 13,532 | 11,569 | 11,423 | 1.000 | 0.730 | 0.730 | 0.533 |
| Saluki human PC1 | Compact human only | 13,532 | 11,569 | 11,423 | 1.000 | 0.729 | 0.729 | 0.530 |

Reader-facing labels are used in the table. Exact setting identifiers and full-precision values are retained in `results/moesm3_saluki_label_prior_benchmark_smoketest/summary.tsv`.

## S13. Supplementary tables for 10-fold robustness, weak ortholog-informed target shrinkage, and uncertainty

Supplementary Tables S9-S15 collect the strict checks. S9 reports 10-fold robustness for the global prior benchmark. S10 compares the 0.10 ortholog-informed shrinkage target with the fixed Saluki-PC1 baseline on the same universe. S11 reports prior ablations. S12-S13 report label distance and study-noise scale. S14 reports human-only cross-target evaluation. S15 provides paired-bootstrap CIs for label auditing, ortholog concordance, and prior gain.

### Supplementary Table S9. 10-fold robustness of the global prior benchmark

| Setting | Prior mode | Folds | Seeds | Genes | Pearson r, mean ± SD | Spearman ρ, mean ± SD | R², mean ± SD |
|:--|:--|--:|--:|--:|:--|:--|:--|
| Compact human plus both mouse priors | Real priors | 10 | 3 | 12,916 | 0.830 ± 0.001 | 0.824 ± 0.001 | 0.689 ± 0.001 |
| Compact human plus both mouse priors | Shuffled priors | 10 | 3 | 12,916 | 0.748 ± 0.001 | 0.741 ± 0.001 | 0.558 ± 0.001 |
| Compact human only | None | 10 | 3 | 12,916 | 0.748 ± 0.001 | 0.740 ± 0.001 | 0.558 ± 0.002 |
| Coverage-aware fallback | Availability-based selection | 10 | 3 | 12,916 | 0.832 ± 0.001 | 0.826 ± 0.001 | 0.690 ± 0.001 |

### Supplementary Table S10. 10-fold target benchmark across three random seeds

| Target and mouse comparator | λ | Genes | Target-Saluki r | Target-mouse r | Real-prior OOF r | Human-only OOF r | Shuffled-prior OOF r | Real-shuffled Δr |
|:--|--:|--:|--:|--:|:--|:--|:--|:--|
| 0.10 shrinkage; target-construction mouse PC1 | 0.10 | 12,307 | 0.987 | 0.827 | 0.855 ± 0.001 | 0.754 ± 0.001 | 0.753 ± 0.001 | 0.102 ± 0.001 |
| Saluki human PC1; released Saluki mouse PC1 | 0 | 12,307 | 1.000 | 0.798 | 0.839 ± 0.001 | 0.754 ± 0.001 | 0.754 ± 0.001 | 0.085 ± 0.000 |

### Supplementary Table S11. 10-fold prior-ablation benchmark for 0.10 ortholog-informed target shrinkage

| Setting | Features | Folds | Seeds | OOF r, target | OOF r, Saluki | Δr vs fixed-target prior baseline |
|:--|--:|--:|--:|:--|:--|:--|
| Human-only features | 1802 | 10 | 3 | 0.7538 ± 0.0010 | 0.7533 ± 0.0013 | -0.0857 ± 0.0010 |
| Availability/confidence indicators only | 1806 | 10 | 3 | 0.7538 ± 0.0007 | 0.7533 ± 0.0009 | -0.0857 ± 0.0007 |
| Released Saluki mouse PC1 only | 1805 | 10 | 3 | 0.8521 ± 0.0005 | 0.8347 ± 0.0006 | 0.0127 ± 0.0005 |
| Strict no-direct: Saluki value plus all indicators | 1807 | 10 | 3 | 0.8524 ± 0.0005 | 0.8351 ± 0.0005 | 0.0130 ± 0.0005 |
| Exact target-construction mouse PC1 only | 1805 | 10 | 3 | 0.8525 ± 0.0003 | 0.8341 ± 0.0004 | 0.0130 ± 0.0003 |
| Exact target-construction plus Saluki mouse PC1 | 1808 | 10 | 3 | 0.8537 ± 0.0004 | 0.8354 ± 0.0004 | 0.0143 ± 0.0004 |
| Primary audit-threshold plus Saluki mouse PC1 | 1808 | 10 | 3 | 0.8551 ± 0.0007 | 0.8366 ± 0.0006 | 0.0156 ± 0.0007 |

The primary row reproduces the manuscript's original target-prediction input. The strict direct-vector comparison is the sixth row minus the fourth row: adding the exact target-construction mouse PC1 changes target correlation by 0.0013 while released Saluki mouse PC1 and all indicators remain fixed.

### Supplementary Table S12. Human-mouse ortholog label distance

| Comparison | Ortholog subset | N | Pearson r | Spearman ρ | MSE | RMSE | MAE | Median absolute error | 95th-percentile absolute error |
|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| Human no-Gejman PC1 vs target-construction mouse PC1 | All one-to-one | 10,768 | 0.789 | 0.782 | 0.440 | 0.663 | 0.512 | 0.413 | 1.330 |
| Shrinkage target vs target-construction mouse PC1 | All one-to-one | 10,768 | 0.827 | 0.820 | 0.361 | 0.601 | 0.464 | 0.373 | 1.209 |
| Shrinkage target vs human no-Gejman PC1 | All one-to-one | 10,768 | 0.998 | 0.998 | 0.004 | 0.065 | 0.050 | 0.041 | 0.129 |
| Saluki human PC1 vs released Saluki mouse PC1 | All one-to-one | 10,768 | 0.798 | 0.788 | 0.418 | 0.646 | 0.503 | 0.411 | 1.285 |
| Shrinkage target vs Saluki human PC1 | All one-to-one | 10,768 | 0.990 | 0.992 | 0.020 | 0.142 | 0.094 | 0.063 | 0.291 |
| Human no-Gejman PC1 vs target-construction mouse PC1 | High confidence | 10,626 | 0.791 | 0.784 | 0.436 | 0.660 | 0.510 | 0.410 | 1.326 |
| Shrinkage target vs target-construction mouse PC1 | High confidence | 10,626 | 0.829 | 0.822 | 0.357 | 0.598 | 0.462 | 0.371 | 1.201 |
| Shrinkage target vs human no-Gejman PC1 | High confidence | 10,626 | 0.998 | 0.998 | 0.004 | 0.064 | 0.050 | 0.041 | 0.129 |
| Saluki human PC1 vs released Saluki mouse PC1 | High confidence | 10,626 | 0.800 | 0.789 | 0.414 | 0.644 | 0.500 | 0.409 | 1.275 |
| Shrinkage target vs Saluki human PC1 | High confidence | 10,626 | 0.990 | 0.992 | 0.020 | 0.142 | 0.094 | 0.063 | 0.290 |

### Supplementary Table S13. Study-level noise versus $\lambda=0.10$ label shift

| Comparison | Cohort | Comparisons | Median Pearson r | Median residual r | Median RMSE | Median MAE |
|:--|:--|--:|--:|--:|--:|--:|
| $\lambda = 0.10$ target vs human no-Gejman PC1 | all one-to-one ortholog genes | 1 | 0.9979 | NA | 0.0646 | 0.0501 |
| study mean label pairs | no-Gejman studies | 153 | 0.5303 | NA | 0.9540 | 0.7203 |
| study-specific residual pairs | no-Gejman studies | 153 | NA | 0.0796 | 0.9540 | 0.7203 |
| study labels vs human no-Gejman PC1 | no-Gejman studies | 18 | 0.7420 | NA | 0.7145 | 0.5429 |
| study mean label pairs | no-Gejman studies with at least two samples | 21 | 0.5756 | NA | 0.8963 | 0.6963 |
| study-specific residual pairs | no-Gejman studies with at least two samples | 21 | NA | 0.0337 | 0.8963 | 0.6963 |
| study mean label pairs | no-Gejman full preprocessing before study averaging | 153 | 0.5937 | NA | 0.9014 | 0.6779 |
| study-specific residual pairs | no-Gejman full preprocessing before study averaging | 153 | NA | 0.1026 | 0.9014 | 0.6779 |

### Supplementary Table S14. Human-only cross-target evaluation of the 0.10 ortholog-informed shrinkage target

Both training targets use the same 12,307 genes, 1802 human-only `compact_all` features, 10 folds, 3 random seeds, and identical hyperparameters. SD denotes across-seed split sensitivity.

| Training target | Evaluation target | Genes | Pearson r, mean ± SD | R², mean ± SD | RMSE, mean ± SD | MAE, mean ± SD |
|:--|:--|--:|:--|:--|:--|:--|
| human no-Gejman PC1 | human no-Gejman PC1 | 12,307 | 0.7476 ± 0.0010 | 0.5572 ± 0.0015 | 0.6687 ± 0.0012 | 0.5260 ± 0.0013 |
| 0.10 ortholog-informed shrinkage target | human no-Gejman PC1 | 12,307 | 0.7480 ± 0.0010 | 0.5582 ± 0.0014 | 0.6680 ± 0.0011 | 0.5248 ± 0.0009 |
| human no-Gejman PC1 | 0.10 ortholog-informed shrinkage target | 12,307 | 0.7532 ± 0.0010 | 0.5651 ± 0.0015 | 0.6627 ± 0.0012 | 0.5215 ± 0.0014 |
| 0.10 ortholog-informed shrinkage target | 0.10 ortholog-informed shrinkage target | 12,307 | 0.7538 ± 0.0010 | 0.5664 ± 0.0014 | 0.6617 ± 0.0010 | 0.5201 ± 0.0009 |
| human no-Gejman PC1 | Saluki human PC1 | 12,307 | 0.7529 ± 0.0014 | 0.5638 ± 0.0021 | 0.6723 ± 0.0016 | 0.5273 ± 0.0017 |
| 0.10 ortholog-informed shrinkage target | Saluki human PC1 | 12,307 | 0.7533 ± 0.0013 | 0.5650 ± 0.0019 | 0.6714 ± 0.0015 | 0.5261 ± 0.0014 |

Against the original human no-Gejman PC1, the seed-averaged OOF ensemble gave Pearson = 0.75233 after human-target training and 0.75259 after shrinkage-target training. The paired difference was +0.00026, with a 95% CI of -0.00061 to 0.00118 from 2,000 gene bootstraps. This interval is conditional on the fitted three-seed ensemble and gene resampling; it does not represent retraining or split uncertainty. Because the interval crosses zero, the analysis provides no evidence that shrinkage-target training either reduces or improves predictability of the original human label at the present resolution; it is not presented as a formal equivalence test with a prespecified margin.

### Supplementary Table S15. Paired-bootstrap uncertainty for key correlation differences

The cross-target analysis above uses 2,000 bootstraps. This table uses 5,000 gene or one-to-one ortholog-pair bootstraps. Prior predictions are seed-averaged OOF ensembles formed as gene-wise means of the three 10-fold OOF vectors. Their intervals condition on the fitted predictions and resampled genes or pairs; they do not include retraining uncertainty. Resampling treats genes or ortholog pairs as units and does not model dependence within pathways or paralog families.

All prediction folds in this study are gene-level random splits rather than paralog-family-grouped splits. The reported OOF values therefore estimate generalization to held-out genes under this split rule; gene-family-held-out generalization remains untested.

| claim | N | candidate A Pearson | candidate B Pearson | A - B | 95% CI for difference |
|:--|--:|--:|--:|--:|:--|
| Saluki agreement: no-Gejman vs full human PC1 | 13,265 | 0.98294 | 0.95152 | 0.03142 | 0.02978 to 0.03303 |
| Ortholog concordance: no-Gejman vs full, all one-to-one | 12,592 | 0.75528 | 0.74062 | 0.01466 | 0.01110 to 0.01814 |
| Ortholog concordance: no-Gejman vs full, high confidence | 12,376 | 0.75829 | 0.74366 | 0.01463 | 0.01099 to 0.01811 |
| Real prior vs human-only seed-averaged OOF prediction | 12,916 | 0.83244 | 0.75253 | 0.07991 | 0.07398 to 0.08617 |
| Real prior vs shuffled-prior seed-averaged OOF prediction | 12,916 | 0.83244 | 0.75306 | 0.07938 | 0.07348 to 0.08554 |

## S14. Claim boundaries

Both the main text and the Supplementary Information use four interpretation rules. First, the pure human route is a transparent human-only fixed-target benchmark. Second, higher fixed-target performance is limited to the cross-species transfer setting because the inputs explicitly contain mouse ortholog priors. Third, 0.10 ortholog-informed target shrinkage is a core label-construction result, but 0.855 describes predictability of the shrinkage target from explicit cross-species inputs. Cross-target evaluation detected no reduction in human-label predictability; because no equivalence margin was prespecified and the confidence interval includes zero, it is not a formal equivalence result. Fourth, complete compact/prior-enhanced prediction depends on precomputed regulatory features and mouse ortholog priors, whereas a new sequence can enter only the base-sequence tier after 5'UTR / CDS / 3'UTR annotation.

## S15. Reproducibility entry points

Tables 1 and 2 in the main text summarize the method modules and analysis universes without listing local absolute paths. Concrete result files, script entry points, and reproduction commands are documented in this section and in the public versioned repository: `https://github.com/wwzdl/mrna-pc1-label`. All paths below are repository-relative so that reviewers and readers can reproduce the analysis across environments.

From the repository root, the staged reproduction entry points are:

```bash
bash scripts/create_env.sh
source scripts/activate_env.sh
bash scripts/reproduce_bmb_key_results.sh labels
bash scripts/reproduce_bmb_key_results.sh models
bash scripts/reproduce_bmb_key_results.sh figures
```

The `labels` stage downloads MOESM2/MOESM3 and runs label reconstruction, study influence, ortholog concordance, label-distance, study-noise, and PCA-imputation analyses. The `models` stage requires the 18.8-GB public Saluki datapack (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)) staged as documented in the repository README and a CUDA-capable XGBoost environment; if the human sequence-feature table is absent, this stage rebuilds it from Ensembl release 115. The `figures` stage redraws all active figures from fixed result tables, rebuilds the DOCX/PDF files, and runs the submission audit. MOESM2, MOESM3, and the Ensembl release 115 mouse-homology file are checksum-verified by the download/analysis code. The repository includes a portable dependency file, a Conda specification, and `requirements-validated.txt`, which records the direct package versions that passed the final pre-submission audit.

### Supplementary Table S16. Main reproduction entry points and result files

| Category | Path or entry point | Purpose |
|:--|:--|:--|
| Public repository | `https://github.com/wwzdl/mrna-pc1-label` | Provides code, result tables, figure scripts, and environment notes; audited two-author manuscript tag: `mRNA-PC1-label-v1.4.3` |
| Staged reproduction driver | `scripts/reproduce_bmb_key_results.sh` | Runs the `labels`, `models`, and `figures` stages with the active result paths |
| Environment and input integrity | `requirements.txt`, `environment.yml`, `requirements-validated.txt`, `scripts/fetch_real_data.sh` | Records portable and validated dependencies and verifies the public supplementary inputs by checksum |
| OOF integrity audit | `scripts/audit_oof_integrity.py` | Checks canonical gene counts, unique IDs, shared universes across settings/seeds, fold coverage, and complete predictions |
| Analysis-universe ledger | `scripts/audit_analysis_universes.py`, `results/bmb_analysis_universe_ledger.tsv` | Recomputes the twelve core sets and derived overlaps from source/result tables and verifies their set relationships |
| Main and supplementary manuscripts | `manuscript/bmb_submission/` | Stores the Chinese/English main text, supplementary material, DOCX drafts, and submission figures |
| Sample-level PCA imputation reproducibility | `scripts/run_sample_pca_imputation_reproducibility.sh`, `scripts/reproduce_sample_pca_imputation.py` | Reproduces low-rank imputation, coordinates, parameter JSON, and diagnostic figures for raw/processed sample PCA |
| Human-only benchmark | `results/human_compact_oof_models.tsv` | Summarizes out-of-fold results for base-sequence and compact-feature models |
| XGBoost parameter provenance | `scripts/reproduce_xgb_parameter_scan.py`, `results/human_xgboost_gpu_saluki_refine.tsv`, `results/human_xgboost_parameter_provenance.json` | Reproduces the preliminary same-compendium five-setting scan and records the frozen selection rule and interpretation limit |
| Global prior benchmark | `results/human_global_mouse_prior_benchmark.tsv`, `results/tenfold_global_prior/summary_by_setting.tsv` | Records mouse-prior-enhanced models and 10-fold robustness |
| Prior control | `results/human_global_mouse_prior_permutation_control.tsv`, `results/prior_residual_analysis/summary.tsv` | Supports the shuffled-prior control and the two-stage residual analysis |
| 0.10 ortholog-informed target shrinkage | `results/ortholog_regularized_label_10fold_main/summary_by_label.tsv`, `results/ortholog_regularized_label_multiseed/summary_by_label.tsv`, `results/ortholog_prior_ablation_lambda0p1_10fold/summary_by_setting.tsv` | Records the $\lambda = 0.10$ target benchmark, λ sensitivity, and prior ablations |
| Human-only cross-target evaluation | `results/ortholog_regularized_cross_target/summary_aggregate.tsv`, `results/ortholog_regularized_cross_target/paired_bootstrap_ci.tsv` | Tests for a change in original-human-label predictability after shrinkage-target training |
| Key-claim uncertainty | `results/key_claim_bootstrap/paired_bootstrap_summary.tsv` | Provides paired-bootstrap CIs for Gejman, ortholog concordance, and real-prior gain |
| Label distance and study noise | `results/ortholog_label_distance/distance_summary.tsv`, `results/study_noise_vs_orthoreg_shift/noise_summary.tsv` | Compares label shift, human-mouse distance, and study-to-study variability |
| Study-influence sensitivity | `src/mrna_half_life_paper/study_influence_sensitivity.py`, `results/study_influence_sensitivity/` | Reproduces conditional same-size deletion, fixed-universe, preprocessing, and study-balanced controls |
| Main-text weak-shrinkage validation figure | `scripts/plot_bmb_ortholog_regularized_target.py`, `manuscript/bmb_submission/figures/main/Fig08_ortholog_regularized_target_cn.*` | Generates the label-geometry, study-noise, and cross-target panels in Fig. 8 |
| Main-text ortholog concordance figure | `scripts/redraw_bmb_ortholog_combined_cn.py`, `manuscript/bmb_submission/figures/main/Fig05_ortholog_summary_cn.*` | Generates the paired scatter, correlation-summary, and gain panels in Fig. 5 directly from the ortholog result tables |
| Label auditing and ortholog concordance | `results/study_influence/**/*.tsv`, `results/ortholog_analysis/*.tsv` | Stores leave-one-study-out, Saluki-label agreement, and human-mouse consistency results |
| Usage-boundary analysis | `results/prior_missingness_analysis/summary_by_coverage.tsv`, `docs/prior_residual_analysis.md` | Documents the performance boundary under missing priors and the residual-analysis implementation |

The main figure and manuscript-generation scripts include `scripts/redraw_bmb_dataset_summary_cn.py`, `scripts/redraw_bmb_human_pca_panel.py`, `scripts/redraw_bmb_study_influence_cn.py`, `scripts/redraw_bmb_ortholog_combined_cn.py`, `scripts/plot_bmb_ortholog_regularized_target.py`, `scripts/redraw_bmb_sequence_progression_cn.py`, `scripts/redraw_bmb_prior_summary_cn.py`, `scripts/redraw_bmb_flow_figures.py`, `scripts/redraw_bmb_supplementary_diagnostics.py`, `scripts/render_markdown_to_docx.py`, and `scripts/render_markdown_to_pdf.py`. `scripts/audit_analysis_universes.py` rebuilds the machine-readable universe ledger and checks the core/derived set relationships used in the text. The flow and diagnostic scripts use native matplotlib vector objects and do not crop or embed AI-generated artwork. Figures are exported as editable SVG/PDF and 600-dpi PNG/TIFF. `ortholog_regularized_cross_target.py` and `bootstrap_bmb_key_claims.py` generate the cross-target and paired-bootstrap results, respectively. Markdown tables are retained as the version-controlled source and rendered in compact DOCX/PDF submission formats.

## S16. References

Agarwal V, Bell GW, Nam JW, Bartel DP (2015) Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://doi.org/10.7554/eLife.05005

Agarwal V, Kelley DR (2022) The genetic and biochemical determinants of mRNA degradation rates in mammals. Genome Biology 23(1):245. https://doi.org/10.1186/s13059-022-02811-x

Bolstad BM, Irizarry RA, Astrand M, Speed TP (2003) A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics 19(2):185-193. https://doi.org/10.1093/bioinformatics/19.2.185

Chen T, Guestrin C (2016) XGBoost: A scalable tree boosting system. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining 785-794. https://doi.org/10.1145/2939672.2939785

Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al (2025) Ensembl 2025. Nucleic Acids Research 53(D1):D948-D957. https://doi.org/10.1093/nar/gkae1071

Gabaldon T, Koonin EV (2013) Functional and evolutionary implications of gene orthology. Nature Reviews Genetics 14(5):360-366. https://doi.org/10.1038/nrg3456

Ghanbari M, Ohler U (2020) Deep neural networks for interpreting RNA-binding protein target preferences. Genome Research 30(2):214-226. https://doi.org/10.1101/gr.247494.118

Park CY, Zhou J, Wong AK, Chen KM, Theesfeld CL, Darnell RB, et al (2021) Genome-wide landscape of RNA-binding protein target site dysregulation reveals a major impact on psychiatric disorder risk. Nature Genetics 53:166-173. https://doi.org/10.1038/s41588-020-00761-3

Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, Blondel M, Prettenhofer P, Weiss R, Dubourg V, et al (2011) Scikit-learn: Machine learning in Python. Journal of Machine Learning Research 12:2825-2830. https://jmlr.org/papers/v12/pedregosa11a.html

Tipping ME, Bishop CM (1999) Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B 61(3):611-622. https://doi.org/10.1111/1467-9868.00196

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R, Botstein D, Altman RB (2001) Missing value estimation methods for DNA microarrays. Bioinformatics 17(6):520-525. https://doi.org/10.1093/bioinformatics/17.6.520
