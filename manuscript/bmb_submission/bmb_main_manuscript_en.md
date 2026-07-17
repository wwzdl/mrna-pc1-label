# Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction

**Article type:** Original research article  
**Authors:** Wenzhuo Wang<sup>1</sup>, Ying Shao*<sup>1</sup>  
**Affiliation:** <sup>1</sup> School of Science, Dalian Maritime University, Dalian 116026, China  
**Corresponding author:** Ying Shao; Email: yshao@dlmu.edu.cn<br>
**ORCID:** Ying Shao: https://orcid.org/0000-0002-4056-5757<br>
**Funding:** No specific funding was received for this work.  
**Keywords:** mRNA half-life; study heterogeneity; label audit; target shrinkage; cross-species transfer; XGBoost

## Abstract

Mammalian mRNA half-life compendia combine heterogeneous protocols, cellular states, and processing pipelines, confounding label quality with model performance. We present a reproducible framework for consensus-label auditing, fixed-target prediction, cross-species transfer, and weak ortholog-informed target shrinkage. On 13,265 common genes, removing `Gejman` increased agreement with the Saluki human first principal component (PC1) label from r = 0.9515 to r = 0.9829. Its PC1 displacement fell within a null from 500 size-matched random removals, and its geometric rank fell from first to third under two study-balanced estimators. However, no random removal matched its increase in Saluki agreement or human-mouse ortholog concordance (one-sided empirical p = 0.002 for each); ortholog correlation increased from r = 0.741 to r = 0.755. With Saluki human PC1 fixed, the human-only model achieved r = 0.748 ± 0.001 in repeated 10-fold evaluation. External mouse ortholog priors raised performance to r = 0.830 ± 0.001; gene-prior permutation returned it to r = 0.748 ± 0.001. Across the 12,307-gene target universe, a target formed by 0.10 ortholog-informed shrinkage remained nearly identical to human no-Gejman PC1 (r = 0.998); within the 10,768 mapped one-to-one pairs, root mean squared error was 0.065. Cross-target evaluation did not detect a reduction in human-label predictability (difference = +0.0003; 95% confidence interval, -0.0006 to 0.0012). These analyses distinguish sample-size effects, cross-species covariate gains, and changes to the supervised target.

## 1. Introduction

mRNA stability regulation is an indispensable layer of gene-expression control. Even at the same transcriptional output, different transcripts can exhibit markedly different steady-state abundance and translational output because of differences in degradation rate (Ross 1995; Garneau et al. 2007; Schoenberg and Maquat 2012; Rambout and Maquat 2024). Classical and recent studies have shown that codon composition, ribosome elongation, 3' untranslated region (3'UTR) cis-elements, RNA-binding protein occupancy, RNA structure, and RNA modification all influence mRNA half-life (Meyer et al. 2012; Wang et al. 2014; Presnyak et al. 2015; Wu et al. 2019; Hia et al. 2019; Grimson et al. 2007; Agarwal et al. 2015; Mauger et al. 2019). Half-life prediction therefore matters for basic biology, variant interpretation, and RNA therapeutic engineering (Leppek et al. 2022; Sample et al. 2019; Cetnar et al. 2024; Musaev et al. 2024).

At the data level, mammalian mRNA half-life measurements have evolved from small-scale experiments to cross-study compendia. Transcriptional shutoff assays, 4sU/5EU metabolic labeling, nucleoside-recoding approaches, and related high-throughput workflows now provide degradation estimates for thousands to tens of thousands of genes in diverse cellular contexts (Schwanhausser et al. 2011; Rabani et al. 2011; Tani et al. 2012; Schwalb et al. 2012; Duffy et al. 2015; Herzog et al. 2017; Schofield et al. 2018; Muhar et al. 2018). Heterogeneity is the main challenge in using these resources. Batch effects and study-specific processing can remain prominent in high-throughput matrices even after routine normalization (Johnson et al. 2007; Leek et al. 2010). Study source, experimental protocol, time-point design, cellular context, and downstream transformation steps all leave structure in the gene by sample matrix. Treating these observations as homogeneous samples may cause a model to absorb technical structure together with biological structure.

Saluki and related work provided an important first-generation solution to this problem by integrating multiple studies into a compendium and extracting a consensus first principal component (PC1) label after unified preprocessing (Agarwal and Kelley 2022). That work provided both a mammalian half-life benchmark and a strong demonstration of sequence-based prediction. At the same time, it exposed a deeper issue. If the consensus label itself still carries study-level bias, downstream performance comparisons reflect both model capacity and label quality. A benchmark-oriented computational biology study must address this issue explicitly.

We place label auditing, rather than model scale, at the center of the paper. Complete leave-one-study-out reconstruction screens influential studies by PC1 stability without using Saluki labels. Saluki human PC1 agreement and human-mouse ortholog concordance then assess the direction of each change. With Saluki human PC1 fixed, we compare a human-only baseline with cross-species transfer using external mouse ortholog priors. We also examine a target formed by weak ortholog-informed shrinkage at $\lambda=0.10$. Input augmentation and target shrinkage are evaluated separately because they change different parts of the prediction problem.

## 2. Related Work

### 2.1 Large-scale measurement and integration of mammalian mRNA half-life

Early reviews and systematic studies established that variation in transcript degradation rate is a major source of gene-expression regulation in mammals (Ross 1995; Garneau et al. 2007; Schoenberg and Maquat 2012). Transcriptome-scale experiments subsequently accumulated across classical shutoff assays, metabolic labeling, and nucleoside-recoding strategies (Schwanhausser et al. 2011; Rabani et al. 2011; Tani et al. 2012; Schwalb et al. 2012; Duffy et al. 2015; Herzog et al. 2017; Schofield et al. 2018; Muhar et al. 2018). Duan et al. quantified interindividual RNA-stability variation in human lymphoblastoid cell lines; this source is represented as the `Gejman` study in the Saluki compendium (Duan et al. 2013). Friedel et al. highlighted conserved principles of mammalian half-life regulation (Friedel et al. 2009), and Agarwal and Kelley advanced the field to a systematic compendium and sequence-based benchmark (Agarwal and Kelley 2022). However, most existing studies treat the construction of a gene-level consensus label from a noisy multi-study matrix as an implicit preprocessing step rather than as an independent methodological problem.

### 2.2 Sequence and biochemical determinants of mRNA stability

Mechanistically, codon optimality, translation efficiency, and degradation coupling have become some of the most influential findings in the field (Presnyak et al. 2015; Wu et al. 2019; Hia et al. 2019). In parallel, microRNA (miRNA) targeting rules, 3'UTR cis-elements, and RNA structure have all been shown to affect half-life systematically (Mayr 2017; Grimson et al. 2007; Agarwal et al. 2015; Mauger et al. 2019). RNA modifications such as N6-methyladenosine (m6A) provide an additional explanatory layer (Meyer et al. 2012; Wang et al. 2014). Together, these findings give sequence-to-stability prediction clear biological grounding.

### 2.3 Computational modeling from tabular learning to deep sequence models

At the modeling level, random forests, gradient boosting, and their modern implementations such as XGBoost and LightGBM remain strong baselines for heterogeneous feature sets (Breiman 2001; Friedman 2001; Chen and Guestrin 2016; Ke et al. 2017). Models closer to raw sequence input include DeepBind, cross-species regulatory sequence activity predictors, Enformer, and DNABERT (Alipanahi et al. 2015; Kelley 2020; Avsec et al. 2021; Ji et al. 2021). Here we use a unified and transparent feature-learning framework to quantify how label diagnosis and cross-species priors affect predictive performance.

## 3. Materials and Methods

Fig. 1 summarizes the analysis, and Table 1 defines the inputs, operations, and validation roles. Label auditing precedes prediction. Fixed-target analyses compare human-only inputs with external mouse ortholog covariates, whereas target-shrinkage analyses modify the supervised target. The two prediction settings are evaluated separately. Reproduction paths and commands are provided in Online Resource 1, Supplementary Section S15, and Data and Code Availability.

![](figures/main/Fig01_workflow.png)

**Fig. 1** Overall study-aware framework. a, Consensus-label reconstruction and auditing: `MOESM2` is processed using a minimum of three observed samples per gene, sample-wise z-scoring, iterative PCA imputation, quantile normalization, and PC1 extraction before leave-one-study-out analysis; `MOESM3` and Ensembl orthologs provide the Saluki processing reference and cross-species comparison axis. The directional gains are shown together with the nonsignificant size-matched geometry-null result. b, Fixed-target prediction: Saluki human PC1 remains the target while human-only and mouse-prior transfer inputs are evaluated, with fold-wise permutation testing gene-prior identity. c, Ortholog-informed target shrinkage: 90% human no-Gejman PC1 and 10% mapped mouse PC1 define the target, which is checked by label distance, study-noise comparison, and cross-target evaluation. Values are generated from result tables by a reproducible vector-graphics script

**Table 1** Method modules, logical inputs, and validation roles. The main text retains only the analytical logic and intended role of each module. Specific filenames, script paths, result tables, and reproduction commands are listed in Supplementary Section S15 and in the versioned GitHub repository specified in Data and Code Availability.

| Method module | Logical input | Core operation | Validation role |
| --- | --- | --- | --- |
| Data organization | Public human/mouse half-life matrices and sample metadata | Parse species, study, method, cell type, and replicate; define human, mouse, and no-Gejman subsets | Clarify sample structure and gene universes, and prevent different task settings from being conflated |
| Consensus-label reconstruction | Gene by sample half-life matrix | Coverage filtering, sample standardization, missing-value imputation, distribution alignment, and PC1 extraction | Generate a reproducible half-life consensus label and align it to the released Saluki PC1 label |
| Study influence | Study-specific sample subsets | Rebuild the label after removing one study at a time; screen by PC1 stability and examine direction by Saluki agreement | Separate Saluki-label-independent screening from the published-processing comparison; downstream scores are not used |
| Ortholog concordance | Human-mouse one-to-one ortholog genes | Compare human full PC1, no-Gejman PC1, Saluki human PC1, and mouse labels | Test whether the label change moves toward cross-species conserved structure |
| Human-only prediction | Human sequence/regulatory features | Gene-level out-of-fold (OOF) XGBoost benchmark | Quantify performance from human-only representation and define a transparent human-only baseline |
| Prior-enhanced prediction | Human features plus mouse ortholog priors | Prior-enhanced OOF benchmark, permutation control, and residual analysis | Test gene-specific cross-species signal in a transfer setting whose mouse-side covariates are external to the human target |
| 0.10 ortholog-informed target shrinkage | Human consensus label plus mapped mouse ortholog label | Mix at 90:10; quantify label shift, study-noise scale, and cross-target performance | Define a human-dominant mammalian stability target, interpreted separately from the fixed-target benchmark |

### 3.1 Data sources, metadata parsing, and common gene universes

We used the public supplementary data from the Saluki study by Agarwal and Kelley (Agarwal and Kelley 2022). `MOESM2` provides transformed human and mouse mRNA half-life gene by sample matrices and serves as the primary input for label reconstruction and leave-one-study-out auditing. Its first two columns are Ensembl gene ID and gene name; the remaining columns contain transformed values from different studies, methods, cell types, and replicates. `MOESM3` provides the released `half-life (PC1)` summary label together with processed sample values. Throughout this manuscript, **Saluki human PC1** denotes the released human PC1 label, whereas **released Saluki mouse PC1** denotes the released mouse PC1 label, which is mapped through orthologs and used here as an external mouse-side prior. The human worksheet in `MOESM3` indicates that its PC1 label was derived after removing `Gejman`, whereas the mouse worksheet uses the full mouse dataset. Accordingly, `MOESM3` is used here mainly as a reference and control for the Saluki processing convention rather than as a replacement training-label source. Precomputed Saluki sequence/regulatory feature blocks were obtained from the associated public dataset (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)). Because `MOESM2` already contains transformed values, we do not reapply a log transform and instead begin from sample-wise standardization.

Sample metadata were parsed from sample names and local curation tables, including `species`, `study`, `method`, `cell type`, and `replicate`. The full human analysis includes 54 samples, 19 studies, and 5 experimental methods. The mouse analysis includes 27 samples, 17 studies, and 3 experimental methods. The no-Gejman human subset includes 39 samples and 18 studies. Label auditing uses species-specific coverage-filtered genes, while the global prior benchmark uses 12,916 common human genes. The 0.10 ortholog-informed shrinkage target uses 12,307 shared genes, of which 10,768 have a one-to-one ortholog pair. Human-mouse ortholog mapping was obtained from Ensembl release 115 comparative genomics resources (Dyer et al. 2025).

![](figures/main/Fig02_dataset_summary_cn.png)

**Fig. 2** Data sources, sample structure, and analysis gene universes. a, Numbers of samples, studies, and experimental methods in the human, mouse, and no-Gejman human subsets. b, Gene coverage retained after preprocessing. c, Gene or ortholog-pair universes used for MOESM3 controls, the global prior benchmark, ortholog concordance, and 0.10 target shrinkage. d, Roles of `MOESM2`, `MOESM3`, and Ensembl ortholog resources as primary reconstruction input, Saluki reference/control, and cross-species comparison/mapping, respectively

### 3.2 Sample-level principal component analysis and diagnosis of compendium heterogeneity

Before constructing supervision labels, we examined the dominant structure of the compendium at the sample level using principal component analysis (PCA). For the raw transformed matrix, we applied iterative low-rank PCA imputation only to missing values for visualization, following the general low-rank treatment of incomplete high-dimensional matrices (Tipping and Bishop 1999; Troyanskaya et al. 2001). For the processed matrix, we applied the same preprocessing steps used in label reconstruction, namely sample-wise z-scoring, iterative PCA imputation, and quantile normalization (Bolstad et al. 2003). We then calculated the first two principal components in sample by gene space and colored the samples by study or method. Supplementary Table S3 gives the matrix definitions, imputation parameters, and intended use of each PCA view.

This exploratory PCA tested whether major variation in the matrix already contained pronounced study or method structure. Such clustering would make label provenance relevant to downstream model performance. In the human compendium, preprocessing changed but did not eliminate this structure, motivating the subsequent leave-one-study-out analysis. Mouse results are retained in Online Resource 1 because the two species show different heterogeneity patterns.

![](figures/main/Fig03_human_pca_panel.png)

**Fig. 3** Human sample-level PCA diagnosis. a, Sample-level PCA of the raw transformed matrix. b, Sample-level PCA after sample-wise z-scoring, iterative PCA imputation, and quantile normalization. Each point represents one experimental sample; color denotes experimental method and marker shape denotes study. Both panels reveal study/method-related structure in the human compendium. Preprocessing changes the PC geometry, but study/method structure remains; the panels use separate PC scales and are interpreted qualitatively

### 3.3 Consensus-label reconstruction

Consensus-label reconstruction was performed on the gene by sample matrix. For each species or study-removal subset, we applied coverage filtering, sample-wise z-scoring, iterative PCA imputation, and quantile normalization in sequence. We then extracted PC1 from the processed matrix and aligned its sign to the gene-wise mean. The resulting PC1 was used as the local consensus label.

We use this local procedure as an auditable and rerunnable reconstruction, while acknowledging that unpublished details of Saluki processing cannot be reproduced exactly. The reconstructed full human PC1 correlates with Saluki human PC1 at 0.948, and the reconstructed mouse PC1 correlates with the mouse reference at 0.958. For 0.10 ortholog-informed target shrinkage, we use coverage thresholds of at least 10 non-missing samples for human and 5 for mouse. Both labels are reconstructed locally from `MOESM2`; "Saluki-like" denotes only these coverage thresholds.

### 3.4 Leave-one-study-out influence analysis

For each study $s$, we removed all of its samples and reran coverage filtering, z-scoring, PCA imputation, quantile normalization, and PC1 extraction. If $L_{\mathrm{full}}$ is the full label and $L_{-s}$ is the reconstructed label after removal, PC1 stability is defined as

```equation
number: 1
S_s = \operatorname{corr}(L_{\mathrm{full}}, L_{-s})
```

Lower values of Eq. (1) indicate greater geometric influence. For a comparison label $R$, agreement change is defined as

```equation
number: 2
A_s(R) = \operatorname{corr}(L_{-s}, R) - \operatorname{corr}(L_{\mathrm{full}}, R)
```

Equation (2) is computed on the paired common genes for each comparison.

The audit has two stages. PC1 stability first screens influential studies without using `MOESM3`, Saluki human PC1, or downstream prediction scores. We then assess the direction of change using agreement with released Saluki human PC1 and human-mouse ortholog concordance. Because Saluki documentation identifies `Gejman`, the first comparison measures recovery of the published processing choice and does not establish biological invalidity. Ortholog concordance provides a cross-species check with a different evidential role.

Because `Gejman` contributes 15 of the 54 human samples, we added a size-matched null. In each of $B=500$ replicates, 15 samples were drawn without replacement from the 39 non-Gejman samples, and the full reconstruction pipeline was rerun. Let $N_{\mathrm{extreme}}$ denote the number of null values at least as extreme as the observed statistic. The one-sided empirical p value was

```equation
number: 3
p_{\mathrm{emp}} = \frac{1 + N_{\mathrm{extreme}}}{B + 1}
```

We used the lower tail for PC1 stability and the upper tail for both agreement gains. This separates the amount of geometric change caused by removing many samples from the direction of that change on the Saluki-processing and ortholog-comparison axes.

Three sensitivity analyses addressed pipeline dependence. First, the primary dynamic coverage analysis was repeated on a fixed gene universe that satisfied the coverage threshold in the full matrix and every leave-one-study-out matrix. Second, we varied the minimum observed samples per gene (3, 5, or 10), PCA-imputation rank (3, 5, or 10), and replaced iterative PCA imputation with sample-wise median imputation. Third, we reduced sample-count imbalance by either multiplying samples from study $s$ by $1/\sqrt{n_s}$ before PCA, which equalizes each study's total squared PCA weight, or collapsing each study to its mean processed profile before PC1 extraction. Complete definitions and results are provided in Supplementary Section S5.

### 3.5 Ortholog concordance as a cross-species label check

To avoid keeping label diagnosis entirely internal to the human measurements, we used human-mouse one-to-one ortholog concordance as a cross-species comparison (Koonin 2005; Gabaldon and Koonin 2013). Specifically, human full PC1, human no-Gejman PC1, and Saluki human PC1 were mapped to one-to-one ortholog genes and compared with reconstructed mouse PC1 or released Saluki mouse PC1 by Pearson and Spearman correlation. The analysis was repeated on all one-to-one ortholog pairs and on a high-confidence subset. High confidence was defined by the binary `is_high_confidence` field in Ensembl Compara release 115. The full analysis retained all `ortholog_one2one` pairs, and the subset analysis retained pairs with `is_high_confidence == 1`.

Saluki agreement and ortholog concordance serve different purposes. An isolated increase in Saluki agreement would be consistent with mechanical alignment to its processing convention. Improvement in both supports movement toward a conserved mRNA-stability axis, without identifying a molecular mechanism or making the underlying data sources independent. We interpret the `Gejman` result from this combined pattern.

### 3.6 Construction and cross-target validation of the 0.10 ortholog-informed shrinkage target

We use a small ortholog constraint to construct a human-dominant target that represents conserved mammalian stability. This step changes the supervised target while leaving the human feature definitions unchanged. We reconstructed human no-Gejman PC1 and mouse PC1 under Saluki-like coverage, z-scored both labels, and mapped the mouse label to human genes through one-to-one orthologs. For genes with an ortholog signal, let $H$ denote the standardized human label and $M$ the mapped standardized mouse label. The general shrinkage target and the primary operating point are

```equation
number: 4
T_{\lambda} = (1-\lambda)H + \lambda M, \qquad T_{0.10} = 0.90H + 0.10M
```

Genes without an ortholog signal retained $H$. The resulting target was re-z-scored and sign-aligned to Saluki human PC1. We refer to it as the **0.10 ortholog-informed shrinkage target**; it is a constructed target rather than a new experimental human half-life measurement. `ortholog_regularized` is retained only in code and result filenames.

A measurement-error model gives the statistical interpretation. Let $Z$ be a shared mammalian stability axis, $\Delta$ a species-specific shift, and $\varepsilon_H$ and $\varepsilon_M$ compendium-specific errors. The observed labels are represented as

```equation
number: 5
H = Z + \varepsilon_H, \qquad M = Z + \Delta + \varepsilon_M
```

Substitution of Eq. (5) into Eq. (4) gives

```equation
number: 6
T_{\lambda} = Z + \lambda\Delta + (1-\lambda)\varepsilon_H + \lambda\varepsilon_M
```

Writing $\sigma_H^2=\operatorname{Var}(\varepsilon_H\mid\Delta)$, $\sigma_M^2=\operatorname{Var}(\varepsilon_M\mid\Delta)$, and $C_{HM}=\operatorname{Cov}(\varepsilon_H,\varepsilon_M\mid\Delta)$, the conditional error variance is

```equation
number: 7
V_{\lambda} = (1-\lambda)^2\sigma_H^2 + \lambda^2\sigma_M^2 + 2\lambda(1-\lambda)C_{HM}
```

Under independent errors, $C_{HM}=0$. If the errors also have zero conditional means, the mean squared error relative to $Z$ is

```equation
number: 8
\operatorname{E}\left[(T_{\lambda}-Z)^2\mid\Delta\right] = V_{\lambda} + \lambda^2\Delta^2
```

Increasing $\lambda$ can therefore reduce the contribution of human-label noise but increases potential species-shift bias by $\lambda\Delta$. Equations (4)-(8) describe the mapped-gene mixture before the final affine re-z-scoring step; re-z-scoring preserves correlations but changes the numerical error scale. This is a working bias-variance model, and neither error independence nor theoretical optimality of $\lambda=0.10$ is assumed.

The target was evaluated in three separate ways. First, within a common one-to-one ortholog universe, we calculated Pearson, Spearman, MSE, RMSE, and MAE against the human and mouse labels. Second, we compared the error in $T_{0.10}-H$ with study-to-study variability among no-Gejman human studies on the same z-scored label scale. Third, we performed cross-target evaluation. Models were trained on either human no-Gejman PC1 or $T_{0.10}$, then evaluated against the original human no-Gejman PC1. Both analyses used the same 12,307 genes, 1802 human-only `compact_all` features, 10 folds by 3 random seeds, and model parameters. A gene-level paired bootstrap provided a 95% CI for the difference in Pearson correlation. This design tests whether target shrinkage sacrifices predictability of the original human target.

We fixed $\lambda=0.10$ for the primary analysis to limit the mouse contribution to 10%. It was not selected by maximizing prior-enhanced cross-validation. Sensitivity results for $\lambda=0.05$, 0.10, and 0.30 are reported in full in the supplement. Larger values produce higher predictability after alignment with mouse priors, but also a larger target displacement, and are therefore not used for the primary claim.

### 3.7 Sequence/regulatory feature construction and evaluation of the human-only baseline

The human-only route uses two feature layers while keeping the target unchanged. The first layer is a 526-dimensional base sequence feature set derived from Ensembl release 115 representative transcripts, covering 5' untranslated region (5'UTR), coding sequence (CDS), and 3'UTR lengths, GC content, codon frequencies, region-level 3-mers, and 3'UTR 4-mers. The second layer is the 1802-dimensional `compact_all` feature set, formed by augmenting base sequence features with precomputed CWCS, SeqWeaver, and DeepRiPe regulatory blocks from the Saluki datapack (Agarwal et al. 2015; Park et al. 2021; Ghanbari and Ohler 2020).

All human-only benchmarks use Saluki human PC1 as the primary supervision target and are evaluated by gene-level out-of-fold prediction. The main results use XGBoost/CUDA with `tree_method=hist`, `device=cuda`, `max_depth=6`, `learning_rate=0.02`, `min_child_weight=4`, `subsample=0.9`, `colsample_bytree=0.75`, and `reg_lambda=1.0`. The human-only benchmark allowed up to 5,000 trees with 120-round early stopping, whereas the global prior benchmark allowed up to 3,000 trees with 100-round early stopping. Within every outer fold, `max(512, 10% of outer-training genes)` was reserved from the outer training set as an inner validation subset. The outer held-out fold was never used for early stopping and was evaluated only after model fitting. The main text reports 5-fold OOF results, and a 10-fold by 3-random-seed evaluation is used as a split-sensitivity robustness check.

We evaluated all models by gene-level OOF prediction rather than recreating Saluki's fixed test split. This unified design compares label processing, prior input, and controls within one feature table and label universe. The human-only route provides a transparent and reproducible baseline.

### 3.8 Cross-species prior-enhanced evaluation

The global cross-species transfer benchmark retains Saluki human PC1 as the target and uses each human gene as the modeling unit. We concatenated six mouse-related columns to the human `compact_all` features: two prior values (`mouse_pc1` and released Saluki mouse PC1), two availability indicators, and two high-confidence ortholog indicators. Here, `mouse_pc1` is our locally reconstructed mouse prior, whereas released Saluki mouse PC1 is a separately processed mouse label used as a prior. The two values carry related but non-identical cross-species stability information. Both priors are derived exclusively from mouse measurements and ortholog mapping, then fixed before human cross-validation. Human target values from training or held-out folds do not enter their construction. The resulting inputs are external cross-species covariates, so this setting is reported as cross-species transfer.

The indicator columns distinguish the prior value from its availability and mapping confidence. If a human gene lacks a mouse prior, its temporary prior value and availability indicator are both set to 0. An available prior has an indicator of 1. The high-confidence columns similarly record whether the ortholog has `is_high_confidence == 1` in Ensembl Compara. We used the full common-gene universe to report performance across the complete benchmark, including genes with missing priors.

Thus, the training and evaluation set for the fixed-target prediction benchmark is the same `global prediction universe`, containing 12,916 human genes with Saluki human PC1 fixed as the supervised target. Human-only `compact_all`, single-prior models, dual-prior models, and shuffled-prior controls are all evaluated on these 12,916 genes by gene-level OOF prediction. The `both-priors` subset contains 11,107 genes with both `mouse_pc1` and released Saluki mouse PC1 available. It is used mainly for prior-coverage interpretation and residual analysis, not as the main training universe. The 0.10 target-shrinkage analysis uses a separate 12,307-gene universe, and its one-to-one ortholog label-distance comparison uses 10,768 pairs. These sets overlap but are not identical.

### 3.9 Control experiments, evaluation metrics, and usage boundaries

We used four controls. First, permutation tests whether prior-enhanced gain depends on gene-specific prior values rather than simply on extra columns or missingness structure. Within each fold, `mouse_pc1` and released Saluki mouse PC1 values were shuffled only among genes that originally had the respective value, while availability and high-confidence indicators were held fixed. Second, residual analysis tests whether the transfer result merely repeats the priors. It was restricted to the `both_priors_available` subset (N = 11,107). Within each outer fold, RidgeCV was fitted on outer-training genes; its fitted values defined outer-training residuals, while held-out prior predictions were generated separately. XGBoost used only the outer-training residuals and human `compact_all` features, so the final held-out predictions remained OOF. Third, no-direct-prior ablation tests whether prediction of the 0.10 ortholog-informed shrinkage target depends entirely on the locally reconstructed mouse PC1 used to construct that target. Fourth, the cross-target evaluation in Section 3.6 tests whether target shrinkage impairs OOF prediction of the original human label.

Model performance is reported by Pearson, Spearman, and coefficient of determination (R2). Label distance is reported by Pearson, Spearman, MSE, RMSE, MAE, median absolute error, and p95 absolute error. Unless otherwise noted, all prediction results are OOF rather than in-sample fits. Mean ± SD across random seeds describes split sensitivity and is not treated as a sampling confidence interval. Paired gene or ortholog bootstrap intervals use 2,000 resamples for the cross-target comparison and 5,000 resamples for the other key claims (Supplementary Tables S14-S15). Online Resource 1 retains the full rankings, multi-seed summaries, prior ablations, label-distance analyses, study-noise comparisons, and bootstrap outputs.

Finally, the prediction-usage boundary in this study has three levels. The base sequence model requires only explicitly partitioned 5'UTR, CDS, and 3'UTR sequence, and can therefore be applied to new transcripts or designed sequences, but it is not equivalent to the full `compact_all` feature setting. `compact_all` applies to Ensembl human genes represented in the precomputed Saluki feature datapack. The prior-enhanced model additionally requires mouse ortholog priors; if both mouse priors are missing, the prediction should fall back to the human-only compact model.

Generative AI tools were used for manuscript-structure and language-drafting support, code review, and figure-layout suggestions; they were not used to generate analytical data or final figures. All implementations, numerical results, references, figures, and final wording were verified by the authors, who take responsibility for the submitted work and its reproducibility.

**Table 2** Prediction settings, targets, inputs, and interpretation boundaries.

| Modeling setting | Target label | Input features | Intended interpretation and boundary |
| --- | --- | --- | --- |
| Base sequence | Saluki human PC1 | 526 sequence-statistics features | Applicable to new sequences, but only at the base-sequence level |
| Compact human-only | Saluki human PC1 | 1802-dimensional `compact_all` | Pure human sequence/regulatory benchmark and transparent human-only baseline |
| Global prior-enhanced | Saluki human PC1 | `compact_all` + two mouse priors + flags | Cross-species transfer using mouse-side covariates external to the human target |
| 0.10 ortholog-informed target shrinkage | `T_0.10` shrinkage target | Human-only or prior-enhanced features | Human-dominant mammalian stability target; cross-target evaluation tests original-human-label predictability |
| Shuffled/no-direct controls | Same as corresponding parent target | Shuffled priors or removal of the direct prior value | Test whether the gain arises from gene-specific cross-species signal |

## 4. Results

### 4.1 Study-aware auditing separates sample-count effects from the direction of label change

The reconstructed full human label correlated with Saluki human PC1 at 0.948, indicating that the local pipeline reproduced its principal structure. In the primary sample-weighted analysis, leave-one-study-out PC1 stability ranked `Gejman` first among 19 studies: removing its 15 samples gave r = 0.955, compared with r = 0.973 for the next-ranked study. `Gejman` remained first under dynamic and fixed gene universes, coverage thresholds of 3, 5, or 10 observed samples, PCA-imputation ranks of 3, 5, or 10, and sample-wise median imputation.

The size-matched analysis changed the interpretation of this geometric ranking. Across 500 complete pipeline reruns that each removed 15 non-Gejman samples, median PC1 stability was 0.944 (95% null interval, 0.932-0.956). The observed Gejman value of 0.955 was not unusually low (one-sided empirical p = 0.964). Moreover, `Gejman` fell from rank 1 to rank 3 under both equal-study PCA weighting and study-mean collapse. Its raw influence on PC1 geometry is therefore partly attributable to sample count; the data do not identify a sample-count-adjusted statistical outlier.

The direction of the change, however, was study-specific. On a fixed set of 13,265 paired common genes, removing `Gejman` increased Pearson agreement with Saluki human PC1 from 0.9515 to 0.9829, a gain of 0.0314 (95% bootstrap confidence interval, 0.0298-0.0330). None of the 500 size-matched removals produced a gain as large (null median, -0.0926; one-sided empirical p = 0.002). Fig. 4b instead shows the primary dynamic-coverage estimate computed on each removal's available paired genes; for `Gejman`, that estimate is +0.0349. This processing-reference comparison was not used to generate the stability ranking, and no downstream model score entered the audit. Because Saluki human PC1 was released after Gejman removal, the result measures recovery of that published processing choice rather than independent validation.

![](figures/main/Fig04_human_study_influence.png)

**Fig. 4** Primary sample-weighted leave-one-study-out audit of human studies. a, Points show PC1 stability for the eight studies with the largest effect on full-label geometry; lower values indicate greater influence, and the dashed line marks r = 1. This ranking does not use Saluki PC1. b, Agreement change is `Pearson(leave-one-study-out PC1, Saluki human PC1) - Pearson(full PC1, Saluki human PC1)`, calculated under the primary dynamic-coverage estimator on the paired genes available after each removal. The `Gejman` value is +0.0349; the fixed 13,265-gene paired estimate used for bootstrap uncertainty in the text is +0.0314. Size-matched, fixed-universe, preprocessing, and study-balanced analyses are reported in Supplementary Section S5 and Supplementary Fig. S4

Human-mouse orthologs supplied a second comparison. Across 12,592 one-to-one pairs, Pearson concordance increased from 0.741 to 0.755 after removing `Gejman`, a gain of 0.0147 (95% bootstrap confidence interval, 0.0111-0.0181); Spearman concordance increased from 0.730 to 0.744. The gain remained 0.0146 in the Ensembl high-confidence subset. None of the 500 size-matched removals matched the full-set ortholog gain (null median, -0.0616; one-sided empirical p = 0.002). Sample count contributes to how far the label moves. Removing `Gejman` specifically moves it toward both the Saluki processing reference and conserved ortholog structure.

![](figures/main/Fig05_ortholog_summary_cn.png)

**Fig. 5** Ortholog concordance analysis of the no-Gejman label. a-b, Each point represents one one-to-one human-mouse ortholog pair; the x-axis is human consensus PC1 and the y-axis is mouse consensus PC1. a uses full human PC1 and b uses no-Gejman human PC1; red lines are ordinary least-squares fits for visual reference. c, Pearson and Spearman comparisons across full PC1, no-Gejman PC1, and their high-confidence ortholog subsets. d, Paired Pearson and Spearman gains after removing `Gejman` in the full one-to-one and high-confidence subsets

### 4.2 A human-only compact sequence/regulatory model defines a stable fixed-target baseline

For the fixed Saluki human PC1 target, the primary global analysis gave Pearson = 0.745, Spearman = 0.737, and R2 = 0.554 for `compact_all` across 12,916 genes. In a matched common-universe comparison, Ridge, histogram gradient boosting, GPU XGBoost, and `compact_all` achieved Pearson values of 0.698, 0.726, 0.737, and 0.743, respectively (Fig. 6a). The 0.002 difference between the matched and primary results reflects their cross-validation and training configurations.

Regulatory-block comparisons on a separate 13,532-gene feature-complete universe showed heterogeneous individual contributions. CWCS and DeepRiPe each improved over the base setting, SeqWeaver alone was nearly unchanged, and their combined `compact_all` representation was highest (Fig. 6b). On the matched 12,916-gene common universe, `compact_all` correlated at 0.721 with full human PC1, 0.737 with human no-Gejman PC1, and 0.743 with Saluki human PC1 (Fig. 6c). Three-fold gene-subset analyses retained the same overall performance range without concentrating the result in one selected subset (Fig. 6d). These analyses establish a reproducible human-only OOF baseline. They do not reproduce Saluki's fixed train-test split.

![](figures/main/Fig06_sequence_progression_cn.png)

**Fig. 6** Evidence for the human-only fixed-target benchmark. a, Matched 5-fold algorithm and feature comparison on 12,916 genes across Ridge, histogram gradient boosting, GPU XGBoost, and `compact_all`. b, Five-fold regulatory-block comparison on a separate 13,532-gene feature-complete universe; individual contributions are heterogeneous, whereas the combined `compact_all` setting is highest. c, Matched 5-fold label-definition sensitivity on 12,916 genes across full human PC1, human no-Gejman PC1, and Saluki human PC1. d, Three-fold gene-subset sensitivity across all-common, all-ortholog, high-confidence ortholog, and top-regulatory-33% subsets; subset values are not extrapolated as global performance

### 4.3 Cross-species transfer produces a gene-specific prior gain under a fixed human target

We next kept Saluki human PC1 fixed and changed only the input information. Human-only `compact_all` reached Pearson = 0.745. Adding the locally reconstructed mouse PC1 increased Pearson to 0.828, adding only the released Saluki mouse PC1 as a prior gave 0.827, and adding both priors gave 0.829 (Fig. 7a). Each mouse-side prior therefore carried substantial, partly redundant cross-species information.

Repeated 10-fold evaluation with three random states gave 0.748 ± 0.001 for human-only `compact_all`, 0.748 ± 0.001 for the shuffled-prior control, and 0.830 ± 0.001 for the real dual-prior model; ± denotes the SD across random states. The coverage-aware fallback rule, which returns genes without priors to the human-only prediction, gave 0.832 ± 0.001. Of the 12,916 benchmark genes, 11,107 had both priors, 462 had only the reconstructed mouse PC1, and 1,347 had neither (Fig. 7b-c). These coverage groups define which model can be used for a given gene.

![](figures/main/Fig07_prior_summary_cn.png)

**Fig. 7** Evidence for cross-species transfer under a fixed human target. a, Global benchmark across human-only, reconstructed-mouse-PC1, released-Saluki-mouse-PC1, and dual-prior inputs. b, Repeated 10-fold controls; error bars are across-seed SD. Shuffling preserves prior availability but disrupts gene-prior identity, and the gain disappears. c, Prior availability in the common human gene universe. d, Absolute model performance in the full one-to-one and high-confidence ortholog subsets

Permutation control identified gene-prior correspondence as the source of the gain. In the primary 5-fold analysis, the real dual-prior model reached Pearson = 0.829 and R2 = 0.688, whereas fold-wise prior shuffling returned it to Pearson = 0.745 and R2 = 0.553, nearly matching the human-only baseline. Prior availability and high-confidence indicators were held fixed, so only the mapping between each gene and its prior value was disrupted. After averaging the three held-out predictions per gene, the real-prior and shuffled-prior correlations were 0.8324 and 0.7531, respectively, for a paired gain of 0.0794 (95% confidence interval, 0.0735-0.0855).

Residual decomposition showed that the prior was strong but incomplete. On the 11,107 genes with both priors, prior-only RidgeCV achieved Pearson = 0.792 and $R^2=0.627$, whereas human `compact_all` alone achieved Pearson = 0.739 and $R^2=0.544$. Adding an OOF human-feature correction learned from outer-training residuals increased the final result to Pearson = 0.826 and $R^2=0.675$. We quantified the fraction of residual variance explained by human features as

```equation
number: 9
R_{\mathrm{remaining}}^2 = \frac{R_{\mathrm{final}}^2 - R_{\mathrm{prior}}^2}{1 - R_{\mathrm{prior}}^2}
```

Equation (9) gives 12.8%. This subset analysis decomposes the signal and does not replace the 12,916-gene benchmark.

### 4.4 Weak ortholog-informed target shrinkage remains human-dominant without a detectable reduction in human-label predictability

This analysis changed the target. Across all 12,307 genes in the target universe, the shrunken target correlated with human no-Gejman PC1 at 0.9982 (95% bootstrap confidence interval, 0.9981-0.9983). Within the 10,768 mapped one-to-one orthologs, human no-Gejman PC1 and locally reconstructed mouse PC1 correlated at 0.789 before shrinkage (RMSE = 0.663; MAE = 0.512); after applying $\lambda=0.10$, the target correlated with the human label at 0.9979 (RMSE = 0.065; MAE = 0.050) and with mouse PC1 at 0.827. The operation therefore makes a small ortholog-informed adjustment to a predominantly human label. It is a target-construction step; no penalty is added to the prediction model.

**Table 3** Same-scale comparison of the 0.10 ortholog-informed shrinkage target with human/mouse labels and study-level noise. All error metrics are computed in z-scored label units. Study-level rows report medians.

| Comparison | N or number of comparisons | Pearson | RMSE | MAE | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| Shrinkage target vs human no-Gejman PC1 | 10,768 genes | 0.998 | 0.065 | 0.050 | Human label structure remains dominant |
| Shrinkage target vs mouse PC1 | 10,768 genes | 0.827 | 0.601 | 0.464 | Mild shift toward a cross-species conserved direction |
| Human no-Gejman PC1 vs mouse PC1 | 10,768 genes | 0.789 | 0.663 | 0.512 | Human and mouse labels already differ before shrinkage |
| No-Gejman study-study label pairs (median) | 153 study pairs | 0.530 | 0.954 | 0.720 | Variability among ordinary studies |
| Single-study vs human no-Gejman PC1 (median) | 18 studies | 0.742 | 0.715 | 0.543 | Deviation of a single study from the human consensus |

The target displacement was also small relative to study-level variability. Pairwise no-Gejman study labels had median Pearson = 0.530, RMSE = 0.954, and MAE = 0.720; individual studies relative to human no-Gejman PC1 had median Pearson = 0.742, RMSE = 0.715, and MAE = 0.543. Cross-target evaluation did not detect a reduction in original-human-label predictability. Direct training on human no-Gejman PC1 achieved r = 0.7476 ± 0.0010 against that label, whereas training on the shrunken target and evaluating against the same human label achieved r = 0.7480 ± 0.0010. The paired difference was +0.0003 (95% confidence interval, -0.0006 to 0.0012); because the interval includes zero and no equivalence margin was prespecified, this result does not establish formal equivalence.

![](figures/main/Fig08_ortholog_regularized_target_cn.png)

**Fig. 8** The 0.10 ortholog-informed target shrinkage retains human-dominant structure and shows no detectable reduction in original-human-label predictability. a, The shrunken target lies almost on the diagonal relative to human no-Gejman PC1. b, Its correlation with mouse PC1 remains lower than its correlation with the human label. c, Its displacement is smaller than study-study and study-consensus differences on the same z-scored scale. d, Models trained on the human or shrunken target use the same 12,307 genes, human-only features, folds, seeds, and parameters, and are both evaluated against the original human target. The point is the Pearson difference after pooling the three OOF prediction vectors, and the error bar is its 95% paired gene-bootstrap confidence interval; values in the text are mean ± SD across seeds

Prediction of this target provided a separate conserved-target test. It remained correlated with Saluki human PC1 at 0.987. On the same 12,307 genes, human-only prediction achieved 0.754 ± 0.001 and dual-prior prediction achieved 0.855 ± 0.001. Removing the locally reconstructed mouse PC1 used in target construction, while retaining the released Saluki mouse PC1 and indicators, still gave 0.852 ± 0.001; retaining only missingness and high-confidence indicators returned performance to 0.753 ± 0.001. These results support a human-dominant mammalian stability target whose conserved component is more predictable with external ortholog information. Because this analysis uses a different target, it is reported separately from the fixed Saluki human PC1 ranking.

## 5. Discussion

This study brings label construction into benchmark evaluation. Study influence, external cross-species covariates, and target shrinkage are evaluated with controls matched to the quantity that changes.

The sample-size analysis refines the interpretation of `Gejman`. Its 15 samples give it substantial leverage in ordinary sample-weighted PCA, and both study-balanced estimators reduced its geometric rank. The data do not establish a sample-count-adjusted outlier or identify a unique experimental cause. Yet removing 15 other samples did not reproduce the positive Saluki-agreement or ortholog-concordance changes. The first comparison measures recovery of the published processing choice, while the second provides a cross-species check. Their joint pattern supports a directional study effect superimposed on sample-count leverage. In other compendia, geometric influence should likewise be distinguished from the direction of label change.

The fixed-target analyses make a second distinction between human representation and cross-species transfer. The compact model defines what can be obtained from the available human sequence/regulatory feature table, whereas the mouse priors add measured ortholog-level covariates. Their gain may reflect conserved determinants of transcript stability, shared measurement structure, or both; the present design does not assign a molecular mechanism. The permutation control does show that availability flags and feature count are insufficient: gene-specific correspondence is required. Residual analysis further shows that human features retain information after the strong mouse prior is modeled.

The $\lambda=0.10$ analysis is a shrinkage estimator for a noisy cross-study target. It moves the human consensus slightly toward a correlated ortholog signal, with no detectable reduction in predictive performance against the original human label. We interpret the resulting score as an explicitly defined mammalian stability target. It is not a new experimental half-life measurement. The same analysis could extend to other multi-study molecular phenotypes with a correlated external signal. Such applications should report label displacement, cross-target performance, and input coverage together.

## 6. Limitations

The study has three limitations. First, all label audits and principal benchmarks use one public compendium. The size-matched null controls sample count but removes independent samples rather than whole correlated pseudo-studies, and study balancing changes the geometric ranking; an independent half-life resource is therefore needed to test portability. Second, mouse priors are external to the human target but are not mechanistically independent measurements of every stability determinant. Cross-species transfer and 0.10 target shrinkage change input information and target definition, respectively, and cannot be placed on one leaderboard. The empirical choice $\lambda=0.10$ also requires independent calibration. Third, we did not retrain Saluki or other large sequence models, and the implementation is not a one-click predictor for arbitrary unannotated sequences because base-sequence, compact-regulatory, and prior-enhanced predictions require different input resources.

## 7. Conclusion

Study-aware auditing showed that `Gejman` has sample-count-sensitive leverage on human consensus-label geometry, while its positive Saluki and ortholog agreement gains were not reproduced by 500 size-matched removals. Under fixed Saluki human PC1, a compact human-only model provided a reproducible baseline and gene-specific mouse priors increased repeated 10-fold performance to Pearson = 0.830; permutation removed the gain. Ortholog-informed target shrinkage at $\lambda=0.10$ remained nearly identical to human no-Gejman PC1, and cross-target evaluation detected no reduction in its out-of-fold predictability. Together, these results provide an auditable framework for mammalian mRNA half-life benchmarking.

## 8. Statements and Declarations

### Funding

No specific funding was received for this work.

### Author Contributions

Wenzhuo Wang contributed to conceptualization, methodology, software, formal analysis, data curation, validation, visualization, and writing of the original draft. Ying Shao contributed to conceptualization, methodology, validation, supervision, project administration, and review and editing of the manuscript. Both authors reviewed and approved the final manuscript.

### Data Availability

The half-life matrices analyzed here are from the public supplementary material of the Saluki study; the precomputed Saluki feature datapack is available from Zenodo (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)); and ortholog mappings are from Ensembl comparative genomics resources (Agarwal and Kelley 2022; Dyer et al. 2025). Result tables, figure source data, and reproduction notes are publicly available in the versioned repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label).

### Code Availability

Code, figure-generation scripts, and result-reproduction workflows are publicly available in [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label); the audited two-author manuscript state is identified by the release tag `mRNA-PC1-label-v1.3`.

### Competing Interests

The authors declare that they have no known competing financial interests or personal relationships that could have influenced the work reported in this paper.

### Ethics Approval

This study is based on public data and computational analysis and does not involve new animal experiments, human participant recruitment, or clinical intervention. Ethics approval is therefore not applicable.

## 9. References

Agarwal V, Bell GW, Nam JW, Bartel DP (2015) Predicting effective microRNA target sites in mammalian mRNAs. eLife 4:e05005. https://doi.org/10.7554/eLife.05005

Agarwal V, Kelley DR (2022) The genetic and biochemical determinants of mRNA degradation rates in mammals. Genome Biology 23(1):245. https://doi.org/10.1186/s13059-022-02811-x

Alipanahi B, Delong A, Weirauch MT, Frey BJ (2015) Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning. Nature Biotechnology 33:831-838. https://doi.org/10.1038/nbt.3300

Avsec Z, Agarwal V, Visentin D, Ledsam JR, Grabska-Barwinska A, Taylor KR, Assael Y, Jumper J, Kohli P, Kelley DR (2021) Effective gene expression prediction from sequence by integrating long-range interactions. Nature Methods 18:1196-1203. https://doi.org/10.1038/s41592-021-01252-x

Bolstad BM, Irizarry RA, Astrand M, Speed TP (2003) A comparison of normalization methods for high density oligonucleotide array data based on variance and bias. Bioinformatics 19(2):185-193. https://doi.org/10.1093/bioinformatics/19.2.185

Breiman L (2001) Random forests. Machine Learning 45:5-32. https://doi.org/10.1023/A:1010933404324

Cetnar DP, Hossain A, Vezeau GE, Salis HM (2024) Predicting synthetic mRNA stability using massively parallel kinetic measurements, biophysical modeling, and machine learning. Nature Communications 15:9601. https://doi.org/10.1038/s41467-024-54059-7

Chen T, Guestrin C (2016) XGBoost: A scalable tree boosting system. Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining 785-794. https://doi.org/10.1145/2939672.2939785

Duan J, Shi J, Ge X, Dolken L, Moy W, He D, Shi S, Sanders AR, Ross J, Gejman PV (2013) Genome-wide survey of interindividual differences of RNA stability in human lymphoblastoid cell lines. Scientific Reports 3:1318. https://doi.org/10.1038/srep01318

Duffy EE, Rutenberg-Schoenberg M, Stark CD, Kitchen RR, Gerstein MB, Simon MD (2015) Tracking distinct RNA populations using efficient and reversible covalent chemistry. Molecular Cell 59(5):858-866. https://doi.org/10.1016/j.molcel.2015.07.023

Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al (2025) Ensembl 2025. Nucleic Acids Research 53(D1):D948-D957. https://doi.org/10.1093/nar/gkae1071

Friedel CC, Dolken L, Ruzsics Z, Koszinowski UH, Zimmer R (2009) Conserved principles of mammalian transcriptional regulation revealed by RNA half-life. Nucleic Acids Research 37(17):e115. https://doi.org/10.1093/nar/gkp542

Friedman JH (2001) Greedy function approximation: A gradient boosting machine. The Annals of Statistics 29(5):1189-1232. https://doi.org/10.1214/aos/1013203451

Gabaldon T, Koonin EV (2013) Functional and evolutionary implications of gene orthology. Nature Reviews Genetics 14(5):360-366. https://doi.org/10.1038/nrg3456

Garneau NL, Wilusz J, Wilusz CJ (2007) The highways and byways of mRNA decay. Nature Reviews Molecular Cell Biology 8:113-126. https://doi.org/10.1038/nrm2104

Ghanbari M, Ohler U (2020) Deep neural networks for interpreting RNA-binding protein target preferences. Genome Research 30(2):214-226. https://doi.org/10.1101/gr.247494.118

Grimson A, Farh KKH, Johnston WK, Garrett-Engele P, Lim LP, Bartel DP (2007) MicroRNA targeting specificity in mammals: determinants beyond seed pairing. Molecular Cell 27(1):91-105. https://doi.org/10.1016/j.molcel.2007.06.017

Herzog VA, Reichholf B, Neumann T, Rescheneder P, Bhat P, Burkard TR, Wlotzka W, von Haeseler A, Zuber J, Ameres SL (2017) Thiol-linked alkylation of RNA to assess expression dynamics. Nature Methods 14:1198-1204. https://doi.org/10.1038/nmeth.4435

Hia F, Yang SR, Shichino Y, Yoshinaga M, Murakawa Y, Vandenbon A, Fukao A, Fujiwara T, Landthaler M, Natsume T, Adachi S, Iwasaki S, Takeuchi O (2019) Codon bias confers stability to human mRNAs. EMBO Reports 20(11):e48220. https://doi.org/10.15252/embr.201948220

Ji Y, Zhou Z, Liu H, Davuluri RV (2021) DNABERT: pre-trained Bidirectional Encoder Representations from Transformers model for DNA-language in genome. Bioinformatics 37(15):2112-2120. https://doi.org/10.1093/bioinformatics/btab083

Johnson WE, Li C, Rabinovic A (2007) Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics 8(1):118-127. https://doi.org/10.1093/biostatistics/kxj037

Ke G, Meng Q, Finley T, Wang T, Chen W, Ma W, Ye Q, Liu TY (2017) LightGBM: A highly efficient gradient boosting decision tree. Advances in Neural Information Processing Systems 30. https://proceedings.neurips.cc/paper_files/paper/2017/hash/6449f44a102fde848669bdd9eb6b76fa-Abstract.html

Kelley DR (2020) Cross-species regulatory sequence activity prediction. PLOS Computational Biology 16(7):e1008050. https://doi.org/10.1371/journal.pcbi.1008050

Koonin EV (2005) Orthologs, paralogs, and evolutionary genomics. Annual Review of Genetics 39:309-338. https://doi.org/10.1146/annurev.genet.39.073003.114725

Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA (2010) Tackling the widespread and critical impact of batch effects in high-throughput data. Nature Reviews Genetics 11(10):733-739. https://doi.org/10.1038/nrg2825

Leppek K, Byeon GW, Kladwang W, Wayment-Steele HK, Kerr CH, Xu AF, Kim DS, Topkar VV, Choe C, Rothschild D, et al (2022) Combinatorial optimization of mRNA structure, stability, and translation for RNA-based therapeutics. Nature Communications 13:1536. https://doi.org/10.1038/s41467-022-28776-w

Mauger DM, Cabral BJ, Presnyak V, Su SV, Reid DW, Goodman B, Link K, Khatwani N, Reynders J, Moore MJ, McFadyen IJ (2019) mRNA structure regulates protein expression through changes in functional half-life. Proceedings of the National Academy of Sciences 116(48):24075-24083. https://doi.org/10.1073/pnas.1908052116

Mayr C (2017) Regulation by 3'-untranslated regions. Annual Review of Genetics 51:171-194. https://doi.org/10.1146/annurev-genet-120116-024704

Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. Cell 149(7):1635-1646. https://doi.org/10.1016/j.cell.2012.05.003

Muhar M, Ebert A, Neumann T, Umkehrer C, Jude J, Wieshofer C, Rescheneder P, Lipp JJ, Herzog VA, Reichholf B, et al (2018) SLAM-seq defines direct gene-regulatory functions of the BRD4-MYC axis. Science 360(6390):800-805. https://doi.org/10.1126/science.aao2793

Musaev D, Abdelmessih M, Vejnar CE, Yartseva V, Weiss LA, Strayer EC, Takacs CM, Giraldez AJ (2024) UPF1 regulates mRNA stability by sensing poorly translated coding sequences. Cell Reports 43:114074. https://doi.org/10.1016/j.celrep.2024.114074

Park CY, Zhou J, Wong AK, Chen KM, Theesfeld CL, Darnell RB, et al (2021) Genome-wide landscape of RNA-binding protein target site dysregulation reveals a major impact on psychiatric disorder risk. Nature Genetics 53:166-173. https://doi.org/10.1038/s41588-020-00761-3

Presnyak V, Alhusaini N, Chen YH, Martin S, Morris N, Kline N, Olson S, Weinberg D, Baker KE, Graveley BR, Coller J (2015) Codon optimality is a major determinant of mRNA stability. Cell 160(6):1111-1124. https://doi.org/10.1016/j.cell.2015.02.029

Rabani M, Levin JZ, Fan L, Adiconis X, Raychowdhury R, Garber M, Gnirke A, Nusbaum C, Hacohen N, Friedman N, Amit I, Regev A (2011) Metabolic labeling of RNA uncovers principles of RNA production and degradation dynamics in mammalian cells. Nature Biotechnology 29:436-442. https://doi.org/10.1038/nbt.1861

Rambout X, Maquat LE (2024) Nuclear mRNA decay: regulatory networks that control gene expression. Nature Reviews Genetics 25:679-697. https://doi.org/10.1038/s41576-024-00712-2

Ross J (1995) mRNA stability in mammalian cells. Microbiological Reviews 59(3):423-450. https://doi.org/10.1128/mr.59.3.423-450.1995

Sample PJ, Wang B, Reid DW, Presnyak V, McFadyen IJ, Morris DR, Seelig G (2019) Human 5' UTR design and variant effect prediction from a massively parallel translation assay. Nature Biotechnology 37:803-809. https://doi.org/10.1038/s41587-019-0164-5

Schoenberg DR, Maquat LE (2012) Regulation of cytoplasmic mRNA decay. Nature Reviews Genetics 13(4):246-259. https://doi.org/10.1038/nrg3160

Schofield JA, Duffy EE, Kiefer L, Sullivan MC, Simon MD (2018) TimeLapse-seq: adding a temporal dimension to RNA sequencing through nucleoside recoding. Nature Methods 15(3):221-225. https://doi.org/10.1038/nmeth.4582

Schwalb B, Schulz D, Sun M, Zacher B, Dumcke S, Martin DE, Cramer P, Tresch A (2012) Measurement of genome-wide RNA synthesis and decay rates with Dynamic Transcriptome Analysis (DTA). Bioinformatics 28(6):884-885. https://doi.org/10.1093/bioinformatics/bts052

Schwanhausser B, Busse D, Li N, Dittmar G, Schuchhardt J, Wolf J, Chen W, Selbach M (2011) Global quantification of mammalian gene expression control. Nature 473:337-342. https://doi.org/10.1038/nature10098

Tani H, Mizutani R, Salam KA, Tano K, Ijiri K, Wakamatsu A, Isogai T, Suzuki Y, Akimitsu N (2012) Genome-wide determination of RNA stability reveals hundreds of short-lived noncoding transcripts in mammals. Genome Research 22(5):947-956. https://doi.org/10.1101/gr.130559.111

Tipping ME, Bishop CM (1999) Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B 61(3):611-622. https://doi.org/10.1111/1467-9868.00196

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R, Botstein D, Altman RB (2001) Missing value estimation methods for DNA microarrays. Bioinformatics 17(6):520-525. https://doi.org/10.1093/bioinformatics/17.6.520

Wang X, Lu Z, Gomez A, Hon GC, Yue Y, Han D, Fu Y, Parisien M, Dai Q, Jia G, et al (2014) N6-methyladenosine-dependent regulation of messenger RNA stability. Nature 505:117-120. https://doi.org/10.1038/nature12730

Wu Q, Medina SG, Kushawah G, DeVore ML, Castellano LA, Hand JM, Wright M, Bazzini AA (2019) Translation affects mRNA stability in a codon-dependent manner in human cells. eLife 8:e45396. https://doi.org/10.7554/eLife.45396
