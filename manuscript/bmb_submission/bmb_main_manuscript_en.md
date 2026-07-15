# Study-aware label auditing, weak ortholog regularization, and cross-species transfer for mammalian mRNA half-life prediction

**Article type:** Original research article  
**Authors:** Xu Jin<sup>1,2</sup>, Wenzhuo Wang<sup>3</sup>, Anhui Wang<sup>1,4</sup>, Yuebin Zhang<sup>1,4</sup>, Yingchen Mao*<sup>2</sup>, Dinglin Zhang*<sup>1,4</sup>  
**Affiliations:** <sup>1</sup> Interdisciplinary Research Center for Biology and Chemistry, Liaoning Normal University, Dalian 116029, China; <sup>2</sup> School of Physics and Electronic Technology, Liaoning Normal University, Dalian 116029, China; <sup>3</sup> School of Science, Dalian Maritime University, Dalian 116026, China; <sup>4</sup> Laboratory of Molecular Modeling and Design, State Key Laboratory of Molecular Reaction Dynamics, Dalian Institute of Chemical Physics, Chinese Academy of Sciences, Dalian 116023, China  
**Corresponding authors:** Yingchen Mao, myc@lnnu.edu.cn; Dinglin Zhang, dlzhang@dicp.ac.cn  
**Funding:** National Natural Science Foundation of China (22373101, 22203089, 22573043)  
**Keywords:** mRNA half-life; study heterogeneity; label audit; ortholog regularization; cross-species transfer; XGBoost

## Abstract

Mammalian mRNA half-life compendia combine heterogeneous protocols, cellular states, and processing pipelines, making label quality difficult to separate from model performance. We present a reproducible study-aware framework for consensus-label auditing, fixed-target prediction, and weak ortholog regularization. A leave-one-study-out screen of first principal component (PC1) stability nominated `Gejman` as the most influential human study; Saluki agreement and human-mouse ortholog concordance then validated its effect. Removing `Gejman` increased Saluki agreement from r = 0.948 to r = 0.983. On a paired common-gene universe, the gain was 0.0314 (95% bootstrap CI, 0.0298-0.0330), while ortholog concordance increased from r = 0.741 to r = 0.755 (gain 0.0147; 95% CI, 0.0111-0.0181). With Saluki human PC1 fixed as the target, the human-only model achieved r = 0.748 ± 0.001 in repeated 10-fold evaluation. Mouse ortholog priors increased performance to r = 0.830 ± 0.001; shuffling gene-prior identities returned it to r = 0.748 ± 0.001, with a paired gain of 0.0794 in seed-averaged OOF predictions (95% CI, 0.0735-0.0855). A weakly ortholog-regularized target at λ = 0.10 remained human-dominant (r = 0.998; RMSE = 0.065). Training on this target changed prediction against the original human label by +0.0003 (95% CI, -0.0006 to 0.0012), with no detectable reduction in human-label predictability; prediction with explicit cross-species inputs reached r = 0.855 ± 0.001. The framework separates label auditing, fixed-target cross-species transfer, and weak ortholog regularization.

## 1. Introduction

mRNA stability regulation is an indispensable layer of gene-expression control. Even at the same transcriptional output, different transcripts can exhibit markedly different steady-state abundance and translational output because of differences in degradation rate (Ross 1995; Garneau et al. 2007; Schoenberg and Maquat 2012; Rambout and Maquat 2024). Classical and recent studies have shown that codon composition, ribosome elongation, 3' untranslated region (3'UTR) cis-elements, RNA-binding protein occupancy, RNA structure, and RNA modification all influence mRNA half-life (Meyer et al. 2012; Wang et al. 2014; Presnyak et al. 2015; Wu et al. 2019; Hia et al. 2019; Grimson et al. 2007; Agarwal et al. 2015; Mauger et al. 2019). This is why half-life prediction is relevant not only to basic biology, but also to variant interpretation and RNA therapeutic engineering (Leppek et al. 2022; Sample et al. 2019; Cetnar et al. 2024; Musaev et al. 2024).

At the data level, mammalian mRNA half-life measurements have evolved from small-scale experiments to cross-study compendia. Transcriptional shutoff assays, 4sU/5EU metabolic labeling, nucleoside-recoding approaches, and related high-throughput workflows now provide degradation estimates for thousands to tens of thousands of genes in diverse cellular contexts (Schwanhausser et al. 2011; Rabani et al. 2011; Tani et al. 2012; Schwalb et al. 2012; Duffy et al. 2015; Herzog et al. 2017; Schofield et al. 2018; Muhar et al. 2018). However, the main challenge of such resources is not simply scale, but heterogeneity. Batch effects and study-specific processing can remain prominent in high-throughput matrices even after routine normalization (Johnson et al. 2007; Leek et al. 2010). Study source, experimental protocol, time-point design, cellular context, and downstream transformation steps all leave structure in the gene by sample matrix. If such observations are treated as if they were homogeneous samples, the trained model may absorb technical structure together with biological structure.

Saluki and related work provided an important first-generation solution to this problem by integrating multiple studies into a compendium and extracting a consensus first principal component (PC1) label after unified preprocessing (Agarwal and Kelley 2022). That work provided both a mammalian half-life benchmark and a strong demonstration of sequence-based prediction. At the same time, it exposed a deeper issue. If the consensus label itself still carries study-level bias, then downstream performance comparisons reflect not only model capacity, but also label quality. For a benchmark-oriented computational biology study, this issue must be addressed explicitly.

We therefore place label auditing, rather than model scale, at the center of the paper and organize the study around three contributions. First, a complete leave-one-study-out reconstruction audits consensus-label stability: PC1 stability, which does not use Saluki labels, screens influential studies, while agreement with Saluki human PC1 and human-mouse ortholog concordance validate the consequences. Second, with Saluki human PC1 fixed as the target, we quantify both a human-only sequence/regulatory baseline and a cross-species transfer setting with external mouse ortholog priors, supported by permutation and residual controls. Third, we construct a weakly ortholog-regularized target at `λ = 0.10`, quantify its displacement, and use cross-target evaluation to test whether it preserves the predictable structure of the original human label.

The resulting evidence chain is deliberately bounded. Source heterogeneity is present; influential studies can be screened from label geometry and validated by two external criteria; the human-only model defines a transparent baseline; cross-species transfer improves only when gene-specific mouse-prior values are preserved; and weak target regularization changes the supervised target without reducing out-of-fold prediction of the original human label. The latter two settings alter input information and target definition, respectively, and are therefore reported separately.

## 2. Related Work

### 2.1 Large-scale measurement and integration of mammalian mRNA half-life

Early reviews and systematic studies established that variation in transcript degradation rate is a major source of gene-expression regulation in mammals (Ross 1995; Garneau et al. 2007; Schoenberg and Maquat 2012). Transcriptome-scale experiments subsequently accumulated across classical shutoff assays, metabolic labeling, and nucleoside-recoding strategies (Schwanhausser et al. 2011; Rabani et al. 2011; Tani et al. 2012; Schwalb et al. 2012; Duffy et al. 2015; Herzog et al. 2017; Schofield et al. 2018; Muhar et al. 2018). Duan et al. quantified interindividual RNA-stability variation in human lymphoblastoid cell lines; this source is represented as the `Gejman` study in the Saluki compendium (Duan et al. 2013). Friedel et al. highlighted conserved principles of mammalian half-life regulation (Friedel et al. 2009), and Agarwal and Kelley advanced the field to a systematic compendium and sequence-based benchmark (Agarwal and Kelley 2022). However, most existing studies treat the construction of a gene-level consensus label from a noisy multi-study matrix as an implicit preprocessing step rather than as an independent methodological problem.

### 2.2 Sequence and biochemical determinants of mRNA stability

Mechanistically, codon optimality, translation efficiency, and degradation coupling have become some of the most influential findings in the field (Presnyak et al. 2015; Wu et al. 2019; Hia et al. 2019). In parallel, microRNA (miRNA) targeting rules, 3'UTR cis-elements, and RNA structure have all been shown to affect half-life systematically (Mayr 2017; Grimson et al. 2007; Agarwal et al. 2015; Mauger et al. 2019). RNA modifications such as N6-methyladenosine (m6A) provide an additional explanatory layer (Meyer et al. 2012; Wang et al. 2014). Together, these results establish that sequence-to-stability prediction is not a purely black-box task, but a supervised learning problem with clear biological grounding.

### 2.3 Computational modeling from tabular learning to deep sequence models

At the modeling level, random forests, gradient boosting, and their modern implementations such as XGBoost and LightGBM remain strong baselines for heterogeneous feature sets (Breiman 2001; Friedman 2001; Chen and Guestrin 2016; Ke et al. 2017). Models closer to raw sequence input include DeepBind, cross-species regulatory sequence activity predictors, Enformer, and DNABERT (Alipanahi et al. 2015; Kelley 2020; Avsec et al. 2021; Ji et al. 2021). For the present study, the key question is not to show once again that deep models can be strong, but to quantify how label diagnosis and cross-species priors change predictive performance within a unified and transparent feature-learning framework.

## 3. Materials and Methods

The study contains three linked modules with distinct interpretation boundaries. The first is **label auditing**, which screens studies that alter consensus-label geometry and validates their effects using the Saluki processing convention and cross-species ortholog structure. The second is **fixed-target prediction**: Saluki human PC1 remains unchanged while human-only inputs are compared with a cross-species transfer setting that adds external mouse ortholog priors. The third is **weak ortholog regularization**, in which a mapped mouse signal is introduced into the human consensus label at `λ = 0.10` and tested for human dominance and preservation of original-human-label predictability. The second module changes input information; the third changes the target, so they are not interpreted on a shared leaderboard. Fig. 1 summarizes the workflow, and Table 1 lists the inputs, operations, and validation roles. Reproduction paths and commands are provided in Online Resource 1, Supplementary Section S15, and Data and Code Availability.

![](figures/main/Fig01_workflow.png)

**Fig. 1** Overall study-aware framework. a, Consensus-label reconstruction and auditing: `MOESM2` is processed using a minimum of three observed samples per gene, sample-wise z-scoring, iterative PCA imputation, quantile normalization, and PC1 extraction before leave-one-study-out analysis; `MOESM3` and Ensembl orthologs provide the Saluki reference and cross-species validation. b, Fixed-target prediction: Saluki human PC1 remains the target while human-only and mouse-prior transfer inputs are evaluated, with fold-wise permutation testing gene-prior identity. c, Weak ortholog regularization: 90% human no-Gejman PC1 and 10% mapped mouse PC1 define the target, which is checked by label distance, study-noise comparison, and cross-target evaluation. Values are generated from result tables by a reproducible vector-graphics script

**Table 1** Method modules, logical inputs, and validation roles. The main text retains only the analytical logic and intended role of each module. Specific filenames, script paths, result tables, and reproduction commands are listed in Supplementary Section S15 and in the versioned GitHub repository specified in Data and Code Availability.

| Method module | Logical input | Core operation | Validation role |
| --- | --- | --- | --- |
| Data organization | Public human/mouse half-life matrices and sample metadata | Parse species, study, method, cell type, and replicate; define human, mouse, and no-Gejman subsets | Clarify sample structure and gene universes, and prevent different task settings from being conflated |
| Consensus-label reconstruction | Gene by sample half-life matrix | Coverage filtering, sample standardization, missing-value imputation, distribution alignment, and PC1 extraction | Generate a reproducible half-life consensus label and align it to the released Saluki PC1 label |
| Study influence | Study-specific sample subsets | Rebuild the label after removing one study at a time; screen by PC1 stability and validate by Saluki agreement | Separate reference-free screening from reference-informed validation; downstream scores are not used |
| Ortholog validation | Human-mouse one-to-one ortholog genes | Compare human full PC1, no-Gejman PC1, Saluki human PC1, and mouse labels | Test whether label cleaning improves cross-species conserved structure |
| Human-only prediction | Human sequence/regulatory features | Gene-level out-of-fold (OOF) XGBoost benchmark | Quantify performance from human-only representation and define a transparent human-only baseline |
| Prior-enhanced prediction | Human features plus mouse ortholog priors | Prior-enhanced OOF benchmark, permutation control, and residual analysis | Test gene-specific cross-species signal in a transfer setting whose mouse-side covariates are external to the human target |
| 0.10 weak ortholog regularization | Human consensus label plus mapped mouse ortholog label | Mix at 90:10; quantify label shift, study-noise scale, and cross-target performance | Define a human-dominant mammalian stability target, interpreted separately from the fixed-target benchmark |

### 3.1 Data sources, metadata parsing, and common gene universes

We used the public supplementary data from the Saluki study by Agarwal and Kelley (Agarwal and Kelley 2022). `MOESM2` provides transformed human and mouse mRNA half-life gene by sample matrices and serves as the primary input for label reconstruction and leave-one-study-out auditing; the first two columns are Ensembl gene ID and gene name, and the remaining columns are transformed values from different studies, methods, cell types, and replicates. `MOESM3` provides the released `half-life (PC1)` summary label together with processed sample values; throughout this manuscript, the human side is denoted Saluki human PC1 and the mouse side is denoted Saluki mouse prior. The human worksheet in `MOESM3` indicates that the PC1 label was derived after removing `Gejman`, whereas the mouse worksheet uses the full mouse dataset. Accordingly, `MOESM3` is used here mainly as a reference and control for the Saluki processing convention rather than as a replacement training-label source. Precomputed Saluki sequence/regulatory feature blocks were obtained from the associated public dataset (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)). Because `MOESM2` already contains transformed values, we do not reapply a log transform and instead begin from sample-wise standardization.

Sample metadata were parsed from sample names and local curation tables, including `species`, `study`, `method`, `cell type`, and `replicate`. The full human analysis includes 54 samples, 19 studies, and 5 experimental methods. The mouse analysis includes 27 samples, 17 studies, and 3 experimental methods. The no-Gejman human subset includes 39 samples and 18 studies. Different analyses use explicitly defined gene universes: label auditing uses species-specific coverage-filtered genes; the global prior benchmark uses 12,916 common human genes; and the 0.10 ortholog-regularized target uses 12,307 shared genes, of which 10,768 have a one-to-one ortholog pair. Human-mouse ortholog mapping was obtained from Ensembl release 115 comparative genomics resources (Martin et al. 2023).

![](figures/main/Fig02_dataset_summary_cn.png)

**Fig. 2** Data sources, sample structure, and analysis gene universes. a, Numbers of samples, studies, and experimental methods in the human, mouse, and no-Gejman human subsets. b, Gene coverage retained after preprocessing. c, Gene or ortholog-pair universes used for MOESM3 controls, the global prior benchmark, ortholog validation, and 0.10 ortholog regularization. d, Roles of `MOESM2`, `MOESM3`, and Ensembl ortholog resources as primary reconstruction input, Saluki reference/control, and cross-species validation/mapping, respectively

### 3.2 Sample-level principal component analysis and diagnosis of compendium heterogeneity

Before constructing supervision labels, we examined the dominant structure of the compendium at the sample level using principal component analysis (PCA). For the raw transformed matrix, we applied iterative low-rank PCA imputation only to missing values for visualization, following the general low-rank treatment of incomplete high-dimensional matrices (Tipping and Bishop 1999; Troyanskaya et al. 2001). For the processed matrix, we applied the same preprocessing steps used in label reconstruction, namely sample-wise z-scoring, iterative PCA imputation, and quantile normalization (Bolstad et al. 2003). We then calculated the first two principal components in sample by gene space and colored the samples by study or method. The matrix definitions, imputation parameters, and interpretation boundary of the two PCA views are given in Supplementary Table S3.

The purpose of this step was not to generate the final label, but to test a key premise: whether the major variation in the matrix already contains pronounced study/method structure. If samples cluster primarily by study, downstream model performance cannot be interpreted only as sequence-to-stability signal; label provenance must also be audited. In the human compendium, preprocessing weakened but did not eliminate study/method structure, which motivated subsequent leave-one-study-out analysis. The mouse results are retained in Online Resource 1 to show that human and mouse heterogeneity are not fully symmetric.

![](figures/main/Fig03_human_pca_panel.png)

**Fig. 3** Human sample-level PCA diagnosis. a, Sample-level PCA of the raw transformed matrix. b, Sample-level PCA after sample-wise z-scoring, iterative PCA imputation, and quantile normalization. Each point represents one experimental sample; color denotes experimental method and marker shape denotes study. Both panels reveal clear study/method-related structure in the human compendium. Preprocessing weakens this structure but does not eliminate it, motivating further study-influence auditing

### 3.3 Consensus-label reconstruction

Consensus-label reconstruction was performed on the gene by sample matrix. For each species or study-removal subset, we applied coverage filtering, sample-wise z-scoring, iterative PCA imputation, and quantile normalization in sequence. We then extracted PC1 from the processed matrix and aligned its sign to the gene-wise mean. The resulting PC1 was used as the local consensus label.

This local procedure does not claim to reproduce every unpublished processing detail of the original Saluki workflow exactly. Instead, it defines an auditable and rerunnable baseline consensus-label reconstruction pipeline that remains sufficiently close to the released Saluki PC1 label. The reconstructed full human PC1 correlates with Saluki human PC1 at 0.948, and the reconstructed mouse PC1 correlates with the mouse reference at 0.958. For 0.10 weak ortholog regularization, we additionally use Saluki-like coverage thresholds of at least 10 non-missing samples for human and 5 for mouse. "Saluki-like" refers only to the coverage filter; both PC1 labels are reconstructed locally from `MOESM2`, not copied from the Saluki mouse prior.

### 3.4 Leave-one-study-out influence analysis

The study-influence analysis does not simply delete one column from a finished label matrix. Instead, for each study we reran the complete label-reconstruction pipeline after removing all samples from that study. For every candidate study, we reapplied coverage filtering, z-scoring, PCA imputation, quantile normalization, and PC1 extraction. We then calculated two metrics. `PC1 stability` is the Pearson or Spearman correlation between the leave-one-study-out label and the full PC1 label. `Saluki-PC1 agreement gain` is defined as `Pearson(leave-one-study-out PC1, Saluki human PC1) - Pearson(full PC1, Saluki human PC1)`.

To reduce circular selection through the Saluki reference, we interpret the audit in two stages. The first stage screens influential studies using `PC1 stability` alone; this metric does not use `MOESM3`, Saluki human PC1, or any downstream prediction score. The second stage examines `Saluki-PC1 agreement gain` as a reference-informed validation and then applies the human-mouse ortholog concordance analysis in Section 3.5 as an orthogonal cross-species validation. Because the Saluki data documentation already identifies `Gejman`, our contribution is reference-free recovery of its dominant influence and quantitative validation of the consequences of removal. Some mouse studies alter PC1 geometry, but none show a comparable convergence of evidence, so the same removal strategy is not applied.

### 3.5 Ortholog concordance as cross-species label validation

To avoid keeping label diagnosis entirely internal to the human compendium, we used human-mouse one-to-one orthologs as an external biological reference (Koonin 2005; Gabaldon and Koonin 2013). Specifically, human full PC1, human no-Gejman PC1, and Saluki human PC1 were mapped to one-to-one ortholog genes and compared with mouse PC1 or Saluki mouse prior by Pearson and Spearman correlation. The analysis was repeated both on all one-to-one ortholog pairs and on a high-confidence ortholog subset. Here, high confidence was not defined manually in this study; instead, we directly used the binary `is_high_confidence` field provided in Ensembl Compara release 115. All `ortholog_one2one` pairs were retained in the full analysis, and the subset with `is_high_confidence == 1` was analyzed separately.

This design provides a validation axis that is orthogonal to agreement with Saluki human PC1. If removing a human study improves agreement with the Saluki label but not human-mouse ortholog concordance, the result would be consistent with mechanical alignment to the Saluki processing convention. Concordant improvement instead supports recovery of a cross-species conserved mRNA-stability structure. The interpretation of `Gejman` therefore rests on the joint support of these two evidence types.

### 3.6 Construction and cross-target validation of the 0.10 ortholog-regularized target

This module changes the target, not the definition of the human features. It does not claim a new pure human half-life ground truth. Instead, it tests whether a small ortholog constraint can define a supervised target that remains human-dominant while representing conserved mammalian stability. We reconstructed human no-Gejman PC1 and mouse PC1 under Saluki-like coverage, z-scored both labels, and mapped the mouse label to human genes through one-to-one orthologs. For genes with an ortholog signal, we defined `T_0.10 = 0.90 × human z-scored label + 0.10 × mapped mouse z-scored label`; genes without an ortholog signal retained the human label. The resulting target was re-z-scored and sign-aligned to Saluki human PC1. We refer to it as the **0.10 ortholog-regularized target**.

The target was evaluated in three separate ways. First, within a common one-to-one ortholog universe, we calculated Pearson, Spearman, MSE, RMSE, and MAE against the human and mouse labels. Second, we compared the error in `T_0.10 - human label` with study-to-study variability among no-Gejman human studies on the same z-scored label scale. Third, we performed cross-target evaluation. On the same 12,307 genes, the same 1802 human-only `compact_all` features, the same 10 folds by 3 random seeds, and identical model parameters, models were trained on either human no-Gejman PC1 or `T_0.10` and then evaluated against the original human no-Gejman PC1. A gene-level paired bootstrap provided a 95% CI for the difference in Pearson correlation. This design directly tests whether regularization gains an apparent advantage by sacrificing predictability of the original human target.

`λ = 0.10` is the conservative operating point used for the primary analysis: the mouse contribution is limited to 10%, rather than selected to maximize prior-enhanced cross-validation. Sensitivity results for `λ = 0.05`, 0.10, and 0.30 are reported in full in the supplement. Larger values produce higher predictability after alignment with mouse priors, but also a larger target displacement, and are therefore not used for the primary claim.

### 3.7 Sequence/regulatory feature construction and evaluation of the human-only baseline

This section returns to the main benchmark route in which the target is unchanged and only human inputs are used. The human-only route uses two feature layers. The first layer is a 526-dimensional base sequence feature set derived from Ensembl release 115 representative transcripts, covering 5' untranslated region (5'UTR), coding sequence (CDS), and 3'UTR lengths, GC content, codon frequencies, region-level 3-mers, and 3'UTR 4-mers. The second layer is the 1802-dimensional `compact_all` feature set, formed by augmenting base sequence features with precomputed CWCS, SeqWeaver, and DeepRiPe regulatory blocks from the Saluki datapack (Alipanahi et al. 2015; Kelley 2020; Avsec et al. 2021).

All human-only benchmarks use Saluki human PC1 as the primary supervision target and are evaluated by gene-level out-of-fold prediction. The main results use XGBoost/CUDA. Tree-model hyperparameters were kept fixed in the main experiments to avoid over-tuning for a single figure. Typical settings included `tree_method=hist`, `device=cuda`, `max_depth=6`, `learning_rate=0.02`, `min_child_weight=4`, `subsample=0.9`, `colsample_bytree=0.75`, `reg_lambda=1.0`, and early stopping. The main text reports 5-fold OOF results, and a 10-fold by 3-random-seed evaluation is used as a split-sensitivity robustness check.

We did not attempt to recreate the exact fixed test split used by Saluki. All main results are gene-level OOF benchmarks. The purpose is to compare the relative effects of label processing, prior input, and control experiments within a unified feature table and a unified label universe. The human-only route is therefore presented as a transparent and reproducible baseline.

### 3.8 Cross-species prior-enhanced evaluation

This section is the mirror image of Section 3.6: the label is not changed, only the input is changed. In the global prior-enhanced benchmark, the supervision target remains Saluki human PC1 and the modeling unit remains each human gene. For each gene, we retained the original human `compact_all` features and concatenated six additional mouse-related columns: two prior values (`mouse_pc1` and Saluki mouse prior), two binary indicators for whether each prior is available, and two binary indicators for whether the ortholog mapping is high confidence. Here, `mouse_pc1` is our locally reconstructed mouse prior, whereas Saluki mouse prior is the released mouse-side prior from Saluki; the two priors carry related but non-identical cross-species stability information. Both priors are derived exclusively from mouse-side measurements and ortholog mapping. They are fixed before human cross-validation, and neither training-fold nor held-out human target values enter their construction. They are therefore external cross-species covariates at prediction time, and this setting is cross-species transfer rather than human-only prediction.

The role of these indicator columns is to tell the model separately about the prior value itself and about whether that prior exists and is high confidence. In implementation, if a human gene lacks a given mouse prior, the corresponding prior value is temporarily set to 0 and the associated availability indicator is set to 0. If the prior exists, the indicator is set to 1. The high-confidence columns are defined analogously and indicate whether the human-mouse ortholog belongs to the high-confidence Ensembl Compara subset, that is, `is_high_confidence == 1`. We used the full common-gene universe rather than only the prior-available subset so that performance could be reported on the complete benchmark universe rather than only on an easier prior-complete subset.

Thus, the training and evaluation set for the fixed-target prediction benchmark is the same `global prediction universe`, containing 12,916 human genes with Saluki human PC1 fixed as the supervised target. Human-only `compact_all`, single-prior models, dual-prior models, and shuffled-prior controls are all evaluated on these 12,916 genes by gene-level OOF prediction. The `both-priors` subset contains 11,107 genes with both `mouse_pc1` and the Saluki mouse prior available. It is used mainly for prior-coverage interpretation and residual analysis, not as the main training universe. The 0.10 ortholog-regularized analysis uses a separate 12,307-gene universe, and its one-to-one ortholog label-distance comparison uses 10,768 pairs. These sets overlap but are not identical.

The difference from the human-only benchmark is therefore clear: the supervision target is unchanged, but cross-species priors are added explicitly to the input. Section 3.7 answers what can be achieved using only human features, whereas Section 3.8 asks how much more can be gained for the same human target after adding ortholog priors.

### 3.9 Control experiments, evaluation metrics, and usage boundaries

We used four controls. First, permutation tests whether prior-enhanced gain depends on gene-specific prior values rather than simply on extra columns or missingness structure. Within each fold, `mouse_pc1` and Saluki mouse prior values were shuffled only among genes that originally had the respective value, while availability and high-confidence indicators were held fixed. Second, residual analysis tests whether the transfer result merely repeats the priors. It was restricted to the `both_priors_available` subset (N = 11,107), where prior-only OOF predictions were generated before human `compact_all` features were used to fit `Saluki human PC1 - prior prediction`. Third, no-direct-prior ablation tests whether prediction of the 0.10 ortholog-regularized target depends entirely on the locally reconstructed mouse PC1 used to construct that target. Fourth, the cross-target evaluation in Section 3.6 tests whether target regularization impairs OOF prediction of the original human label.

Model performance is reported by Pearson, Spearman, and coefficient of determination (R2). Label distance is reported by Pearson, Spearman, MSE, RMSE, MAE, median absolute error, and p95 absolute error. Unless otherwise noted, all prediction results are OOF rather than in-sample fits. Mean ± SD across random seeds describes split sensitivity and is not treated as a sampling confidence interval. Paired gene or ortholog bootstrap intervals use 2,000 resamples for the cross-target comparison and 5,000 resamples for the other key claims (Supplementary Tables S13-S14). Online Resource 1 retains the full rankings, multi-seed summaries, prior ablations, label-distance analyses, study-noise comparisons, and bootstrap outputs.

Finally, the prediction-usage boundary in this study has three levels. The base sequence model requires only explicitly partitioned 5'UTR, CDS, and 3'UTR sequence, and can therefore be applied to new transcripts or designed sequences, but it is not equivalent to the full `compact_all` feature setting. `compact_all` applies to existing Ensembl human genes because the regulatory blocks come from precomputed Saluki datapack features. The prior-enhanced model additionally requires mouse ortholog priors; if both mouse priors are missing, the prediction should fall back to the human-only compact model.

Generative AI tools were used for manuscript-structure and language-drafting support, code review, and figure-layout suggestions; they were not used to generate analytical data or final figures. All implementations, numerical results, references, figures, and final wording were verified by the authors, who take responsibility for the submitted work and its reproducibility.

**Table 2** Prediction settings, targets, inputs, and interpretation boundaries.

| Modeling setting | Target label | Input features | Intended interpretation and boundary |
| --- | --- | --- | --- |
| Base sequence | Saluki human PC1 | 526 sequence-statistics features | Applicable to new sequences, but only at the base-sequence level |
| Compact human-only | Saluki human PC1 | 1802-dimensional `compact_all` | Pure human sequence/regulatory benchmark and transparent human-only baseline |
| Global prior-enhanced | Saluki human PC1 | `compact_all` + two mouse priors + flags | Cross-species transfer using mouse-side covariates external to the human target |
| 0.10 weak ortholog regularization | 0.10 ortholog-regularized target | Human-only or prior-enhanced features | Human-dominant mammalian stability target; cross-target evaluation tests original-human-label predictability |
| Shuffled/no-direct controls | Same as corresponding parent target | Shuffled priors or removal of the direct prior value | Test whether the gain arises from gene-specific cross-species signal |

The Results are organized around three claims: study-aware label auditing; human-only and cross-species transfer under a fixed Saluki human PC1 target; and 0.10 weak ortholog regularization validated by label distance, study-noise comparison, and cross-target evaluation.

## 4. Results

### 4.1 Study-aware auditing identifies a dominant study influence supported by cross-species validation

The reconstructed full human label correlated with Saluki human PC1 at 0.948, indicating that the local pipeline reproduced its major structure. When studies were ranked by leave-one-study-out PC1 stability alone, removal of `Gejman` gave the lowest stability (r = 0.955; the next study was 0.973). Thus, `Gejman` was already the strongest influence in the screening stage, which used neither Saluki labels nor downstream prediction scores. Reference-informed validation then showed that removing it increased Pearson agreement with Saluki human PC1 from 0.948 to 0.983 on the respective available-gene sets and increased Spearman agreement from 0.950 to 0.986.

To avoid comparing correlations over different coverage sets, we repeated the analysis on 13,265 common genes with paired bootstrap resampling. Pearson agreement was 0.9515 for the full label and 0.9829 for the no-Gejman label, a paired gain of 0.0314 (95% CI, 0.0298-0.0330). The reference-free stability screen thus recovers the `Gejman` effect reported in the Saluki data documentation, while the common-gene analysis quantifies its magnitude and uncertainty. No human-only, prior-enhanced, or ortholog-regularized model score contributed to the study ranking.

![](figures/main/Fig04_human_study_influence.png)

**Fig. 4** Leave-one-study-out auditing recovers the dominant influence of `Gejman`. a, Points show `PC1 stability`, defined as the correlation between each leave-one-study-out PC1 and the full human PC1, for the eight studies with the largest effect on full-label geometry; the dashed line marks perfect agreement at r = 1. This screen does not use Saluki PC1. b, Reference-informed agreement change, defined as `Pearson(leave-one-study-out PC1, Saluki human PC1) - Pearson(full PC1, Saluki human PC1)`. `Gejman` has both the lowest stability and the only clearly positive agreement gain. The full ranking is described in Supplementary Section S4; its repository-relative result path is listed in Supplementary Table S15

To test whether the gain from removing `Gejman` extends beyond human-label agreement, we used human-mouse orthologs as an orthogonal cross-species validation. Across 12,592 one-to-one ortholog pairs, Pearson correlation between full human PC1 and mouse PC1 was 0.741 and increased to 0.755 after removing `Gejman`, a paired gain of 0.0147 (95% bootstrap CI, 0.0111-0.0181); Spearman correlation increased from 0.730 to 0.744. In the Ensembl high-confidence subset, the Pearson gain remained 0.0146 (95% CI, 0.0110-0.0181).

This result indicates that the influence of `Gejman` extends to the cross-species structure of the label. Ortholog concordance cannot by itself define a uniquely optimal human label, because biological and measurement differences between species remain. It nevertheless supplies evidence beyond agreement with Saluki human PC1: after removal, the human consensus label moves toward a more conserved cross-species stability axis rather than a direction dominated by one study.

![](figures/main/Fig05_ortholog_summary_cn.png)

**Fig. 5** Ortholog validation shows that removing `Gejman` improves cross-species concordance at multiple levels. a-b, Each point represents one one-to-one human-mouse ortholog pair; the x-axis is human consensus PC1 and the y-axis is mouse consensus PC1. a uses the full human PC1 and b uses the no-Gejman human PC1. The fitted lines and correlation values show that removing `Gejman` moves the point cloud closer to the shared principal axis. c, Pearson and Spearman comparisons across full PC1, no-Gejman PC1, and their respective high-confidence ortholog subsets show that the improvement is not tied to one correlation metric. d, Direct gains after removing `Gejman` in the full one-to-one and high-confidence ortholog sets

### 4.2 The 0.10 ortholog-regularized target remains human-dominant and passes cross-target evaluation

We first asked whether weak ortholog regularization with `λ = 0.10` still preserves dominant human signal. Before regularization, the reconstructed human no-Gejman PC1 and the locally reconstructed mouse PC1 under Saluki-like coverage correlated at 0.789 across 10,768 one-to-one orthologs, with RMSE = 0.663 and MAE = 0.512, indicating substantial but incomplete agreement.

Table 3 and Fig. 6 locate the regularized target directly. Within the same one-to-one ortholog universe, it nearly overlaps with human no-Gejman PC1, with Pearson = 0.9982 (95% bootstrap CI, 0.9981-0.9983), RMSE = 0.065, and MAE = 0.050. By contrast, its correlation with mouse PC1 is 0.827, with RMSE = 0.601 and MAE = 0.464. The primary effect of `λ = 0.10` is therefore a mild shift toward a conserved direction rather than replacement of the human label.

**Table 3** Same-scale comparison of the 0.10 ortholog-regularized target with human/mouse labels and study-level noise. All error metrics are computed in z-scored label units. Study-level rows report medians.

| Comparison | N or number of comparisons | Pearson | RMSE | MAE | Interpretation |
| --- | ---: | ---: | ---: | ---: | --- |
| Regularized target vs human no-Gejman PC1 | 10,768 genes | 0.998 | 0.065 | 0.050 | The dominant human label structure is preserved |
| Regularized target vs mouse PC1 | 10,768 genes | 0.827 | 0.601 | 0.464 | Mild shift toward a cross-species conserved direction |
| Human no-Gejman PC1 vs mouse PC1 | 10,768 genes | 0.789 | 0.663 | 0.512 | Human and mouse labels already differ before regularization |
| No-Gejman study-study label pairs (median) | 153 study pairs | 0.530 | 0.954 | 0.720 | Noise scale among ordinary studies |
| Single-study vs human no-Gejman PC1 (median) | 18 studies | 0.742 | 0.715 | 0.543 | Deviation of a single study from the human consensus |

![](figures/main/Fig06_ortholog_regularized_target_cn.png)

**Fig. 6** The 0.10 ortholog-regularized target retains human-dominant structure and shows no detectable reduction in original-human-label predictability. a, The regularized target lies almost exactly on the diagonal relative to human no-Gejman PC1. b, It remains correlated with mouse PC1 but is not a mouse-like label. c, The target shift is substantially smaller than study-study and study-consensus differences on the same z-scored scale. d, Models trained on the human or regularized target use the same 12,307 genes, human-only features, and 10 folds by 3 seeds, and are both evaluated against the original human target; the paired Pearson difference is +0.0003 with a 95% CI spanning zero

This shift is also small relative to ordinary study-level variability. After excluding `Gejman`, pairwise comparisons among study mean labels gave median Pearson = 0.530, RMSE = 0.954, and MAE = 0.720. Single-study labels relative to human no-Gejman PC1 gave median Pearson = 0.742, RMSE = 0.715, and MAE = 0.543. By comparison, the regularized-target shift had RMSE = 0.065 and MAE = 0.050. Fig. 6c places these quantities on the same z-scored scale; full calculations are given in Supplementary Table S12.

More importantly, cross-target evaluation showed no detectable trade-off in human-label predictability. Against the original human no-Gejman PC1, direct training on the human target achieved r = 0.7476 ± 0.0010, whereas training on the regularized target and returning to the same human evaluation achieved r = 0.7480 ± 0.0010. The paired difference after combining the three OOF predictions was +0.0003 (95% CI, -0.0006 to 0.0012). At the present data resolution, `λ = 0.10` therefore preserves the predictable structure of the original human label. We treat it as a distinct human-dominant mammalian stability target, not as a higher-ranked replacement for the fixed Saluki human PC1 benchmark.

### 4.3 A pure human compact sequence/regulatory model defines a stable baseline

Along the pure human route, the primary global fixed-target analysis gave Pearson = 0.745, Spearman = 0.737, and R2 = 0.554 for `compact_all` across 12,916 common genes. As an independent matched comparison, Fig. 7a evaluates Ridge, HGB, GPU XGBoost, and `compact_all` under one common-universe workflow, where the highest value is 0.743; the small difference from 0.745 reflects cross-validation and training configuration rather than model blending. Fig. 7b is a feature-block comparison rather than a sequential ablation. CWCS and DeepRiPe improve over the base setting when added individually, SeqWeaver alone is nearly unchanged, and the combined `compact_all` setting is highest. The evidence therefore supports complementary feature integration, not monotonic gain from every block.

Human-only performance also depends on label definition. Fig. 7c shows that within the same common-gene universe, the Pearson correlation of `compact_all` with `human_full_pc1` is 0.721, compared with 0.737 for human no-Gejman PC1 and 0.743 for Saluki human PC1. Thus, the no-Gejman label has greater downstream predictability under the matched feature setting. Fig. 7d further shows an orderly pattern across all-common, ortholog, and higher-signal subsets, without concentrating the conclusion in one unusually easy subset.

![](figures/main/Fig07_sequence_progression_cn.png)

**Fig. 7** Four evidence panels for the human-only benchmark. a, Matched algorithm/feature comparison across Ridge, HGB, GPU XGBoost, and `compact_all`. b, Regulatory-block comparison; individual-block contributions are heterogeneous, while the combined `compact_all` setting is highest. c, Label-definition sensitivity across full human PC1, human no-Gejman PC1, and Saluki human PC1 on the same common-gene universe. d, Gene-subset sensitivity across all-common, all-ortholog, high-confidence ortholog, and top-regulatory-33% subsets; subset values are not extrapolated as global performance

### 4.4 Cross-species transfer produces a gene-specific prior-enhanced gain

When the reconstructed `mouse PC1` and Saluki mouse prior were added to `compact_all`, performance increased sharply. This is a cross-species transfer setting: the mouse-side covariates are external to the human target, which remains Saluki human PC1. Fig. 8a compares four global settings directly. Human-only `compact_all` reached Pearson = 0.745. Adding only the reconstructed `mouse PC1` increased Pearson to 0.828. Adding only the Saluki mouse prior gave 0.827. Adding both priors together increased Pearson to 0.829. Thus, each single prior already provides a large gain, and the dual-prior model yields the best single-model result, indicating that the two mouse priors contain related but not fully redundant cross-species stability information.

To exclude dependence on a single favorable split, we performed repeated 10-fold gene-level OOF evaluation with three random states. Human-only `compact_all` achieved Pearson = 0.748 ± 0.001, the shuffled-prior control achieved 0.748 ± 0.001, the real dual-prior model achieved 0.830 ± 0.001, and the coverage-aware fallback achieved 0.832 ± 0.001. Here, ± denotes the SD across random states, not a confidence interval. The results show that the prior-enhanced gain is stable across splits and that genes without priors should explicitly fall back to the human-only predictor rather than be interpreted under a prior-complete input condition.

Prior coverage itself was also non-uniform. Fig. 8c shows that among 12,916 common human genes, 11,107 had both mouse priors available, 462 had only the reconstructed `mouse PC1`, and 1,347 lacked both priors. No gene had only the Saluki mouse prior without the reconstructed `mouse PC1`. This distribution defines the operational boundary of the model. Prior-enhanced results belong to the cross-species transfer setting, whereas missing-prior genes are handled by coverage-aware fallback rather than interpreted under the same input condition as prior-complete genes.

![](figures/main/Fig08_prior_summary_cn.png)

**Fig. 8** Four evidence panels for the prior-enhanced benchmark. a, Global fixed-target benchmark across human-only, reconstructed-mouse-PC1, Saluki-mouse-prior, and dual-prior inputs. b, Repeated 10-fold controls; error bars are across-seed SD, and shuffling retains prior availability while disrupting gene identity mapping. The gain disappears when gene-prior identity is disrupted. c, Prior availability in the common human gene universe, which defines the applicability boundary. d, Ortholog-conditioned gains in the full one-to-one and high-confidence subsets

Permutation control directly tests whether prior-enhanced gain depends on gene-specific prior values. The real `compact_all_plus_both_mouse_priors_global` model reached Pearson = 0.829, whereas the fold-wise prior-shuffled version `compact_all_plus_both_mouse_priors_shuffled_global` dropped back to 0.745. R2 likewise fell from 0.688 to 0.553, nearly matching the human-only `compact_all_global` baseline. Prior availability and high-confidence indicators were held fixed during shuffling. The disappearance of the gain therefore reflects the broken mapping between prior value and gene identity rather than a simple feature-count or missingness effect.

The repeated 10-fold results show the same pattern. After averaging the three OOF predictions for each gene, the seed-averaged real-prior prediction achieved r = 0.8324 and the seed-averaged shuffled-prior prediction achieved r = 0.7531, a paired gain of 0.0794 (95% CI, 0.0735-0.0855). The gain over the seed-averaged human-only prediction was 0.0799 (95% CI, 0.0740-0.0862). This averaging summarizes split sensitivity and is not an additional model-blending step. Fig. 8d further shows the same direction in both the full one-to-one ortholog subset and the high-confidence subset. The gain therefore depends on **gene-specific cross-species prior information**. Disrupting gene identity mapping removes it. This is controlled transfer from mouse ortholog labels, not performance available from human sequence alone.

Residual decomposition further clarifies that the prior is strong, but not sufficient by itself. This analysis was performed only on the `both_priors_available` subset (N = 11,107) to decompose prior signal, not to replace the main 12,916-gene benchmark. Prior-only RidgeCV achieved Pearson = 0.792 and R2 = 0.627, whereas human `compact_all` alone achieved Pearson = 0.739 and R2 = 0.544. If the two mouse priors are first used to produce a prior prediction and human `compact_all` is then used only to model the residual `Saluki human PC1 - prior prediction`, the final result reaches Pearson = 0.826 and R2 = 0.675. Using `(R2_final - R2_prior) / (1 - R2_prior)`, human features explain about 12.8% of the variance left unexplained by the prior. The two-stage residual model is therefore best interpreted as a combination of mouse-prior signal and human-feature correction.

### 4.5 Weak ortholog regularization improves conserved-target predictability without changing the fixed-target ranking

After establishing human dominance and cross-target preservation, we evaluated the 0.10 ortholog-regularized target under different input conditions. It remained strongly correlated with Saluki human PC1 (r = 0.987). On the same 12,307-gene universe, human-only prediction achieved 0.754 ± 0.001; adding both mouse priors increased prediction of the regularized target to 0.855 ± 0.001, whereas the prior-enhanced baseline against fixed Saluki human PC1 was 0.839 ± 0.001. These values do not form a common target ranking. They show that weak regularization aligns the target more closely with a conserved component captured by explicit ortholog covariates.

Prior ablation supports this interpretation. Removing the locally reconstructed mouse PC1 used to construct the target, while retaining the Saluki mouse prior and indicators, still gave 0.852 ± 0.001. Retaining only missingness and high-confidence indicators returned performance to 0.753 ± 0.001, close to the human-only level. Together with Fig. 6d, these results make `λ = 0.10` a bounded core label-construction result: it defines a human-dominant mammalian stability target and improves conserved-target predictability when external ortholog information is available, but it does not replace the Saluki human PC1 or pure-sequence benchmark ranking.

## 5. Discussion

This study separates three issues that are often conflated in mRNA half-life benchmarks: how to audit a multi-study consensus label; what can be claimed from human-only and cross-species-transfer inputs under a fixed human target; and whether a weak ortholog constraint can define a mammalian stability target that remains human-dominant. The three claims are supported by stability and ortholog validation, permutation and residual controls, and label-distance plus cross-target evaluation, respectively, rather than by a single maximum score.

The `Gejman` result is supported at two evidence levels. PC1 stability identifies it as the most influential study without using the Saluki reference; agreement gain is reference-informed validation, and ortholog concordance supplies an orthogonal cross-species validation. This design establishes dominant influence on the present consensus-label geometry without implying that `Gejman` should be removed for every biological question. The human-only `compact_all` model then defines a fixed-target baseline, while permutation, paired bootstrap, and residual analysis show that prior-enhanced gain depends on gene-specific cross-species signal. `MOESM3`-derived labels serve as reference and control rather than replacement supervision targets (Supplementary Tables S6-S7).

The third core contribution is the `λ = 0.10` ortholog-regularized target. Geometry shows that it remains nearly identical to human no-Gejman PC1, study-noise comparison places its displacement below ordinary inter-study variation, and cross-target evaluation shows no detectable loss of original-human-label predictability. Its higher prior-enhanced score reflects increased predictability of a conserved component from explicit ortholog covariates. The resulting quantity is an auditable, human-dominant mammalian stability score rather than a direct experimental half-life measurement or a pure human sequence benchmark.

## 6. Limitations

The study has three limitations. First, label auditing and the principal benchmarks use one public compendium and require validation in a second independent half-life resource; the current analysis also treats study influence rather than a hierarchy of samples, cell types, and studies. Leave-one-study-out influence combines the number of samples removed with directional differences in their profiles. Without a size-matched null analysis, `Gejman` is therefore interpreted as the dominant influence in this compendium rather than a sample-count-adjusted statistical outlier. Second, cross-species transfer and 0.10 ortholog regularization alter input information and target definition, respectively, and cannot be ranked jointly with the pure human fixed-target task. The value `λ = 0.10` is a conservative empirical weight whose portability and biological scale require independent testing. Third, we did not retrain Saluki or other large sequence models, and the implementation is not a one-click predictor for arbitrary unannotated sequences because the base-sequence, compact-regulatory, and prior-enhanced levels require different inputs.

## 7. Conclusion

This study-aware reanalysis shows that a reference-free PC1-stability screen recovers the dominant influence of `Gejman` on the human consensus label, while Saluki agreement and ortholog concordance quantify and validate that influence. Under a fixed Saluki human PC1 target, a compact human-only model provides a transparent baseline; mouse priors raise repeated 10-fold cross-species-transfer performance to Pearson = 0.830, and disrupting gene-prior identity removes the gain. Weak ortholog regularization at `λ = 0.10` defines a target that remains nearly identical to human no-Gejman PC1 and preserves original-human-label predictability in cross-target evaluation. The framework therefore separates label auditing, fixed-target transfer, and weak ortholog regularization into reproducible and independently testable claims.

## 8. Statements and Declarations

### Funding

This work was supported by the National Natural Science Foundation of China (22373101, 22203089, and 22573043).

### Author Contributions

Xu Jin performed data curation, computational analysis, reproducibility checks, figure preparation, and manuscript drafting. Dinglin Zhang and Yingchen Mao contributed to study conception, method design, result interpretation, project supervision, and manuscript revision. Wenzhuo Wang, Anhui Wang, and Yuebin Zhang contributed to result discussion and manuscript review. All authors reviewed and approved the final manuscript.

### Data Availability

The half-life matrices analyzed here are from the public supplementary material of the Saluki study; the precomputed Saluki feature datapack is available from Zenodo (DOI: [10.5281/zenodo.6326409](https://doi.org/10.5281/zenodo.6326409)); and ortholog mappings are from Ensembl comparative genomics resources (Agarwal and Kelley 2022; Martin et al. 2023). Result tables, figure source data, and reproduction notes are publicly available in the versioned repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label).

### Code Availability

Code, figure-generation scripts, and result-reproduction workflows are publicly available in [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label); the manuscript-reviewed state is identified by the release tag `mRNA-PC1-label-v1.0`.

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

Friedel CC, Dolken L, Ruzsics Z, Koszinowski UH, Zimmer R (2009) Conserved principles of mammalian transcriptional regulation revealed by RNA half-life. Nucleic Acids Research 37(17):e115. https://doi.org/10.1093/nar/gkp542

Friedman JH (2001) Greedy function approximation: A gradient boosting machine. The Annals of Statistics 29(5):1189-1232. https://doi.org/10.1214/aos/1013203451

Gabaldon T, Koonin EV (2013) Functional and evolutionary implications of gene orthology. Nature Reviews Genetics 14(5):360-366. https://doi.org/10.1038/nrg3456

Garneau NL, Wilusz J, Wilusz CJ (2007) The highways and byways of mRNA decay. Nature Reviews Molecular Cell Biology 8:113-126. https://doi.org/10.1038/nrm2104

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

Martin FJ, Amode MR, Aneja A, Austine-Orimoloye O, Azov AG, Barnes I, Becker A, Bennett R, Berry A, Bhai J, et al (2023) Ensembl 2023. Nucleic Acids Research 51(D1):D933-D941. https://doi.org/10.1093/nar/gkac958

Mauger DM, Cabral BJ, Presnyak V, Su SV, Reid DW, Goodman B, Link K, Khatwani N, Reynders J, Moore MJ, McFadyen IJ (2019) mRNA structure regulates protein expression through changes in functional half-life. Proceedings of the National Academy of Sciences 116(48):24075-24083. https://doi.org/10.1073/pnas.1908052116

Mayr C (2017) Regulation by 3'-untranslated regions. Annual Review of Genetics 51:171-194. https://doi.org/10.1146/annurev-genet-120116-024704

Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. Cell 149(7):1635-1646. https://doi.org/10.1016/j.cell.2012.05.003

Muhar M, Ebert A, Neumann T, Umkehrer C, Jude J, Wieshofer C, Rescheneder P, Lipp JJ, Herzog VA, Reichholf B, et al (2018) SLAM-seq defines direct gene-regulatory functions of the BRD4-MYC axis. Science 360(6390):800-805. https://doi.org/10.1126/science.aao2793

Musaev D, Abdelmessih M, Vejnar CE, Yartseva V, Weiss LA, Strayer EC, Takacs CM, Giraldez AJ (2024) UPF1 regulates mRNA stability by sensing poorly translated coding sequences. Cell Reports 43:114074. https://doi.org/10.1016/j.celrep.2024.114074

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
