# Cover Letter

Dear Editor,

We are pleased to submit our manuscript entitled **"Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction"** for consideration as an Original Research Article in *Bulletin of Mathematical Biology*.

The manuscript addresses a benchmark problem in mammalian mRNA half-life prediction. Widely used supervision labels are consensus scores constructed from heterogeneous multi-study measurements. We audit how study composition shapes these labels, evaluate fixed-target prediction with external mouse ortholog covariates, and examine weak ortholog-informed target shrinkage.

Leave-one-study-out PC1 stability identifies Gejman as the largest influence in the ordinary sample-weighted analysis. A conditional same-size comparison, in which 15 non-Gejman samples are deleted while Gejman is retained, shows that sample count explains much of the geometric displacement. These comparator deletions do not reproduce recovery of the published Saluki processing choice or the positive human-mouse ortholog change. With Saluki human PC1 fixed, external mouse ortholog priors improve repeated 10-fold prediction from Pearson = 0.748 to 0.830. The gain disappears when partition-wise permutation within each outer fold disrupts the gene-prior mapping. Finally, 0.10 ortholog-informed target shrinkage remains nearly identical to human no-Gejman PC1; cross-target evaluation detects no reduction in original-human-label predictability but is not a formal equivalence test. A strict ablation attributes only 0.0013 of target-prediction correlation to direct reuse of the exact mouse vector contributing 10% of the target. We report this target analysis separately from the fixed-target comparison.

The manuscript develops quantitative procedures to audit study influence, assess cross-species priors, and measure changes introduced by target shrinkage. This emphasis on transparent statistical validation and reproducible biological prediction fits the scope of *Bulletin of Mathematical Biology*.

This manuscript is original, has not been published previously, and is not under consideration elsewhere. All authors have reviewed the manuscript and approved its submission. The authors declare no known competing financial interests or personal relationships that could have influenced the work reported in this manuscript. Code, result tables, figure-generation scripts, and reproduction notes are publicly available in the versioned repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label).

Thank you for considering our manuscript.

Sincerely,

Ying Shao  
Corresponding author<br>
Email: yshao@dlmu.edu.cn
