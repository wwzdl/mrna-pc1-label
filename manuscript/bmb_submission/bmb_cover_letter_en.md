# Cover Letter

Dear Editor,

We are pleased to submit our manuscript entitled **"Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction"** for consideration as an Original Research Article in *Bulletin of Mathematical Biology*.

The manuscript addresses a benchmark problem in mammalian mRNA half-life prediction: widely used supervision labels are consensus scores constructed from heterogeneous multi-study measurements. We therefore separate three estimands: study-aware label auditing, fixed-human-target prediction with or without external mouse ortholog covariates, and ortholog-informed shrinkage of the target itself.

The study provides three main results. First, a Saluki-label-independent PC1-stability screen identifies `Gejman` as the largest influence in the ordinary sample-weighted analysis. A new size-matched null shows that removing 15 samples explains much of the geometric displacement, but does not reproduce either recovery of the published Saluki processing choice or the positive human-mouse ortholog-concordance change (one-sided empirical p = 0.002 for each comparison). Second, with Saluki human PC1 fixed, external mouse ortholog priors improve repeated 10-fold prediction from Pearson = 0.748 to 0.830; the gain disappears when fold-wise permutation disrupts gene-prior identity. Third, a target defined by 0.10 ortholog-informed shrinkage remains nearly identical to human no-Gejman PC1, and cross-target evaluation detects no loss of original-human-label predictability. This target-construction analysis is kept separate from the fixed-target ranking.

The work fits *Bulletin of Mathematical Biology* by treating consensus-label geometry, study influence, cross-species covariates, and target shrinkage as explicit quantitative objects with distinct validation rules. It offers an auditable framework rather than a model-score comparison alone.

This manuscript is original, has not been published previously, and is not under consideration elsewhere. All authors have reviewed the manuscript and approved its submission. The authors declare no known competing financial interests or personal relationships that could have influenced the work reported in this manuscript. Code, result tables, figure-generation scripts, and reproduction notes are publicly available in the versioned repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label).

Thank you for considering our manuscript.

Sincerely,

Yingchen Mao and Dinglin Zhang  
Corresponding authors
