# Cover Letter

Dear Editor,

We are pleased to submit our manuscript entitled **"Study-aware label auditing, weak ortholog regularization, and cross-species transfer for mammalian mRNA half-life prediction"** for consideration as an Original Research Article in *Bulletin of Mathematical Biology*.

This manuscript addresses a benchmark problem in mammalian mRNA half-life prediction: commonly used supervision labels are processed consensus scores derived from heterogeneous multi-study compendia rather than direct experimental ground truth. We first audit label stability, then separate a fixed-target human-only benchmark from cross-species transfer using mouse-side ortholog priors external to the human target, and finally test whether weak ortholog regularization at λ = 0.10 can introduce conserved mammalian stability signal while retaining a human-dominant target.

The main contributions are threefold. First, a reference-free PC1-stability ranking recovers the dominant influence of `Gejman` without using Saluki labels or downstream model scores; Saluki agreement, ortholog concordance, and paired bootstrap analysis then quantify its consequences. Second, with Saluki human PC1 held fixed, compact human sequence/regulatory features define a transparent human-only baseline, whereas mouse priors enter only as mouse-side covariates in a cross-species transfer setting. No human target value is used to construct these priors, and the gain disappears when fold-wise permutation disrupts gene-prior identity. Third, the weakly ortholog-regularized target remains nearly identical to human no-Gejman PC1 and preserves original-human-label predictability in cross-target evaluation. This core label-construction result is evaluated separately from the fixed-target performance ranking.

We believe the manuscript fits the scope of *Bulletin of Mathematical Biology* because it treats a biologically meaningful phenotype as a quantitative object whose labels, study structure, cross-species conservation, and predictive boundaries can be audited mathematically. The work is not simply a report of a higher score; it provides a reproducible framework for separating label diagnosis, fixed-target baseline modeling, and controlled cross-species prior enhancement in mRNA stability prediction.

This manuscript is original, has not been published previously, and is not under consideration elsewhere. All authors have reviewed the manuscript and approved its submission. The authors declare no known competing financial interests or personal relationships that could have influenced the work reported in this manuscript. Code, result tables, figure-generation scripts, and reproduction notes are publicly available in the versioned repository [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label).

Thank you for considering our manuscript.

Sincerely,

Yingchen Mao and Dinglin Zhang  
Corresponding authors
