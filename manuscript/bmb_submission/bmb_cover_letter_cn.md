# 给编辑的附信（中文内部版）

尊敬的编辑：

我们谨提交题为 **“Study-aware label auditing, weak ortholog regularization, and cross-species transfer for mammalian mRNA half-life prediction”** 的稿件，申请作为 *Bulletin of Mathematical Biology* 的 `Original research article` 审阅。

本文关注哺乳动物 mRNA 半衰期预测中的一个基础 benchmark 问题：常用监督标签并不是天然干净的实验真值，而是由异质 multi-study compendium 处理得到的共识分数。我们先审计标签稳定性，再在固定 Saluki human PC1 目标下分开评估 human-only benchmark 与使用 human target 外部 mouse-side ortholog priors 的 cross-species transfer，最后检验 `λ = 0.10` 的弱同源基因正则化能否在保持 human 主导性的同时引入保守的哺乳动物稳定性信号。

本文的主要贡献可概括为三点：

1. Reference-free PC1 stability 排序在不使用 Saluki 标签或下游模型分数的情况下重现了 `Gejman` 对 human 共识标签的主导影响；Saluki agreement、ortholog concordance 和 paired bootstrap 随后量化并验证了其后果。
2. 在 Saluki human PC1 固定不变时，compact human sequence/regulatory 模型给出透明的 human-only 基线；mouse priors 仅作为外部协变量进入 cross-species transfer。Fold 内打乱 gene-prior 对应关系后增益消失，说明提升依赖 gene-specific cross-species signal。
3. `λ = 0.10` 弱同源正则化产生与 human no-Gejman PC1 高度一致的 human-dominant target。标签距离、study-noise comparison、cross-target evaluation 和 prior ablation 共同界定该核心标签构建结果；其性能与固定 Saluki human PC1 benchmark 分开报告。

我们认为该工作与 *Bulletin of Mathematical Biology* 的定位契合，原因在于它并不只是报告一个更高分数，而是围绕一个具有明确生物学意义的表型，系统讨论了标签构建、异质性诊断、跨物种保守性验证以及受约束的预测建模边界。本文的核心结果涉及：公开哺乳动物转录后调控数据在何种条件下可以被整合为可审计的数学对象；以及在清楚区分输入增强与目标改变之后，mRNA stability prediction 能得出哪些可信结论。我们相信，这种把 biological phenotype、ortholog conservation 和 transparent benchmark design 结合起来的研究，会对期刊读者中的计算生物学、系统生物学和定量建模研究者具有吸引力。

本稿件为原创工作，未曾发表，也未同时投往其他期刊。所有作者均已参与当前稿件内容的整理与审阅；正式投稿前，将由通讯作者在系统中完成最终作者确认。作者为 Xu Jin、Wenzhuo Wang、Anhui Wang、Yuebin Zhang、Yingchen Mao 和 Dinglin Zhang；共同通讯作者为 Yingchen Mao（`myc@lnnu.edu.cn`）和 Dinglin Zhang（`dlzhang@dicp.ac.cn`）。本工作受国家自然科学基金（22373101、22203089、22573043）资助。作者声明不存在任何已知的、可能对本文所报告工作产生影响的竞争性经济利益或个人关系。代码、结果表、图件生成脚本和复现说明已通过版本化仓库 [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label) 公开。当前已公开确认的 ORCID 包括 Xu Jin（0009-0007-7776-0541）、Anhui Wang（0000-0003-4041-9095）、Yingchen Mao（0000-0002-2469-054X）和 Dinglin Zhang（0000-0002-9674-0577）；Wenzhuo Wang 和 Yuebin Zhang 的 ORCID 仍待作者本人最后核对。

感谢您的审阅与考虑。

此致  
敬礼

作者团队
