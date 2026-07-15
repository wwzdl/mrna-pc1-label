# 给编辑的附信（中文内部版）

尊敬的编辑：

我们谨提交题为 **“Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction”** 的稿件，申请作为 *Bulletin of Mathematical Biology* 的 `Original research article` 审阅。

本文关注哺乳动物 mRNA 半衰期预测中的一个基础 benchmark 问题：常用监督标签是由异质 multi-study measurements 构建的共识分数。为避免混淆，我们分开处理三个 estimands：study-aware 标签审计、固定 human target 下有无外部 mouse ortholog covariates 的预测，以及对监督目标本身进行同源信息收缩。

本文得到三项主要结果。第一，不依赖 Saluki 标签的 PC1 stability screen 在普通 sample-weighted 分析中将 `Gejman` 识别为最大影响来源。新增的同样本量零分布表明，删除 15 个样本解释了较大部分几何位移，但无法复制对 Saluki 已发布处理选择的恢复，也无法复制 human-mouse ortholog concordance 的正向变化（两项比较的单侧经验均为 `p=0.002`）。第二，在固定 Saluki human PC1 时，外部 mouse ortholog priors 将重复 10-fold Pearson 从 0.748 提高到 0.830；fold 内打乱 gene-prior identity 后增益消失。第三，`λ = 0.10` 同源基因目标收缩仍与 human no-Gejman PC1 几乎一致，cross-target evaluation 未检测到原 human 标签可预测性下降。该目标构建结果与固定目标排行榜始终分开。

我们认为该工作与 *Bulletin of Mathematical Biology* 的定位契合，因为它把 consensus-label geometry、study influence、cross-species covariates 和 target shrinkage 作为具有不同验证规则的定量对象，而不是只报告一次模型得分比较。

本稿件为原创工作，未曾发表，也未同时投往其他期刊。所有作者均已参与当前稿件内容的整理与审阅；正式投稿前，将由通讯作者在系统中完成最终作者确认。作者为 Xu Jin、Wenzhuo Wang、Anhui Wang、Yuebin Zhang、Yingchen Mao 和 Dinglin Zhang；共同通讯作者为 Yingchen Mao（`myc@lnnu.edu.cn`）和 Dinglin Zhang（`dlzhang@dicp.ac.cn`）。本工作受国家自然科学基金（22373101、22203089、22573043）资助。作者声明不存在任何已知的、可能对本文所报告工作产生影响的竞争性经济利益或个人关系。代码、结果表、图件生成脚本和复现说明已通过版本化仓库 [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label) 公开。当前已公开确认的 ORCID 包括 Xu Jin（0009-0007-7776-0541）、Anhui Wang（0000-0003-4041-9095）、Yingchen Mao（0000-0002-2469-054X）和 Dinglin Zhang（0000-0002-9674-0577）；Wenzhuo Wang 和 Yuebin Zhang 的 ORCID 仍待作者本人最后核对。

感谢您的审阅与考虑。

此致  
敬礼

作者团队
