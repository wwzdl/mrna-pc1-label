# 给编辑的附信（中文内部版）

尊敬的编辑：

我们谨提交题为 **“Study-aware label auditing, cross-species transfer, and ortholog-informed target shrinkage for mammalian mRNA half-life prediction”** 的稿件，申请作为 *Bulletin of Mathematical Biology* 的 `Original research article` 审阅。

本文关注哺乳动物 mRNA 半衰期预测中的一个基础 benchmark 问题。常用监督标签是由异质 multi-study measurements 构建的共识分数。我们审计 study 组成如何影响这些标签，评估加入外部 mouse ortholog covariates 后的固定目标预测，并考察弱同源基因信息目标收缩。

Leave-one-study-out PC1 stability 在普通 sample-weighted 分析中将 `Gejman` 识别为最大影响来源。同样本量零分布表明，删除 15 个样本解释了较大部分几何位移，但随机删除不能复制对 Saluki 已发布处理选择的恢复或 human-mouse ortholog concordance 的正向变化（两项比较的单侧经验均为 `p=0.002`）。在固定 Saluki human PC1 时，外部 mouse ortholog priors 将重复 10-fold Pearson 从 0.748 提高到 0.830；fold 内打乱 gene-prior mapping 后增益消失。`λ = 0.10` 同源基因目标收缩仍与 human no-Gejman PC1 几乎一致，cross-target evaluation 未检测到原 human 标签可预测性下降。该目标分析与固定目标比较分开报告。

本文建立了审计 study influence、评估 cross-species priors 和量化 target shrinkage 所引入变化的定量方法。其重点是透明的统计验证和可复现的生物学预测，与 *Bulletin of Mathematical Biology* 的定位契合。

本稿件为原创工作，未曾发表，也未同时投往其他期刊。两位作者均已参与当前稿件内容的整理与审阅，并将在正式投稿前完成最终确认。作者为 Wenzhuo Wang 和 Ying Shao；通讯作者为 Ying Shao。本研究未获得专项资助。作者声明不存在任何已知的、可能对本文所报告工作产生影响的竞争性经济利益或个人关系。代码、结果表、图件生成脚本和复现说明已通过版本化仓库 [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label) 公开。两位作者的 ORCID 尚未提供。

感谢您的审阅与考虑。

此致  
敬礼

Ying Shao  
通讯作者
