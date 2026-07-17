# BMB 投稿要点摘要（中文内部版，非系统必传）

本文件不是 BMB 明确要求上传的正式材料，而是给后续英文 `cover letter`、投稿系统备注和回复编辑时统一口径用的内部摘要。

## 1. 一句话定位

这不是“又一个分数更高的 mRNA half-life 预测模型”，而是一篇把 **标签诊断、任务边界和跨物种先验增强** 明确拆开的 benchmark 论文。

## 2. 编辑最该快速看到的三项主贡献和两个解释边界

1. **标签先于模型。** 本文首先问的是公开 half-life compendium 的监督标签是否可靠，而不是直接堆更大的预测器。
2. **高影响 study 先由 Saluki-label-independent 稳定性筛查。** Leave-one-study-out PC1 stability 在不使用 Saluki 标签或下游模型分数时将 `Gejman` 排为 primary sample-weighted analysis 的第一位；size-matched null 与 study balancing 限定几何解释，Saluki agreement 与 ortholog concordance 随后检查改变方向。
3. **human-only 与 prior-enhanced 明确分账。** Mouse-side priors 独立于 human target 构建；加入后属于 cross-species transfer setting，与 fixed-target human-only benchmark 分开解释。
4. **prior 增益有控制实验支撑。** Permutation control 和 residual analysis 共同支持：提升来自 gene-specific cross-species signal，而非“多几列特征”或简单复述 prior。
5. **0.10 弱同源正则化是有边界的核心标签结果。** 在 12,307 个 model-eligible prediction genes 上，该 target 与 human no-Gejman PC1 的 `r=0.9982`；在单独定义的 10,768 个 mapped one-to-one ortholog pairs 上，shift RMSE 为 `0.065`。Cross-target evaluation 的 Pearson 变化为 `+0.0003`（95% CI `-0.0006` 至 `0.0012`），未检出原 human 标签可预测性下降，但因未预设等效界值，不构成正式等效性证明。

## 3. 推荐 claim 边界

- 可以说：本文提供了一个可复现的 study-aware label audit 和 cross-species prior-enhanced benchmark。
- 可以说：prior-enhanced 设定在明确边界下稳定优于 human-only baseline。
- 可以说：0.10 弱同源正则化构建了 human-dominant mammalian stability target，并通过标签距离、study-noise 和 cross-target 检查限定其解释。
- 固定目标 benchmark、cross-species transfer 和 ortholog-informed shrinkage target 应分开表述。
- 新序列预测应按 base sequence、compact regulatory 和 prior-enhanced 三个输入层级说明适用条件。

## 4. 为什么适合 BMB

- 研究对象是明确的生物学表型：mammalian mRNA stability / half-life。
- 研究方法强调定量标签重建、PCA 共识结构、留一研究影响分析、ortholog 一致性和受控预测比较。
- 文章的价值主要在于：把一个常被默认的标签问题变成可审计的数学与计算生物学问题。

## 5. 后续英文化时建议保留的三个关键词组

- `study-aware label audit`
- `human-only versus prior-enhanced benchmark`
- `weak ortholog-informed target shrinkage with cross-target validation`
