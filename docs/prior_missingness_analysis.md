# Mouse-Prior Missingness Analysis

## Purpose

评估 global prior-enhanced 模型在测试基因缺少 mouse prior 时的真实 OOF 表现。
该分析不重新训练模型，而是复用 10-fold multi-seed OOF predictions，并按 mouse-prior 覆盖状态分层计算指标。

## Key Results

- 全部 12,916 个基因上，real both-prior 模型 Pearson=0.8300±0.0007。
- coverage-aware fallback（有 mouse prior 用 prior-enhanced；两种 prior 都缺失时退回 human-only）Pearson=0.8318±0.0008。
- 两种 mouse prior 都缺失的 1347 个基因上，real both-prior 模型 Pearson=0.7555±0.0013，human-only 模型 Pearson=0.7695±0.0006。
- 两种 mouse prior 都可用的 11107 个基因上，real both-prior 模型 Pearson=0.8433±0.0010，human-only 模型 Pearson=0.7456±0.0014。

## Interpretation

- 缺少 mouse prior 的基因仍可预测，但不能期待 prior-enhanced 主结果的性能水平。
- 实际部署时可使用 coverage-aware fallback，而不是强行让 missing-prior genes 走 prior-enhanced 模型。
- missing_both_mouse_priors 子集样本量较小，指标方差和选择偏差需要谨慎解释。
- 该结果支持在论文中把模型应用范围写清楚：没有 mouse prior 时可退回 human sequence/regulatory 预测；最高性能依赖可用 ortholog prior。
