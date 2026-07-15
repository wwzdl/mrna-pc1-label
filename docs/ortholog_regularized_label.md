# Ortholog-Informed Target-Shrinkage Benchmark

## Purpose

测试在 human no-Gejman PC1 标签中加入轻度 mouse ortholog-informed target shrinkage 后，
是否能在 prior-enhanced 预测设定下超过直接训练 Saluki human PC1。

## Key Results

- 最高 prior-enhanced Pearson(target)：`orthoreg_reconstructed_mouse_pc1_0.1` = 0.8551±0.0007。
- 最保守、仍高度接近 Saluki human PC1 的候选：`Saluki human PC1`，target-vs-Saluki Pearson=1.0000±0.0000，ortholog Pearson=0.7979±0.0000。
- 对照：Saluki human PC1 在同一 prior-enhanced benchmark 中 Pearson(target)=0.837。

## Interpretation

- 这些标签明确超过了 Saluki human PC1 的 prior-enhanced target prediction，但标签本身包含 mouse ortholog-informed target shrinkage。
- 因此它们适合被描述为 cross-species shrinkage targets，而不是 pure human sequence-only targets。

## Controls

- `Saluki human PC1`：real prior=0.8394±0.0007，pure human=0.7542±0.0007，shuffled prior=0.7544±0.0011，real-shuffled delta=0.0850±0.0004。
- `orthoreg_reconstructed_mouse_pc1_0.1`：real prior=0.8551±0.0007，pure human=0.7538±0.0010，shuffled prior=0.7532±0.0011，real-shuffled delta=0.1018±0.0012。

- real prior 明显高于 pure-human 和 shuffled-prior，说明提升依赖真实 gene-specific mouse prior，而不是额外列数或 missingness pattern。

## Figures

- `figures/ortholog_regularized_label/orthoreg_lambda_tradeoff.png`
- `figures/ortholog_regularized_label/orthoreg_lambda_tradeoff.pdf`
- `figures/ortholog_regularized_label/orthoreg_real_vs_controls.png`
- `figures/ortholog_regularized_label/orthoreg_real_vs_controls.pdf`
