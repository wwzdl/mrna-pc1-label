# Ortholog-Informed Target-Shrinkage Prior Ablation

## Purpose

测试 ortholog-informed shrinkage 标签的预测优势是否只来自模型直接读取同源 mouse prior 值。

## Results

### orthoreg_reconstructed_mouse_pc1_0.1

- `human_only`：Pearson(target)=0.7538±0.0010；Pearson(Saluki human PC1)=0.7533±0.0013。
- `mask_only`：Pearson(target)=0.7533±0.0005；Pearson(Saluki human PC1)=0.7530±0.0007。
- `alternate_saluki_mouse_prior`：Pearson(target)=0.8521±0.0005；Pearson(Saluki human PC1)=0.8347±0.0006。
- `no_direct_reconstructed_mouse_pc1`：Pearson(target)=0.8521±0.0006；Pearson(Saluki human PC1)=0.8348±0.0007。
- `direct_reconstructed_mouse_pc1`：Pearson(target)=0.8542±0.0004；Pearson(Saluki human PC1)=0.8359±0.0004。
- `both_mouse_priors`：Pearson(target)=0.8551±0.0007；Pearson(Saluki human PC1)=0.8366±0.0006。
- no-direct setting 保留替代 Saluki mouse PC1 后，target Pearson=0.8521±0.0006；完整 both-prior setting 为 0.8551±0.0007。
- no-direct setting 相对 Saluki human PC1 prior baseline 的增益为 0.0153±0.0006；both-prior 相对 no-direct 只再增加 0.0029。
- direct reconstructed mouse PC1 相对 no-direct setting 的增益为 0.0021。

## Interpretation

- 如果 mask-only 接近 human-only，说明 missingness/confidence 指示本身不能解释提升。
- 如果 no-direct setting 低于 both-prior setting，说明同源 mouse prior 值确实贡献了标签可预测性。
- 该结果用于限定论文主张：ortholog-informed shrinkage target 的优势属于 cross-species prior-enhanced setting，而不是 pure human-only sequence setting。

## Figures

