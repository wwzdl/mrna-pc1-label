# Ortholog-Informed Target-Shrinkage Prior Ablation

## Purpose

测试 ortholog-informed shrinkage 标签的预测优势是否只来自模型直接读取同源 mouse prior 值。

## Results

### orthoreg_reconstructed_mouse_pc1_0.1

- `human_only`：Pearson(target)=0.7538±0.0010；Pearson(Saluki human PC1)=0.7533±0.0013。
- `mask_only`：Pearson(target)=0.7538±0.0007；Pearson(Saluki human PC1)=0.7533±0.0009。
- `alternate_saluki_mouse_prior`：Pearson(target)=0.8521±0.0005；Pearson(Saluki human PC1)=0.8347±0.0006。
- `no_direct_reconstructed_mouse_pc1`：Pearson(target)=0.8524±0.0005；Pearson(Saluki human PC1)=0.8351±0.0005。
- `direct_reconstructed_mouse_pc1`：Pearson(target)=0.8525±0.0003；Pearson(Saluki human PC1)=0.8341±0.0004。
- `exact_target_and_saluki_priors`：Pearson(target)=0.8537±0.0004；Pearson(Saluki human PC1)=0.8354±0.0004。
- `primary_audit_and_saluki_priors`：Pearson(target)=0.8551±0.0007；Pearson(Saluki human PC1)=0.8366±0.0006。
- 严格 no-direct setting 保留 Saluki mouse PC1 和全部 indicators 后，target Pearson=0.8524±0.0005；加入精确 target-construction mouse PC1 后为 0.8537±0.0004。
- 严格 no-direct setting 相对 Saluki human PC1 prior baseline 的增益为 0.0130±0.0005；精确 target-construction prior 相对 no-direct 再增加 0.0013。

## Interpretation

- 如果 mask-only 接近 human-only，说明 missingness/confidence 指示本身不能解释提升。
- 严格 direct-vector 对比固定 Saluki mouse PC1 和全部 indicators，只检验精确 target-construction mouse PC1 的增量。
- 该结果用于限定论文主张：ortholog-informed shrinkage target 的优势属于 cross-species prior-enhanced setting，而不是 pure human-only sequence setting。

## Figures

- `figures/ortholog_regularized_label/prior_ablation/orthoreg_reconstructed_mouse_pc1_0p1_prior_ablation.png`
- `figures/ortholog_regularized_label/prior_ablation/orthoreg_reconstructed_mouse_pc1_0p1_prior_ablation.pdf`
