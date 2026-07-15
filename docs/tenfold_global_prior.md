# Ten-Fold Robustness Benchmark

## Purpose

将全局 cross-species prior-enhanced 主结果从 5-fold OOF 扩展到 10-fold gene-level OOF，
并用多个 random seeds 检查结果是否依赖特定 gene split。

## Key Results

- `compact_all_global`: Pearson=0.7478±0.0012, R2=0.5576±0.0017.
- `compact_all_plus_both_mouse_priors_global`: Pearson=0.8300±0.0007, R2=0.6887±0.0012.
- shuffled both-prior control: Pearson=0.7483±0.0010, R2=0.5583±0.0014.

## Interpretation

- 10-fold 多 seed 结果用于量化 gene split sensitivity。
- 该结果属于固定 Saluki human PC1 目标下的 cross-species transfer，不是 human-only prediction。
- shuffled control 回到 human-only 附近，支持增益依赖 gene-specific mouse prior，而不是特征列数或 missingness pattern。
- 各 seed 的 OOF prediction 仅用于 split-sensitivity 汇总，不称为模型融合。
