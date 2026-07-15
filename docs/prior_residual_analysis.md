# Prior Residual Analysis

## Purpose

该分析把 mouse-prior 可解释部分与 human compact_all 特征的额外修正分开估计，
用于回答 dual-prior transfer model 是否只是读取 mouse prior。

## Reproducible command

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.prior_residual_analysis
```

## Algorithm

1. 仅保留同时具有重建 mouse prior 和 Saluki mouse prior 的 `both_priors_available` genes。
2. 在每个 outer fold 中，用两列 mouse priors 训练 RidgeCV，并对 held-out fold 得到 `prior_only_prediction`。
3. 在训练折内计算 `Saluki human PC1 - prior_prediction`，作为 residual target。
4. 用 human `compact_all` features 训练 XGBoost/CUDA 去预测 residual target。
5. held-out fold 的最终预测为 `prior_only_prediction + compact_residual_prediction`。

## Key results

- N=11107, 5-fold OOF, random_state=42.
- prior-only linear: Pearson=0.7919, R2=0.6271.
- compact_all only: Pearson=0.7386, R2=0.5436.
- prior + compact residual: Pearson=0.8262, R2=0.6748.
- Remaining variance explained after prior: 12.79%.

## Interpretation

该结果支持的边界是：mouse prior 本身很强，但 human sequence/regulatory features 仍能在 prior 之外提供额外修正；
它不能被写成 pure sequence-only result，也不应被解释为 prior 完全没有贡献。
