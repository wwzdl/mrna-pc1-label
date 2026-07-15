# Ortholog Label Distance

## Purpose

量化 ortholog regularization 之前 human no-Gejman PC1 与 mouse reconstructed prior 的距离，
并检查 `lambda=0.10` 正则化后人鼠标签距离缩小多少。

## Key Results

- 正则化前：human no-Gejman PC1 vs mouse reconstructed prior，N=10768，Pearson=0.7893，Spearman=0.7825，MSE=0.4400，RMSE=0.6633，MAE=0.5124。
- `lambda=0.10` 后：orthoreg target vs mouse reconstructed prior，Pearson=0.8275，Spearman=0.8205，MSE=0.3607，RMSE=0.6006，MAE=0.4640。
- 相对正则化前，MSE 下降 18.0%，RMSE 下降 9.5%，MAE 下降 9.4%。
- `lambda=0.10` 标签相对 human no-Gejman PC1 的改变量较小：Pearson=0.9979，RMSE=0.0646，MAE=0.0501。
- 作为参考，Saluki human PC1 vs Saluki mouse prior：Pearson=0.7979，RMSE=0.6464，MAE=0.5026。

## Interpretation

- 所有距离指标都在 z-scored PC1 单位下计算，因此 MSE/MAE/RMSE 可直接解释为标准化标签尺度上的差异。
- `lambda=0.10` 没有把 human 标签大幅改成 mouse 标签，而是让 human 标签向 mouse ortholog 信号轻度移动。
- 该结果适合作为 ortholog-regularized target 的补充防守：正则化前人鼠标签存在明显但非随机的差距，轻度正则化提高相关性并降低距离，同时仍保持接近 human 标签。
