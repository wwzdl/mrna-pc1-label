# Study Noise Versus Ortholog-Regularized Label Shift

## Purpose

比较 `lambda=0.10` ortholog-regularized target 相对原 human no-Gejman PC1 的微小改变量，
与 human compendium 中不同 study 之间的实验差异是否处在同一量级。

## Key Results

- `lambda=0.10` 混合前后标签几乎不改变 human 标签主体结构：N=10768，Pearson=0.9979，RMSE=0.0646，MAE=0.0501。
- 在排除 Gejman 后的 normal-study 对照中，study mean label 两两比较的中位 Pearson 为 0.5303（IQR 0.3243-0.6110），中位 RMSE=0.9540，中位 MAE=0.7203。
- 因此，normal studies 之间的中位 label RMSE 约为 `lambda=0.10` 标签 shift 的 14.8 倍，中位 MAE 约为 14.4 倍。
- 若把每个 study label 相对 human no-Gejman PC1 的差异视为 study-specific residual，不同 study residual 两两相关的中位 Pearson 为 0.0796，说明误差结构并不高度一致。
- 单个 normal study 相对 human no-Gejman PC1 的中位 RMSE=0.7145，MAE=0.5429；分别约为 `lambda=0.10` shift 的 11.1 倍和 10.8 倍。
- 更保守地只保留 no-Gejman 且样本数不少于 2 的 7 个 study 时，study-pair 中位 RMSE=0.8963，MAE=0.6963，残差两两相关中位 Pearson=0.0337；结论仍不依赖单样本 study。
- 作为预处理敏感性检查，先对 no-Gejman 样本执行完整 Saluki-like preprocessing 后再按 study 取均值，study-pair 中位 Pearson 升至 0.5937，但中位 RMSE=0.9014，MAE=0.6779，仍远大于 `lambda=0.10` shift。

## Interpretation

- 可以说：`lambda=0.10` 是轻度 target regularization，而不是重写 human 标签。
- 更准确的写法不是“完全可以忽略”，而是“其幅度远小于 normal study-to-study variability，因此不太可能改变 human 标签的主要排序结构”。
- 这组结果可作为补充材料防守：ortholog regularization 的收益来自对 noisy consensus target 的轻度跨物种收缩，而不是把 human target 换成 mouse target。
