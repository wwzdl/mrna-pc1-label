# Candidate Label Competition

## Goal

在统一评估口径下比较 Saluki human PC1、Saluki 对齐重建标签，以及 study-balanced 标签，
看我们自己的标签能否在外部一致性、study 稳定性和下游可学习性上拿到更扎实的优势。

## Headline

- 综合最靠前的候选是 `Saluki human PC1`。
- ortholog 一致性最好的是 `Saluki human PC1`，Pearson=0.798。
- leave-one-study-out 平均稳定性最好的是 `Saluki human PC1`，Mean Pearson=0.999。
- pure human OOF 最好学的是 `Saluki human PC1`，Pearson(target)=0.750。
- prior-enhanced OOF 最好学的是 `Saluki human PC1`，Pearson(target)=0.837。

## How To Read

- 如果某个标签只是在“与 Saluki human PC1 的一致性”上更高，它更像是在逼近 Saluki 定义。
- 如果某个标签同时提升了 ortholog 一致性、study 稳定性和 pure/prior OOF，可更有底气说它是“更好的标签”。
