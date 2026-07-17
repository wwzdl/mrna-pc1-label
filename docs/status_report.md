# mRNA PC1 Label 项目状态报告

更新时间：2026-07-17

## 当前阶段

Paper1 已进入 *Bulletin of Mathematical Biology* 预投稿定稿阶段。活跃稿件只位于 `manuscript/bmb_submission/`；旧 TCBB 分支、旧 BMB 快照、Pro 样图和临时 QA 输出已移入工作区级 `workspace_archive/2026-07-13_bmb_release_cleanup/`。

中文主文/补充材料用于作者内部核对，英文主文/补充材料为当前投稿源稿。两种语言的结构、数值、图号和 claim boundaries 已同步。

## 三项核心贡献

1. **Study-aware label audit**

   Leave-one-study-out PC1 stability 在不使用 Saluki label 或下游模型分数的情况下，将 `Gejman` 排为 primary sample-weighted analysis 中的第一影响来源。保留 Gejman、删除 15 个非 Gejman 样本的条件性比较与两类 study-balanced estimators 表明，其几何影响部分来自样本数杠杆；Saluki human PC1 agreement 和 human-mouse ortholog concordance 的同步正向变化未被 500 次比较性删除复制。

2. **Fixed-target human-only baseline and cross-species transfer**

   Saluki human PC1 始终固定为目标。Human-only `compact_all` 是透明基线；增加两类 mouse ortholog priors 的模型属于 cross-species transfer。Fold-wise permutation 打断 gene-prior identity 后增益消失，支持增益来自 gene-specific cross-species signal。

3. **Ortholog-informed target shrinkage at `lambda = 0.10`**

   `0.10 ortholog-informed shrinkage target` 是一个独立核心标签构建结果，而非辅助调参。它与 human no-Gejman PC1 几乎重合，其偏移远小于普通 study-level variability；cross-target evaluation 未检出原 human label 可预测性下降。该目标是 human-dominant mammalian stability target，不是 pure human ground truth 或 Saluki human PC1 的替代排名。

## 关键数值

### Label audit

- full reconstructed human PC1 vs Saluki human PC1：`r = 0.948`
- no-Gejman reconstructed PC1 vs Saluki human PC1：`r = 0.983`
- 13,265 个共同基因上的配对增益：`0.0314`（95% CI `0.0298-0.0330`）
- 12,592 对 one-to-one ortholog 的 human-mouse Pearson：`0.741 -> 0.755`
- ortholog 配对增益：`0.0147`（95% CI `0.0111-0.0181`）
- 500 次保留 Gejman、删除 15 个非 Gejman 样本的条件性同样本量比较：PC1 stability 分布中位数 `0.9438`（中央 95% comparator range `0.9319-0.9561`），Gejman observed `0.9550`，经验低尾比例 `0.964`
- Saluki agreement gain 与 ortholog gain 均超过 500 次条件性比较结果，对应最小经验高尾比例 `1/501=0.002`；该比例不是经 study 选择校正的 whole-study p 值
- Dynamic/fixed universe、coverage 3/5/10、PCA rank 3/5/10 与 median-imputation 共 12 个设定均将 Gejman 排第 1；equal-study weighting 与 study-mean collapse 均将其排第 3

### Fixed-target prediction

- global prediction universe：`N = 12,916`
- repeated 10-fold x 3 seeds human-only：`r = 0.748 +/- 0.001`
- real dual-prior transfer：`r = 0.830 +/- 0.001`
- shuffled-prior control：`r = 0.748 +/- 0.001`
- seed-averaged real vs shuffled paired gain：`0.0794`（95% CI `0.0735-0.0855`）
- both-priors subset：`N = 11,107`，仅用于 coverage/residual analysis，不是主训练宇宙
- cross-fitted residual analysis：prior-only `r = 0.792`, `R2 = 0.627`；prior + human residual correction `r = 0.826`, `R2 = 0.675`；human features 解释 prior 剩余方差的 `12.9%`

### Ortholog-informed shrinkage target

- target-construction universe：`N = 12,644`
- model-eligible prediction universe：`N = 12,307`
- separately defined one-to-one label-geometry set：`N = 10,768 pairs`
- 12,307-gene prediction universe：target vs human no-Gejman PC1 `r = 0.9982`（95% CI `0.9981-0.9983`）
- 10,768 mapped one-to-one orthologs：target shift `RMSE = 0.065`, `MAE = 0.050`；target vs human no-Gejman PC1 `r = 0.9979`
- cross-target human-trained：`r = 0.7476 +/- 0.0010`
- cross-target regularized-trained：`r = 0.7480 +/- 0.0010`
- paired difference：`+0.0003`（95% CI `-0.0006 to 0.0012`）
- target prediction, human-only：`r = 0.7538 +/- 0.0010`
- primary audit-threshold + Saluki mouse priors：`r = 0.8551 +/- 0.0007`
- strict exact target-construction + Saluki mouse priors：`r = 0.8537 +/- 0.0004`
- strict no-direct control（保留 Saluki mouse prior 和全部 indicators）：`r = 0.8524 +/- 0.0005`
- exact direct-vector increment：`+0.0013`
- mask-only control：`r = 0.7538 +/- 0.0007`

## 分析宇宙账本

全文共有 12 个固定定义的核心集合，按三条支路组织：标签重建/审计 6 个（L1-L3、A1-A3），fixed-target prediction 3 个（P1-P3），target shrinkage 3 个（S1-S3）。另有 8 个只用于 coverage/sensitivity 的派生集合（D1-D8），因此机器账本共 20 行。其中只有 P2（12,916 genes）和 S2（12,307 genes）是两条预测任务的主训练/评估集合；P3（11,107 genes）只用于 prior coverage 与 residual analysis。S2 与 S3（10,768 pairs）不是包含关系，交集为 10,682 human genes。逐集合定义见主文表 2 和补充表 S6，机器可读版本为 `results/bmb_analysis_universe_ledger.tsv`，校验入口为 `python3 scripts/audit_analysis_universes.py --check`。

没有把 10,682-gene S2/S3 overlap 强制作为全篇共同集合：它会从 P2 移除 2,234 genes（17.3%），并使 human-only benchmark 以 mouse mapping availability 为纳入条件。公平性由“同一任务内部使用相同 gene IDs 和 folds”保证，而不是强迫不同 estimand 共用一个 ortholog-restricted overlap。

## 任务边界

- base-sequence 模型只需已分区的 human 5'UTR/CDS/3'UTR 序列，可用于新序列，但不等于 `compact_all`。
- `compact_all` 依赖 Saluki datapack 中预计算的 CWCS、SeqWeaver 和 DeepRiPe blocks，适用于已注释 Ensembl human genes。
- prior-enhanced 模型另需 mouse ortholog prior；两类 prior 都缺失时应回退到 human-only compact model。
- Cross-species transfer 改变输入，ortholog-informed shrinkage 改变 target，两类数值不构成同一 leaderboard。

## 稿件与复现

- 中文主文：`manuscript/bmb_submission/bmb_main_manuscript_cn.md`
- 英文主文：`manuscript/bmb_submission/bmb_main_manuscript_en.md`
- 中文补充：`manuscript/bmb_submission/bmb_supplementary_material_cn.md`
- 英文补充：`manuscript/bmb_submission/bmb_supplementary_material_en.md`
- 定量 claim 配对预检：`scripts/audit_bmb_submission_package.py`
- 分析宇宙与集合关系审计：`scripts/audit_analysis_universes.py --check`
- 中英文一致性检查：`scripts/check_bmb_cn_consistency.py --lang both`
- 文档重建：`scripts/render_bmb_submission_bundle.sh`
- 整体预检：`scripts/preflight_bmb_release.sh`
- 上传候选目录组装：`scripts/assemble_bmb_upload_candidate.sh`
- 分阶段结果复现：`scripts/reproduce_bmb_key_results.sh labels|models|figures`

主文 Fig. 1-8 和补图 S1-S4 均由脚本绘制。当前活跃源保留 600 dpi PNG、SVG 和 600 dpi lossless TIFF；PDF 在稿件冻结后由同一脚本重建。活跃投稿链路不使用 Pro/AI 位图作为最终图件。

MOESM2/MOESM3 由脚本直接下载；`compact_all` 所需预计算 Saluki feature datapack 来自 DOI `10.5281/zenodo.6326409`（约 18.8 GB），因体积较大不由默认脚本自动下载，README 已给出放置结构。

## 当前预检状态

- BMB DOCX-only 自动审计：`0 failure / 0 warning / 165 checks`
- English abstract：`228 words`（目标 `150-250 words`）
- keywords：`6`
- 英文主文含 46 条参考文献，补充材料含 11 条；唯一并集为 47 条，其中 45 条带 DOI
- 45 条 DOI 均已通过 Crossref 解析；44 条通过题名/年份/卷自动比对，XGBoost 因 Crossref 仅返回缩写题名而按 ACM 正式记录人工核对。LightGBM 与 scikit-learn 两条无 DOI 记录已用 NeurIPS/JMLR 官方页面核对
- BibTeX/RIS 均仅保留上述 47 条实际使用记录，可由脚本同步导出
- 中英文 DOCX 均已从当前 Markdown 重建；英文主文含连续行号和页码字段，英文补充含页码字段
- 中英文正文各含 9 个右侧编号的 Word 原生可编辑 OMML 公式，双语补充材料各含补充式（S1）-（S2）；行内数学同样以 OMML 输出，PDF 路径使用同源矢量 SVG 公式
- 12 张 PNG 已完成逐图验证，图号、panel label、图注、分辨率、边界和可读性未见阻断问题；Fig. 1 的条件性统计已明确标为 geometry lower-tail proportion，Fig. 3 图例已重排并放大，Fig. 6-8 不再使用会夸大差异的截断柱状表达；最终 SI PDF 的逐页检查留到版本冻结后执行
- GitHub 独立仓库当前发布候选包含 `299` 个文件、`26.17 MiB`，最大单文件 `2.29 MB`；大矩阵、逐基因 OOF arrays、feature packages、DOCX、PDF 与 TIFF 投稿包继续忽略。当前 270 个 release 文本文件及约 21.95 MB Git patch history 的高置信凭据扫描均为 0 命中
- 关键 OOF 完整性审计通过：9 个 global、3 个 cross-target、18 个 target-shrinkage 结果文件均覆盖规范基因宇宙、gene ID 唯一且预测无缺失；residual 11,107 genes 覆盖 folds 1-5

## 发布产物状态

- 2026-07-13 的旧 SI PDF、上传候选目录、ZIP 和图件 PDF 已移入 `../workspace_archive/2026-07-14_presubmission_stale_exports/`，不再作为活动投稿版本。
- 当前活动 DOCX 已重建；最终英文 SI PDF、上传候选目录和 ZIP 尚未生成，避免在作者信息与 GitHub 版本冻结前反复打包。
- 最终源码冻结后的完整项目备份固定为 `../workspace_archive/mrna-pc1-label_paper_pca_presubmission_20260715.tar.zst`，其 SHA-256 保存在同名 `.sha256` 边车文件；完整性结果在项目源码冻结后记录。

## 仍需作者确认

1. 确认 Wenzhuo Wang 的 ORCID（如无 ORCID，可不填）；Ying Shao 的通讯邮箱 `yshao@dlmu.edu.cn` 与 ORCID `0000-0002-4056-5757` 已确认。
2. 两位作者确认 Author Contributions、作者顺序、no-funding statement 与正式英文单位名。
3. 选择软件许可证；当前已提供 `CITATION.cff`，但尚未替作者决定法律许可条款。
4. 生成最终英文 Online Resource 1 PDF，逐页检查后再组装上传目录和最终 ZIP。

GitHub 仓库 [wwzdl/mrna-pc1-label](https://github.com/wwzdl/mrna-pc1-label) 已公开，Data/Code Availability 已改为当前时态；本轮严格科学与投稿编辑审计版本标记为 `mRNA-PC1-label-v1.4.3`。

## 科学局限

- 当前 label audit 基于一个公开 compendium，缺少第二个独立 half-life resource 验证。
- Leave-one-study-out influence 同时包含 sample count 和 profile direction 差异。条件性同样本量删除比较和 study balancing 已表明 Gejman 的几何影响并非 sample-count-adjusted statistical outlier；主张限定为“样本数杠杆叠加有方向的 study effect”。Saluki agreement 量化对已发布处理选择的恢复，ortholog concordance 则提供作用不同的跨物种检查。
- `lambda = 0.10` 是保守经验权重，不是为了最大化 CV 而选定；它的可迁移性仍需独立数据检验。
- 本文没有重训 Saluki 的固定 test split，不声称 human-only model 已在同口径上超过 Saluki。
- Cross-validation 采用 gene-level 随机分折而非 paralog-family 分组；相关基因可能跨 fold，gene-family-held-out 泛化尚未检验。
- Gene/ortholog-pair bootstrap 以 genes 或 pairs 为重采样单位，未显式建模 pathway 或 paralog family 内相关性。
