# mRNA Half-Life 实验记录

## 2026-07-14

### Figure 1 connector spacing and reference-base expansion

- 调整 `scripts/redraw_bmb_flow_figures.py`：所有流程箭头端点与圆角框之间保留一致空隙；Fig. 1a 六个节点整体左移，使节点组相对外框的左右留白对称；Fig. 1c 右侧三个节点统一下移并共用同一中心线，两条标签输入箭头改为严格水平且上下对称。
- 为正文 Fig. 3 的 human raw/processed PCA 补充与补图 S2 一致的 `a/b` panel labels，并同步中英文图注；中间修图阶段支持以 `BMB_SKIP_PDF=1` 跳过图件 PDF，最终预检时再统一刷新。
- 修复正文 Fig. 7 的 panel label 与长标题相碰问题：将 `a-d` 提到标题上方独立留白行，并增加 2 x 2 面板间距；数据、坐标范围和数值标注保持不变。
- 中英文主文参考文献统一为 44 条实际引用文献，重点补充 `Gejman` 原始研究、mRNA half-life 测量技术、跨研究异质性、低秩缺失值估计、quantile normalization 与 orthology 方法依据。
- 中英文补充材料参考文献由 5 条扩展为 10 条，并在 PCA 插补、quantile normalization 和 Ensembl ortholog 定义处建立明确引用。
- 项目 BibTeX/RIS 库由 40 条扩展为 52 条唯一记录；RIS 由 `scripts/export_bib_to_ris.py` 重新生成。

## 2026-07-13

### BMB pre-submission convergence, reproducibility, and release audit

- 将中英文主文收敛为三项边界清楚的核心贡献：
  - Saluki-label-independent study-influence screening + reference/ortholog validation
  - fixed Saluki human PC1 target 下的 human-only baseline 和 cross-species transfer
  - 通过 label geometry、study-noise comparison 与 cross-target evaluation 验证的 `0.10 ortholog-informed shrinkage target`
- 对审稿高风险点做了定稿处理：
  - `Gejman` 称为 dominant study influence，不称为样本数校正后的 statistical outlier
  - 筛选只用 PC1 stability，Saluki agreement 和 ortholog concordance 仅用于验证，不使用下游模型分数
  - mouse priors 明确为不参与 human target 构建的外部 mouse-side covariates
  - `lambda = 0.10` 保留为独立核心亮点，同时明确它改变 target，不参与 fixed-target leaderboard
- 新增 paired bootstrap 定量验证并写入主文：
  - Gejman common-gene agreement gain `0.0314` (95% CI `0.0298-0.0330`)
  - one-to-one ortholog gain `0.0147` (95% CI `0.0111-0.0181`)
  - real-vs-shuffled prior gain `0.0794` (95% CI `0.0735-0.0855`)
  - cross-target delta `+0.0003` (95% CI `-0.0006 to 0.0012`)
- 完成 30 条 DOI 文献的 Crossref 题名核对，LightGBM 条目用 NeurIPS 官方页验证并补入 URL。
- BibTeX 去除重复 key；新增 `scripts/export_bib_to_ris.py`，将当时的 40 条项目文献完整导出为 EndNote-compatible RIS（2026-07-14 扩展为 52 条）。
- 将正式英文单位名统一为：
  - `Interdisciplinary Research Center for Biology and Chemistry`
  - `Laboratory of Molecular Modeling and Design`
- `render_bmb_submission_bundle.sh`、`preflight_bmb_release.sh` 和 `assemble_bmb_upload_candidate.sh` 改为从脚本位置定位仓库根目录，支持将 `paper_pca/` 独立克隆为 `mrna-pc1-label`。
- 投稿预检扩展为同时验证核心数值、正式单位名、AI 披露、BibTeX/RIS 一致性、DOCX 行号/页码、SI metadata 和 600 dpi lossless TIFF。
- 主文和补充材料逐页视觉 QA 已完成；活跃 Fig. 1-8 和 Supplementary Fig. S1-S4 均为脚本生成图件。
- 旧 Pro 图、旧 zip、旧 BMB 快照、TCBB 分支和临时 QA 文件已移入 `workspace_archive/2026-07-13_bmb_release_cleanup/`。
- 清理公开代码中的早期结果备注：不再把 cross-species transfer 写成“超过 Saluki”，不再把 seed/weight sensitivity 称为模型融合；统一强调 fixed target、gene-prior identity control 和适用边界。历史日志中的旧措辞仅保留为分析演进记录，不代表当前稿件口径。

## 2026-07-04

### BMB terminology cleanup and package resync

- 目的：
  - 将 `paper1` BMB 路线下 Saluki 相关旧命名彻底统一，避免正文、补充材料、脚本和结果表之间出现 author 风格旧命名与开发期混杂术语。
  - 让英文 Markdown、英文 docx 和最终 submission zip 保持同步。
- 主要处理：
  - 活跃代码与结果表统一到 `saluki_human_pc1`、`saluki_mouse_prior`、`reconstructed mouse prior` 口径。
  - BMB 英文主文和补充材料中，将开发期说法改为提交版表述，例如：
    - `mouse Saluki-like PC1` -> `the reconstructed mouse prior`
    - `Saluki-like sequence models` -> `sequence models such as Saluki`
  - 重新导出：
    - `manuscript/bmb_submission/bmb_main_manuscript.docx`
    - `manuscript/bmb_submission/bmb_supplementary_material.docx`
  - 重新打包：
    - `manuscript/bmb_submission/bmb_submission_package.zip`
- 验证：
  - `rg` 检查 `paper_pca/manuscript/bmb_submission` 与 `paper_pca/docs`，活跃提交材料中不再命中旧的 Saluki author 风格命名。
  - 新 zip 不包含中文内部稿或 `~$*.docx` 临时锁文件。

## 2026-07-03

### Prior residual analysis reproducibility cleanup for BMB

- 目的：
  - 将主文 4.6 中 “prior 很强但 human 特征仍提供额外信息” 的 residual analysis 从文稿叙述整理为可复现脚本。
  - 明确区分 `mouse-prior component` 与 `human compact_all residual correction`，避免 strongest prior-enhanced result 被误读为简单读取 mouse prior。
- 新增代码：
  - `src/mrna_half_life_paper/prior_residual_analysis.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.prior_residual_analysis
```

- 新输出：
  - `results/prior_residual_analysis/summary.tsv`
  - `results/prior_residual_analysis/oof_predictions.tsv`
  - `results/prior_residual_analysis/summary.json`
  - `docs/prior_residual_analysis.md`
  - `manuscript/bmb_submission/tables/Supplementary_Table_S17_prior_residual_analysis.md`
- 关键结果：
  - `both_priors_available` genes：N=`11107`
  - prior-only RidgeCV OOF：Pearson=`0.791895`，R2=`0.627096`
  - human `compact_all` only：Pearson=`0.738645`，R2=`0.543551`
  - prior + compact residual：Pearson=`0.826224`，R2=`0.674774`
  - human features 在 prior 之外解释的剩余方差比例：`12.79%`
- 文稿更新：
  - BMB 中文主文 3.9 与 4.6 改为明确两阶段 OOF residual analysis。
  - BMB 中文补充材料 S8 新增算法、复现命令和补充表 S17。

### 10-fold robustness and ortholog label-distance update for BMB revision

- 目的：
  - 回应 Saluki split / test-fold 可比性质疑，将 global prior 和主推 `orthoreg_reconstructed_mouse_pc1_0.10` 从 5-fold OOF 补强为 10-fold × 多 seed 稳健性证据。
  - 量化 `lambda=0.10` 正则化前后 human-mouse ortholog 标签距离，解释“轻度混合”实际改变量。
- 新增/更新代码：
  - `src/mrna_half_life_paper/fold_robustness.py`
  - `src/mrna_half_life_paper/ortholog_regularized_label.py`
  - `src/mrna_half_life_paper/ortholog_prior_ablation.py`
  - `src/mrna_half_life_paper/ortholog_label_distance.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_regularized_label \
  --sources reconstructed_mouse_pc1 \
  --lambdas 0.1 \
  --include-saluki-baseline \
  --include-controls \
  --random-states 13 42 101 \
  --cv-splits 10 \
  --results-root paper_pca/results/ortholog_regularized_label_10fold_main

PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_prior_ablation \
  --source reconstructed_mouse_pc1 \
  --lambdas 0.1 \
  --random-states 13 42 101 \
  --cv-splits 10 \
  --saluki-prior-baseline 0.839426 \
  --saluki-pure-baseline 0.754235 \
  --results-root paper_pca/results/ortholog_prior_ablation_lambda0p1_10fold

PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.fold_robustness \
  --random-states 13 42 101 \
  --cv-splits 10 \
  --results-root paper_pca/results/tenfold_global_prior

PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_label_distance
```

- 新输出：
  - `results/ortholog_regularized_label_10fold_main/summary.tsv`
  - `results/ortholog_regularized_label_10fold_main/summary_by_label.tsv`
  - `results/ortholog_prior_ablation_lambda0p1_10fold/summary.tsv`
  - `results/ortholog_prior_ablation_lambda0p1_10fold/summary_by_setting.tsv`
  - `results/tenfold_global_prior/summary_by_setting.tsv`
  - `results/tenfold_global_prior/blend_by_weight.tsv`
  - `results/ortholog_label_distance/distance_summary.tsv`
  - `docs/ortholog_label_distance.md`
- 关键结果：
  - 10-fold global prior:
    - `compact_all_global` Pearson=`0.747830±0.001195`
    - `compact_all_plus_both_mouse_priors_global` Pearson=`0.829958±0.000748`
    - shuffled both-prior Pearson=`0.748302±0.001018`
    - best OOF blend Pearson=`0.832369±0.000651`
  - 10-fold orthoreg target:
    - `saluki_human_pc1` prior-enhanced Pearson=`0.839426±0.000735`
    - `orthoreg_reconstructed_mouse_pc1_0.10` prior-enhanced Pearson=`0.855074±0.000651`
    - `orthoreg_reconstructed_mouse_pc1_0.10` shuffled-prior Pearson=`0.753243±0.001079`
  - 10-fold no-direct ablation:
    - `no_direct_reconstructed_mouse_pc1` Pearson(target)=`0.852126±0.000575`
    - 相对同口径 Saluki-prior baseline 增益=`+0.012700±0.000575`
  - label distance:
    - 正则化前 human no-Gejman PC1 vs 采用 Saluki-like coverage 阈值重建的 mouse prior：Pearson=`0.789265`，RMSE=`0.663318`，MAE=`0.512389`
    - `lambda=0.10` 后 orthoreg target vs mouse PC1：Pearson=`0.827461`，RMSE=`0.600551`，MAE=`0.463980`
    - `lambda=0.10` target vs 原 human no-Gejman PC1：Pearson=`0.997891`，RMSE=`0.064563`，MAE=`0.050059`
- 文稿更新：
  - BMB 主文加入 10-fold robustness、orthoreg 10-fold result、label-distance 和 consensus-label 边界解释。
  - BMB 补充材料新增 `Supplementary Table S12-S15`。

## 2026-07-02

### Paper1 BMB submission package preparation

- 目的：将 `paper1` 从 `TCBB` 路线扩展到 *Bulletin of Mathematical Biology*，按 Springer 官方要求建立新的 submission package。
- 官方核对页面：
  - `https://link.springer.com/journal/11538/submission-guidelines`
  - `https://link.springer.com/journal/11538/aims-and-scope`
  - `https://link.springer.com/journal/11538/how-to-publish-with-us`
  - `https://link.springer.com/journal/11538`
- 新增/更新文件：
  - `manuscript/bmb_submission/README.md`
  - `manuscript/bmb_submission/bmb_main_manuscript.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_en.md`
  - `manuscript/bmb_submission/bmb_main_manuscript.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_en.md`
  - `manuscript/bmb_submission/bmb_supplementary_material.docx`
  - `manuscript/bmb_submission/bmb_title_page_en.md`
  - `manuscript/bmb_submission/bmb_title_page_en.docx`
  - `manuscript/bmb_submission/bmb_cover_letter_en.md`
  - `manuscript/bmb_submission/bmb_cover_letter_en.docx`
  - `manuscript/bmb_submission/bmb_submission_requirements.md`
  - `manuscript/bmb_submission/bmb_submission_checklist.md`
  - `manuscript/bmb_submission/bmb_submission_package.zip`
  - `scripts/render_bmb_markdown_with_references.py`
  - `scripts/render_markdown_to_docx.py`
  - `docs/paper1_bmb_submission_audit_2026-07-02.md`
- 处理内容：
  - 主稿切换为 BMB 标题、作者-年份引用、`Keywords`、`Statements and Declarations`。
  - supplement 切换为 Springer supplementary 口径。
  - 新建 BMB requirements/checklist/README。
  - 生成可交接 zip 包。
  - 修补 `markdown -> docx` 渲染脚本，自动去掉 YAML front matter。
- 验证：
  - 主稿 abstract=`223` 词，符合 BMB `150-250` 词要求。
  - keywords=`6`。
  - `bmb_submission/` 未再命中 `TCBB` / `IEEE` 残留短语。
  - `render_bmb_markdown_with_references.py`、`render_markdown_to_docx.py` 已 `py_compile` 通过。

### Paper1 TCBB submission package text/status cleanup

- 目的：基于交接文档继续推进 paper1 投稿包，将主文和 checklist 从内部 long-draft/投稿策略状态推进到更接近正式投稿状态。
- 读取：
  - `docs/codex_handoff.md`
  - `paper_pca/docs/status_report.md`
  - `paper_pca/docs/experiment_log.md`
  - `workspace_archive/2026-05-17_submission_cleanup/paper_pca/nouse/AGENTS.md`
- 更新文件：
  - `manuscript/tcbb_submission/tcbb_main_manuscript_en.md`
  - `manuscript/tcbb_submission/tcbb_main_manuscript.md`
  - `manuscript/tcbb_submission/tcbb_main_manuscript_en.docx`
  - `manuscript/tcbb_submission/tcbb_main_manuscript.docx`
  - `manuscript/tcbb_submission/tcbb_submission_checklist.md`
  - `manuscript/tcbb_submission/README.md`
  - `manuscript/tcbb_submission/tcbb_submission_package.zip`
  - `scripts/prepare_tcbb_submission_package.py`
  - `scripts/prepare_tcbb_author_materials.py`
  - `docs/paper1_submission_audit_2026-07-02.md`
- 处理内容：
  - 将英文主文中的内部投稿策略段落改为正式论文语气。
  - 同步作者、单位、基金、利益冲突和 ORCID 状态。
  - 重建主文 Word 草稿。
  - 更新投稿 zip 包中的主文、README 和 checklist。
- 验证：
  - 当前 Markdown 和 zip 内英文主文未再命中 `Before submission`、`TCBB Positioning`、`A reviewer may ask`、`should not oversell`、`remain placeholders` 等高风险内部短语。

## 2026-05-16

### MOESM3 Saluki matrix check

- 目的：单独验证作者处理后的 `MOESM3` 矩阵与官方 `saluki_human_pc1` 的一致性，确认 human 侧的主要差异是否来自 `Gejman`
- 代码：
  - `src/mrna_half_life_paper/moesm3_saluki_matrix_check.py`
- 命令：

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.moesm3_saluki_matrix_check
```

- 新输出：
  - `results/moesm3_saluki_matrix_check/summary.tsv`
  - `docs/moesm3_saluki_matrix_check.md`
- 关键结果：
  - human `all_samples + as_is_pc1` vs `saluki_human_pc1`：Pearson=`0.963045`
  - human `no_gejman + as_is_pc1` vs `saluki_human_pc1`：Pearson=`1.000000`
  - mouse `all_samples + as_is_pc1` vs `saluki_human_pc1`：Pearson=`1.000000`
- 解释：
  - `MOESM3` 本身已经非常接近Saluki 官方标签
  - human 侧剩余差异的核心来源就是 `Gejman`

### MOESM3 Saluki label downstream benchmark

- 目的：测试如果直接把 `MOESM3` Saluki 矩阵衍生标签当训练目标，下游模型是否更容易学，以及最终是否更能恢复官方 `saluki_human_pc1`
- 代码：
  - `src/mrna_half_life_paper/moesm3_saluki_label_benchmark.py`
- 命令：

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.moesm3_saluki_label_benchmark
```

- 新输出：
  - `results/moesm3_saluki_label_benchmark/summary.tsv`
  - `results/moesm3_saluki_label_benchmark/oof_predictions/*.tsv`
  - `docs/moesm3_saluki_label_benchmark.md`
- 共同评估口径：
  - 共同基因宇宙 `N=13532`
  - `base_sequence`：`526` 个特征
  - `compact_all`：`1802` 个特征
- 关键结果：
  - `compact_all + saluki_human_pc1`：Pearson(target)=`0.729390`，Pearson(Saluki ref)=`0.729390`
  - `compact_all + moesm3_all_samples_as_is_pc1`：Pearson(target)=`0.707227`，Pearson(Saluki ref)=`0.721250`
  - `compact_all + moesm3_no_gejman_rerun_pipeline_pc1`：Pearson(target)=`0.730410`，Pearson(Saluki ref)=`0.730275`
  - `base_sequence + saluki_human_pc1`：Pearson(target)=`0.719535`，Pearson(Saluki ref)=`0.719535`
  - `base_sequence + moesm3_no_gejman_zscore_only_pc1`：Pearson(target)=`0.719755`，Pearson(Saluki ref)=`0.719759`
- 解释：
  - 保留 `Gejman` 的 `MOESM3 all_samples` 三个变体整体都低于直接训练官方 `saluki_human_pc1`
  - `MOESM3 no_gejman as_is_pc1` 在共同子集上与 `saluki_human_pc1` 等价
  - 唯一出现轻微增益的是 `moesm3_no_gejman_rerun_pipeline_pc1`，但提升只有 `+0.000885`

### MOESM3 Saluki label prior benchmark

- 目的：检查在明确的 cross-species transfer setting 下，`MOESM3` 派生标签是否会比直接训练 `saluki_human_pc1` 更有优势
- 代码：
  - `src/mrna_half_life_paper/moesm3_saluki_label_prior_benchmark.py`
- 命令：

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.moesm3_saluki_label_prior_benchmark \
  --targets saluki_human_pc1 moesm3_no_gejman_rerun_pipeline_pc1 \
  --results-root results/moesm3_saluki_label_prior_benchmark_smoketest
```

- 新输出：
  - `results/moesm3_saluki_label_prior_benchmark_smoketest/summary.tsv`
  - `results/moesm3_saluki_label_prior_benchmark_smoketest/oof_predictions/*.tsv`
  - `docs/moesm3_saluki_label_prior_benchmark.md`
- 共同评估口径：
  - 共同基因宇宙 `N=13532`
  - `compact_all_global / +mouse_pc1 / +saluki_mouse_prior / +both priors`
- 关键结果：
  - `saluki_human_pc1 + compact_all_global`：Pearson(Saluki ref)=`0.729390`
  - `saluki_human_pc1 + compact_all_plus_mouse_pc1_global`：Pearson(Saluki ref)=`0.813902`
  - `saluki_human_pc1 + compact_all_plus_saluki_mouse_prior_global`：Pearson(Saluki ref)=`0.811956`
  - `saluki_human_pc1 + compact_all_plus_both_mouse_priors_global`：Pearson(Saluki ref)=`0.815917`
  - `moesm3_no_gejman_rerun_pipeline_pc1 + compact_all_plus_both_mouse_priors_global`：Pearson(Saluki ref)=`0.815588`
- 解释：
  - 带上 `mouse prior` 后，最优组合仍然是直接训练 `saluki_human_pc1`
  - `moesm3_no_gejman_rerun_pipeline_pc1` 只是在 human-only 下略优，在 strongest prior setting 下反而略低于 `saluki_human_pc1`
  - 这进一步说明：`MOESM3` 派生标签适合做控制实验，但不值得替换主线监督目标

## 2026-05-15

### Global prior permutation control

- 目的：验证 `global mouse prior` 的提升是否来自 gene-specific 跨物种信息，而不是“多加几列特征”或 missingness pattern
- 代码：
  - `src/mrna_half_life_paper/global_mouse_prior.py`
- 命令：

```bash
source scripts/activate_env.sh
python -m mrna_half_life_paper.global_mouse_prior
```

- 新输出：
  - `results/human_global_mouse_prior_permutation_control.tsv`
  - `manuscript/human_global_mouse_prior_note_cn.md`
- 关键结果：
  - `compact_all_global`: Pearson=`0.745425`，R2=`0.553911`
  - `compact_all_plus_both_mouse_priors_global`: Pearson=`0.829493`，R2=`0.687856`
  - `compact_all_plus_both_mouse_priors_shuffled_global`: Pearson=`0.744925`，R2=`0.553149`
- 解释：
  - shuffled 后性能几乎退回 `compact_all_global`
  - 说明真实提升主要来自 gene-specific mouse prior，而不是缺失指示特征本身

### TCBB package refresh

- 脚本：
  - `scripts/prepare_tcbb_submission_package.py`
- 命令：

```bash
source scripts/activate_env.sh
python scripts/prepare_tcbb_submission_package.py
```

- 更新内容：
  - 主文结果段加入 permutation control
  - 补充材料新增 `Supplementary Table S7`
  - 补充材料里原先“未来可加 shuffled-prior control”的表述，改成“当前已纳入”

## 2026-07-02

### Candidate label competition

- 目的：
  - 不再只问“我们的标签像不像 `saluki_human_pc1`”，而是在统一口径下比较多个候选标签是否能在外部一致性、study 稳定性和下游可学习性上体现优势。
- 代码：
  - `src/mrna_half_life_paper/label_competition.py`
  - `src/mrna_half_life_paper/pc1_consensus.py`
  - `src/mrna_half_life_paper/global_mouse_prior.py`
- 额外实现：
  - 新增 `compute_study_weighted_pc1()`，通过对每个 study 的样本列施加 `1/sqrt(n_study_samples)` 权重，降低大 study 对 PC1 的主导。
  - 修复 `global_mouse_prior.py` 中调控特征目录使用相对路径的问题，避免 benchmark 依赖特定工作目录。
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.label_competition
```

- 候选标签：
  - `saluki_human_pc1`
  - `saluki_like_full_pc1`
  - `saluki_like_no_gejman_pc1`
  - `study_weighted_full_pc1`
  - `study_weighted_no_gejman_pc1`
- 统一评估口径：
  - shared genes `N=12307`
  - shared ortholog pairs `N=10768`
  - pure-human setting：`compact_all` (`1802` features)
  - prior-enhanced setting：`compact_all_plus_both_mouse_priors_global` (`1808` features)
- 新输出：
  - `results/label_competition/overall_summary.tsv`
  - `results/label_competition/human_saluki_alignment.tsv`
  - `results/label_competition/stability_summary.tsv`
  - `results/label_competition/ortholog_summary.tsv`
  - `results/label_competition/pure_human_summary.tsv`
  - `results/label_competition/prior_enhanced_summary.tsv`
  - `results/label_competition/oof_predictions/*.tsv`
  - `docs/label_competition.md`
- 关键结果：
  - `saluki_human_pc1` 仍然是 overall 最强标签：
    - ortholog Pearson=`0.797936`
    - pure human OOF Pearson(target)=`0.749565`
    - prior-enhanced OOF Pearson(target)=`0.836801`
  - 最强的非官方下游候选是 `saluki_like_no_gejman_pc1`：
    - pure human OOF Pearson(target)=`0.743907`
    - prior-enhanced OOF Pearson(target)=`0.829931`
  - `study_weighted_full_pc1` 的主要亮点是全量稳健性：
    - 相比 `saluki_like_full_pc1`，worst-case leave-one-study-out Pearson 从 `0.958250` 提升到 `0.983442`
    - ortholog Pearson 从 `0.772333` 提升到 `0.794972`
    - 说明它在**不删 Gejman 的前提下**，显著削弱了 full-data 重建对异常 study 的脆弱性
  - 但 `study_weighted_full_pc1` 的下游 OOF 仍低于 `saluki_human_pc1`：
    - pure human Pearson(target)=`0.740260`
    - prior-enhanced Pearson(target)=`0.826990`
- 当前解释：
  - 还不能声称“我们的标签全面优于 `saluki_human_pc1`”
  - 但可以更有底气地说：`study_weighted_full_pc1` 是一个**有明确稳健性优势的候选标签**
  - 如果 paper1 要体现我们在标签处理上的独立贡献，这条“study-balanced full label improves robustness without manual pruning” 是目前最有价值的角度

### Custom-label expansion quick screen

- 目的：
  - 在 `study_weighted_full_pc1` 之后，继续尝试更稳健的 custom label，优先找“非官方但更平衡”的候选，而不是直接重复 `saluki_human_pc1`。
- 代码：
  - `src/mrna_half_life_paper/pc1_consensus.py`
- 新增原型：
  - `compute_group_collapsed_pc1()`：先按 `study` 折叠样本，再做 PC1
  - `compute_jackknife_pc1()`：对 leave-one-study-out 标签做 jackknife 平均
- 快筛结果（human no-Gejman / mouse full）：
  - `group_collapsed_no_gejman`：
    - human vs Saluki Pearson=`0.968766`
    - ortholog Pearson=`0.785752`
    - mean stability Pearson=`0.998858`
    - min stability Pearson=`0.998191`
  - `jackknife_no_gejman`：
    - human vs Saluki Pearson=`0.987959`
    - ortholog Pearson=`0.789284`
    - mean stability Pearson=`0.998499`
    - min stability Pearson=`0.991137`
- 融合筛选：
  - 以 z-score 后的非官方标签做线性混合，最佳简单两路混合为：
    - `0.6 * saluki_like_no_gejman + 0.4 * group_collapsed_no_gejman`
  - 该混合标签的快筛表现：
    - ortholog Pearson 约=`0.7918`
    - stability mean Pearson 约=`0.9988`
    - stability min Pearson 约=`0.9949`
- 对该混合标签的 OOF 追测：
  - pure human：
    - Pearson(target)=`0.742096`
    - Pearson(Saluki)=`0.749188`
  - prior-enhanced：
    - Pearson(target)=`0.827474`
    - Pearson(Saluki)=`0.834919`
- 当前解释：
  - 这个 hybrid 是当前“我们自己的标签”里最平衡的一个：
    - 比 `saluki_like_no_gejman` 更稳
    - ortholog 也略高
    - pure-human 下对 Saluki human PC1 的恢复几乎持平
  - 但它仍未超过 `saluki_human_pc1`
  - 因此下一轮若还要冲“超越 Saluki 标签”，应考虑更强的标签模型，而不是只在线性 PC1 家族里微调

### Ortholog-informed shrinkage label benchmark

- 目的：
  - 直接追求“我们的标签作为预测目标优于 `saluki_human_pc1`”。
  - 不再只在人类 PC1 家族里微调，而是在 human no-Gejman PC1 中加入 mouse ortholog-informed target shrinkage。
- 代码：
  - `src/mrna_half_life_paper/ortholog_regularized_label.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_regularized_label
```

- 新输出：
  - `results/ortholog_regularized_label/summary.tsv`
  - `results/ortholog_regularized_label/targets/*.tsv`
  - `results/ortholog_regularized_label/oof_predictions/*.tsv`
  - `docs/ortholog_regularized_label.md`
- 设计：
  - 基础 human label：`saluki_like_no_gejman_pc1`
  - mouse regularizer source：
    - `reconstructed_mouse_pc1`
    - `reconstructed_mouse_pc1`
  - label 定义：
    - 对 human 和 mouse label 做 z-score
    - 对有 one-to-one ortholog 的 human genes：`(1-lambda) * human_pc1 + lambda * mouse_pc1`
    - 对无 ortholog 的 human genes：保留 human PC1
  - lambdas：`0.05, 0.10, 0.15, 0.20, 0.30`
  - 下游评估：`compact_all_plus_both_mouse_priors_global` prior-enhanced OOF
- 关键对照：
  - 在同一 shared-gene 口径下，`saluki_human_pc1` prior-enhanced Pearson(target)=`0.836801`
- 关键结果：
  - 保守增强候选：
    - `orthoreg_reconstructed_mouse_pc1_0.05`
    - target-vs-Saluki Pearson=`0.988113`
    - ortholog Pearson=`0.808637`
    - prior-enhanced Pearson(target)=`0.841760`
  - 中等增强候选：
    - `orthoreg_reconstructed_mouse_pc1_0.10`
    - target-vs-Saluki Pearson=`0.987384`
    - ortholog Pearson=`0.827461`
    - prior-enhanced Pearson(target)=`0.853598`
  - upper-bound cross-species target:
    - `orthoreg_reconstructed_mouse_pc1_0.30`
    - target-vs-Saluki Pearson=`0.975043`
    - ortholog Pearson=`0.895449`
    - prior-enhanced Pearson(target)=`0.896171`
- 当前解释：
  - 这是目前第一条真正支持“我们的标签带来更好预测结果”的结果线
  - 但它不是 pure human sequence-only 标签，而是 **cross-species shrinkage target**
  - 最适合论文中的定位：
    - `saluki_human_pc1` 是官方 pure benchmark target
    - `orthoreg_reconstructed_mouse_pc1_0.05/0.10` 是保守 cross-species shrinkage target
    - `orthoreg_reconstructed_mouse_pc1_0.30` 是 upper-bound prior-enhanced target
  - 后续必须补充：
    - pure-human 对照，证明不是所有提升都来自 mouse prior 直接可见
    - permutation control，证明不是只因为目标里混入了可预测的 mouse feature
    - 多 seed，确认不是单次 OOF 波动

### Ortholog-informed shrinkage label controls

- 目的：
  - 验证 `ortholog_regularized_label` 的胜利是否依赖真实 gene-specific mouse prior。
  - 同时检查 pure-human 特征是否也能预测这些标签。
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_regularized_label \
  --sources reconstructed_mouse_pc1 \
  --lambdas 0.05 0.1 0.3 \
  --include-controls \
  --results-root paper_pca/results/ortholog_regularized_label_controls
```

- 新输出：
  - `results/ortholog_regularized_label_controls/summary.tsv`
  - `results/ortholog_regularized_label_controls/targets/*.tsv`
  - `results/ortholog_regularized_label_controls/oof_predictions/*pure_human.tsv`
  - `results/ortholog_regularized_label_controls/oof_predictions/*prior_enhanced.tsv`
  - `results/ortholog_regularized_label_controls/oof_predictions/*prior_shuffled.tsv`
- 关键结果：
  - `orthoreg_reconstructed_mouse_pc1_0.05`：
    - real prior Pearson(target)=`0.841760`
    - pure human Pearson(target)=`0.747794`
    - shuffled prior Pearson(target)=`0.746802`
  - `orthoreg_reconstructed_mouse_pc1_0.10`：
    - real prior Pearson(target)=`0.853598`
    - pure human Pearson(target)=`0.750708`
    - shuffled prior Pearson(target)=`0.749881`
  - `orthoreg_reconstructed_mouse_pc1_0.30`：
    - real prior Pearson(target)=`0.896171`
    - pure human Pearson(target)=`0.755918`
    - shuffled prior Pearson(target)=`0.754076`
- 当前解释：
  - ortholog-informed shrinkage 标签的高预测性能主要来自真实 gene-specific mouse prior。
  - pure-human 特征和 shuffled-prior 都不能复现 real-prior 的大幅提升。
  - 这强化了文章边界：该结果应作为 cross-species regularized / prior-enhanced 标签结果，而不是 pure sequence-only 标签结果。

### Ortholog-informed shrinkage label multi-seed controls and figures

- 目的：
  - 将 ortholog-informed shrinkage target 的正结果从单 seed 升级为多 seed 稳健性证据。
  - 同时输出可直接进入论文或补充材料的 lambda tradeoff 图和 real-vs-control 图。
- 代码更新：
  - `src/mrna_half_life_paper/ortholog_regularized_label.py`
  - 新增 `--random-states`
  - OOF 文件名加入 `seed`
  - 新增 `summary_by_label.tsv`
  - 新增自动图件输出到 `figures/ortholog_regularized_label/`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_regularized_label \
  --sources reconstructed_mouse_pc1 \
  --lambdas 0.05 0.1 0.3 \
  --include-controls \
  --random-states 13 42 101 \
  --results-root paper_pca/results/ortholog_regularized_label_multiseed
```

- 新输出：
  - `results/ortholog_regularized_label_multiseed/summary.tsv`
  - `results/ortholog_regularized_label_multiseed/summary_by_label.tsv`
  - `results/ortholog_regularized_label_multiseed/summary.json`
  - `results/ortholog_regularized_label_multiseed/targets/*.tsv`
  - `results/ortholog_regularized_label_multiseed/oof_predictions/*seed*prior_enhanced.tsv`
  - `results/ortholog_regularized_label_multiseed/oof_predictions/*seed*pure_human.tsv`
  - `results/ortholog_regularized_label_multiseed/oof_predictions/*seed*prior_shuffled.tsv`
  - `figures/ortholog_regularized_label/orthoreg_lambda_tradeoff.png`
  - `figures/ortholog_regularized_label/orthoreg_lambda_tradeoff.pdf`
  - `figures/ortholog_regularized_label/orthoreg_real_vs_controls.png`
  - `figures/ortholog_regularized_label/orthoreg_real_vs_controls.pdf`
- 多 seed 关键结果：
  - `orthoreg_reconstructed_mouse_pc1_0.05`
    - target-vs-Saluki Pearson=`0.988113`
    - ortholog Pearson=`0.808637`
    - real prior Pearson(target)=`0.842083±0.000420`
    - pure human Pearson(target)=`0.748026±0.000213`
    - shuffled prior Pearson(target)=`0.747997±0.001043`
    - real-shuffled delta=`0.094087±0.000877`
  - `orthoreg_reconstructed_mouse_pc1_0.10`
    - target-vs-Saluki Pearson=`0.987384`
    - ortholog Pearson=`0.827461`
    - real prior Pearson(target)=`0.853833±0.000216`
    - pure human Pearson(target)=`0.751078±0.000575`
    - shuffled prior Pearson(target)=`0.750853±0.001010`
    - real-shuffled delta=`0.102980±0.000896`
  - `orthoreg_reconstructed_mouse_pc1_0.30`
    - target-vs-Saluki Pearson=`0.975043`
    - ortholog Pearson=`0.895449`
    - real prior Pearson(target)=`0.896018±0.000147`
    - pure human Pearson(target)=`0.755289±0.000591`
    - shuffled prior Pearson(target)=`0.754531±0.000451`
    - real-shuffled delta=`0.141488±0.000598`
- 当前解释：
  - 三个 random seed 下，全部 ortholog-informed shrinkage targets 的 real-prior Pearson(target) 均稳定超过 `saluki_human_pc1` prior-enhanced baseline `0.836801`。
  - `0.05` 是最保守胜利：与 Saluki human PC1 相关最高，增益较小但稳定。
  - `0.10` 是当前最适合主推的平衡点：仍高度接近 Saluki human PC1，同时 ortholog 一致性和可预测性明显提高。
  - `0.30` 可作为 upper-bound cross-species target：预测性能最高，但与 Saluki human PC1 的距离也更大。
  - pure-human 和 shuffled-prior 控制仍显著低于 real-prior，说明增益依赖真实 gene-specific mouse prior。
  - 仍需在论文中明确边界：这是 cross-species regularized / prior-enhanced label，不是 pure human sequence-only target。

### Ortholog-informed shrinkage prior ablation for the main lambda=0.10 label

- 目的：
  - 回答审稿人可能提出的关键质疑：`orthoreg_reconstructed_mouse_pc1_0.10` 的提升是否只是因为模型直接读取了构造标签时使用的同源 `reconstructed_mouse_pc1` prior 值。
  - 区分 human-only、missingness/confidence 指示、替代 Saluki-mouse prior、direct reconstructed mouse PC1 prior 和 both-prior setting 的贡献。
- 代码：
  - `src/mrna_half_life_paper/ortholog_prior_ablation.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.ortholog_prior_ablation \
  --source reconstructed_mouse_pc1 \
  --lambdas 0.1 \
  --random-states 13 42 101 \
  --results-root paper_pca/results/ortholog_prior_ablation_lambda0p1
```

- 新输出：
  - `results/ortholog_prior_ablation_lambda0p1/summary.tsv`
  - `results/ortholog_prior_ablation_lambda0p1/summary_by_setting.tsv`
  - `results/ortholog_prior_ablation_lambda0p1/summary.json`
  - `results/ortholog_prior_ablation_lambda0p1/oof_predictions/*.tsv`
  - `figures/ortholog_regularized_label/prior_ablation/orthoreg_reconstructed_mouse_pc1_0p1_prior_ablation.png`
  - `figures/ortholog_regularized_label/prior_ablation/orthoreg_reconstructed_mouse_pc1_0p1_prior_ablation.pdf`
  - `docs/ortholog_prior_ablation.md`
- 关键结果：
  - `human_only`：Pearson(target)=`0.751078±0.000575`
  - `mask_only`：Pearson(target)=`0.751310±0.000077`
  - `alternate_saluki_mouse_prior`：Pearson(target)=`0.851225±0.000609`
  - `no_direct_reconstructed_mouse_pc1`：Pearson(target)=`0.851328±0.000357`
  - `direct_reconstructed_mouse_pc1`：Pearson(target)=`0.853249±0.000410`
  - `both_mouse_priors`：Pearson(target)=`0.853833±0.000216`
- 当前解释：
  - `mask_only` 几乎等于 `human_only`，因此提升不是 ortholog availability、missingness 或 confidence 指示造成的。
  - `no_direct_reconstructed_mouse_pc1` 在不使用 direct `reconstructed_mouse_pc1` value 的情况下仍达到 `0.851328±0.000357`，相对 Saluki-prior baseline `0.836801` 仍有 `+0.014527±0.000357` 增益。
  - direct `reconstructed_mouse_pc1` prior 相对 no-direct 只增加约 `+0.0019`，both-prior 相对 no-direct 只增加约 `+0.0025`。
  - 这显著增强了论文防守：主推标签 `orthoreg_reconstructed_mouse_pc1_0.10` 的可预测性优势主要来自可替代的 cross-species gene-specific signal，而不是单纯把构造标签的同一列 mouse prior 交给模型。
  - 仍需保留边界：它是 cross-species prior-enhanced label，不是 pure human-only 标签。

### BMB manuscript integration of ortholog-informed shrinkage targets

- 目的：
  - 将 ortholog-informed shrinkage target 与 prior-ablation 结果从开发记录推进到 BMB 投稿材料。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_en.md`
  - `manuscript/bmb_submission/bmb_main_manuscript.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_en.md`
  - `manuscript/bmb_submission/bmb_supplementary_material.docx`
  - `manuscript/bmb_submission/README.md`
  - `manuscript/bmb_submission/bmb_submission_package.zip`
- 新增补充图：
  - `manuscript/bmb_submission/figures/supplement/FigS12_orthoreg_lambda_tradeoff.{png,tiff}`
  - `manuscript/bmb_submission/figures/supplement/FigS13_orthoreg_real_vs_controls.{png,tiff}`
  - `manuscript/bmb_submission/figures/supplement/FigS14_orthoreg_prior_ablation.{png,tiff}`
- 新增补充表：
  - `manuscript/bmb_submission/tables/Supplementary_Table_S10_orthoreg_multiseed_targets.md`
  - `manuscript/bmb_submission/tables/Supplementary_Table_S11_orthoreg_prior_ablation.md`
- 渲染命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python paper_pca/scripts/render_bmb_markdown_with_references.py
python paper_pca/scripts/render_markdown_to_docx.py \
  paper_pca/manuscript/bmb_submission/bmb_main_manuscript_en.md \
  paper_pca/manuscript/bmb_submission/bmb_main_manuscript.docx
python paper_pca/scripts/render_markdown_to_docx.py \
  paper_pca/manuscript/bmb_submission/bmb_supplementary_material_en.md \
  paper_pca/manuscript/bmb_submission/bmb_supplementary_material.docx
```

- zip 更新：
  - `bmb_submission_package.zip` 已重建为英文投稿材料包。
  - 验证 zip 内已包含 `FigS12-FigS14` 和 `Supplementary_Table_S10-S11`。
  - 验证 zip 内不再包含旧的 `_cn` 内部检查稿或 Word 临时文件。
- 当前解释：
  - 主文只简短呈现 ortholog-informed shrinkage target 的结论，完整 multi-seed、real-vs-shuffled、pure-human 和 prior-ablation 证据放入补充材料。
  - 这让“我们的标签更好”进入投稿稿件，同时保留清晰边界：cross-species prior-enhanced target，不是 pure human-only sequence target。

### Prior-missingness analysis for deployment fallback

- 目的：
  - 回答实际应用问题：如果新 human gene 没有相关 mouse prior，global prior-enhanced 模型还能不能可靠预测。
  - 复用 10-fold multi-seed OOF predictions，按 mouse-prior 覆盖状态分层评估，不重新训练模型。
- 代码：
  - `src/mrna_half_life_paper/prior_missingness_analysis.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.prior_missingness_analysis
```

- 新输出：
  - `results/prior_missingness_analysis/summary_by_seed.tsv`
  - `results/prior_missingness_analysis/summary_by_coverage.tsv`
  - `docs/prior_missingness_analysis.md`
- 关键结果：
  - 全部 `12916` 个基因：real both-prior Pearson=`0.8300±0.0007`
  - coverage-aware fallback（有 mouse prior 用 prior-enhanced；两种 prior 都缺失时退回 human-only）：Pearson=`0.8318±0.0008`
  - 两种 mouse prior 都缺失的 `1347` 个基因：real both-prior Pearson=`0.7555±0.0013`，human-only Pearson=`0.7695±0.0006`
  - 两种 mouse prior 都可用的 `11107` 个基因：real both-prior Pearson=`0.8433±0.0010`，human-only Pearson=`0.7456±0.0014`
  - 只有 reconstructed mouse PC1、没有 Saluki mouse PC1 的 `462` 个基因：real both-prior Pearson=`0.7100±0.0074`，human-only Pearson=`0.6800±0.0046`
- 当前解释：
  - 最高性能主要来自有 mouse prior 的基因；没有任何 mouse prior 时，强行使用 prior-enhanced 模型反而略低于 human-only。
  - 实际部署建议采用 coverage-aware fallback：有 mouse prior 时用 prior-enhanced；两种 prior 都缺失时使用 human-only compact sequence/regulatory model。
  - 这个结果应写入应用范围/局限性：模型不要求每个新基因都有 mouse prior，但没有 prior 时不能宣称达到 `r≈0.83` 的 prior-enhanced 水平。

### Study noise versus lambda=0.10 ortholog-informed shrinkage target shift

- 目的：
  - 回答 `lambda=0.10` 轻度混合后的标签虽然与原 human no-Gejman PC1 极高相关，但这点变化在 multi-study experimental noise 中到底是否可忽略。
  - 将 `lambda=0.10` target shift 与不同 human studies 之间的 label difference / residual correlation 放在同一 z-scored label 单位中比较。
- 代码：
  - `src/mrna_half_life_paper/study_noise_vs_orthoreg_shift.py`
- 命令：

```bash
. .venv/bin/activate
PYTHONPATH=paper_pca/src python -m mrna_half_life_paper.study_noise_vs_orthoreg_shift
```

- 新输出：
  - `results/study_noise_vs_orthoreg_shift/lambda0p10_shift_metrics.tsv`
  - `results/study_noise_vs_orthoreg_shift/noise_summary.tsv`
  - `results/study_noise_vs_orthoreg_shift/study_pair_label_metrics.tsv`
  - `results/study_noise_vs_orthoreg_shift/study_residual_vs_reference_metrics.tsv`
  - `results/study_noise_vs_orthoreg_shift/study_residual_pair_metrics.tsv`
  - `figures/study_noise_vs_orthoreg_shift/study_noise_vs_orthoreg_shift.png`
  - `docs/study_noise_vs_orthoreg_shift.md`
- 关键结果：
  - `lambda=0.10` mixed target vs 原 human no-Gejman PC1：N=`10768`，Pearson=`0.9979`，RMSE=`0.0646`，MAE=`0.0501`
  - 排除 `Gejman` 后，study mean labels 两两比较：中位 Pearson=`0.5303`，中位 RMSE=`0.9540`，中位 MAE=`0.7203`
  - normal-study pair 的中位 RMSE 和 MAE 分别约为 `lambda=0.10` target shift 的 `14.8×` 和 `14.4×`
  - study-specific residual 两两相关中位 Pearson=`0.0796`，说明不同 study 的误差结构并不高度一致
  - 更保守地只保留 no-Gejman 且样本数不少于 2 的 7 个 study，study-pair 中位 RMSE=`0.8963`，MAE=`0.6963`，残差相关中位 Pearson=`0.0337`，结论仍成立
  - 预处理敏感性检查：先对 no-Gejman 样本执行完整 Saluki-like preprocessing 后再按 study 取均值，study-pair 中位 Pearson 升至 `0.5937`，但 RMSE=`0.9014`、MAE=`0.6779`，仍远大于 `lambda=0.10` shift
- 当前解释：
  - 可以说 `lambda=0.10` 是轻度 target regularization，而不是重写 human 标签。
  - 更稳妥的论文措辞是“the perturbation introduced by lambda=0.10 is far smaller than normal study-to-study variability and therefore unlikely to alter the dominant human label structure”，而不是简单说“完全可以忽略”。

### BMB Chinese manuscript flow-figure revision and DOCX refresh

- 时间：
  - `2026-07-03 19:42:14 CST`
- 目的：
  - 按 BMB 投稿方向重新整理中文正文和补充材料，使材料与方法、任务边界、使用边界和复现链条更直观。
  - 清理中文稿中的内部讨论口吻，移除 `paper1`、`strongest`、`防守`、`审稿` 等不适合提交稿的表达。
- 新增/重绘图件：
  - `manuscript/bmb_submission/figures/main/Fig04_task_boundary_map.{png,tiff,pdf,svg}`
  - `manuscript/bmb_submission/figures/supplement/FigS01_reproducibility_flow.{png,tiff,pdf,svg}`
  - `manuscript/bmb_submission/figures/supplement/FigS02_prediction_boundary_flow.{png,tiff,pdf,svg}`
- 新增脚本：
  - `scripts/redraw_bmb_flow_figures.py`
- 稿件更新：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 主文新增图 4，用 A-D 子图总结 label QC、human-only benchmark、prior-enhanced benchmark 和 ortholog-informed shrinkage target extension 的输入、目标与解释边界。
  - 补充材料新增补充图 S1 和 S2，分别说明数据到图件的可复现链条，以及不同输入条件下的预测使用边界。
  - 补充材料图号整理为补充图 S1-S11，表号整理为补充表 S1-S10。
  - Code Availability 保留公开仓库链接：`wwzdl/mrna-pc1-label`。
- 验证：
  - Markdown 图片路径均存在。
  - 中文主文：`12` 张图、`12` 个图注；DOCX 内 `12` 个 media 文件。
  - 中文补充材料：`11` 张图、`11` 个图注；DOCX 内 `11` 个 media 文件。
  - DOCX 内容扫描未发现 `paper1|审稿|防守|答审稿|strongest|正式投稿前|三线表源稿|S8b|S8c|S8d|S10b` 等残留内部词。

### BMB Chinese manuscript redundancy cleanup after Fig. 1 / Fig. 4 overlap check

- 时间：
  - `2026-07-05 11:28:00 CST`
- 目的：
  - 删除中文主文中与总体流程图功能重叠的旧任务边界主图，缩短方法部分概念图密度。
  - 保留“任务边界”逻辑，但并回表 1 后的文字说明，避免主文重复解释同一结构。
- 稿件更新：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/README.md`
- 脚本更新：
  - `scripts/redraw_bmb_flow_figures.py`
- 关键修改：
  - 主文删除 `Fig04_task_boundary_map` 及其图注。
  - 原图 5-9 顺次前移为图 4-8，并同步修正正文内所有对应引用。
  - `redraw_bmb_flow_figures.py` 不再生成旧主图 `Fig04_task_boundary_map`，仅保留补充图 `FigS01` 与 `FigS02` 的生成逻辑。
  - “label audit / human-only / prior-enhanced / 0.10 mixed target” 的解释边界改写为表 1 后的一段方法文字。
- 验证：
  - 中文主文图注现顺序为图 1-8，无断号。
  - `rg` 检查未发现 `Fig04_task_boundary_map` 在中文主文中的残留引用。

### BMB Chinese figure filename alignment cleanup

- 时间：
  - `2026-07-05 11:44:00 CST`
- 目的：
  - 使中文主文中的主图文件名与当前图号一致，避免“正文写图 6，文件名仍是 Fig12”这类历史残留。
- 更新内容：
  - `scripts/redraw_bmb_ortholog_combined_cn.py` 输出改为 `Fig05_ortholog_summary_cn`
  - `scripts/plot_bmb_mixed_target_validation.py` 输出改为 `Fig06_mixed_target_validation_cn`
  - `scripts/redraw_bmb_prior_summary_cn.py` 输出改为 `Fig08_prior_summary_cn`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 清理：
  - 删除旧文件 `Fig06_ortholog_summary_cn.*`、`Fig09_prior_summary_cn.*`、`Fig12_mixed_target_validation.*`
- 验证：
  - 中文主文 `8` 个图片引用全部存在。
  - 中文补充材料 `11` 个图片引用全部存在。

### BMB Chinese main-manuscript structure reorder for mixed-target logic

- 时间：
  - `2026-07-04 CST`
- 目的：
  - 让 `0.10` mixed target 的“定义是否仍属 human 标签”先于所有 benchmark 结果出现，避免读者先看到更高 CV，再回头才知道 target 已经变了。
- 稿件更新：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - 方法部分把原 `3.8` 前移为 `3.6`，先交代 mixed target 的构造、边界和验证口径，再进入 `3.7` human-only 与 `3.8` prior-enhanced benchmark。
  - 结果部分把原 `4.7` 拆成两节：
    - `4.3` 先证明 `lambda=0.10` mixed target 仍是 human-dominant target，且 label shift 远小于 normal study-level noise。
    - `4.8` 再报告 mixed target 在 prior-enhanced setting 下的更高可预测性，并明确这不应被写成 pure human sequence rule 的跃迁。
  - 中文主文图号同步重排：mixed-target 验证图前移为 `图 8`，原 pure-human / prior-enhanced 图依次顺延。
- 当前写作判断：
  - 这种顺序更符合审稿人的阅读路径：先确认任务没有被偷换，再看模型性能。
  - `0.10` mixed target 的核心 claim 应拆成两句分别陈述：“它仍主要是 human label” 和 “它在 matched prior-enhanced setting 下更容易被预测”。

### BMB ortholog concordance figure split between main text and supplement

- 时间：
  - `2026-07-05 CST`
- 目的：
  - 避免主文中的 ortholog 柱状图被读成“谁的跨物种相关性最高”的排行榜，突出真正需要支撑的结论：`Gejman` 剔除后，提升在 Pearson/Spearman 和 high-confidence 子集中都保持同方向。
- 新增脚本：
  - `scripts/redraw_bmb_ortholog_summary_figure.py`
- 关键修改：
  - 主文 `figures/main/Fig06_ortholog_bars.*` 现在只保留 4 组：`Full PC1`、`No-Gejman PC1`、`High-conf full`、`High-conf no-Gejman`。
  - 补充图 `figures/supplement/FigS_ortholog_bars.*` 继续保留包含 Saluki 参考对的完整版本。
  - 主文中文/英文图注同步改为“full, pruned, and high-confidence settings”。
  - 英文主文 `Table_3.md` 去掉 `Saluki human PC1 vs Saluki mouse PC1` 行，避免正文表和正文图传递不同重点。
- 当前写作判断：
  - 主文只需要证明 `no-Gejman > full` 在不同指标和子集中都成立。
  - Saluki 参考对更适合放在补充材料或文字说明中，作为口径参照，而不是正文主图中心元素。

### BMB Chinese supplementary deduplication against main text

- 时间：
  - `2026-07-05 12:42:02 CST`
- 目的：
  - 清理中文补充材料中与正文重复的图件和重复解释，让补充更像“参数、控制和扩展证据”而不是第二份正文。
- 更新内容：
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 删除补充材料中与正文完全重复的两张图引用：
    - human PCA 重复图
    - human leave-one-study-out 重复图
  - 补充材料改为：
    - S3 保留 sample-level PCA 的算法口径与 mouse 扩展，不再重复 human PCA 图
    - S4 保留 mouse leave-one-study-out 全排序，不再重复 human 图
  - 补图编号从 `S1-S11` 收紧为连续的 `S1-S9`。
  - 压缩 S2、S9、S10 的开头说明，减少与主文 `3.1`、`4.3`、`4.8` 的重复，保留参数口径、误差尺度和控制实验边界。
- 验证：
  - 中文补充材料图片引用数从 `11` 变为 `9`。
  - `rg` 检查补图编号连续为 `S1-S9`。
  - 所有补图路径存在，`missing=0`。
  - `render_markdown_to_docx.py` 已成功重建 `bmb_supplementary_material_cn.docx`。

### BMB Chinese manuscript source-of-truth consistency check

- 时间：
  - `2026-07-05 12:42:02 CST`
- 目的：
  - 把中文主文/补充作为当前唯一维护主线，并用脚本固定检查图号、图片引用和正文/补充重复图问题。
- 新增脚本：
  - `scripts/check_bmb_cn_consistency.py`
- 检查内容：
  - 中文主文图片引用是否全部存在
  - 中文补充材料图片引用是否全部存在
  - 中文主文图号是否连续为 `图 1..N`
  - 中文补充图号是否连续为 `补充图 S1..SN`
  - 中文主文与中文补充之间是否仍引用 byte-identical 的重复图
- 当前结果：
  - 主文图片引用 `8` 个，补充图片引用 `9` 个
  - 主文图号连续为 `图 1-8`
  - 补图编号连续为 `补充图 S1-S9`
  - 中文链路中不再引用 byte-identical 的正文/补充重复图

### BMB English draft removal under CN-first workflow

- 时间：
  - `2026-07-05 12:42:02 CST`
- 目的：
  - 停止中英文双线同步维护，把 BMB 路线正式收束为中文主文/补充材料唯一主线。
- 删除文件：
  - `manuscript/bmb_submission/bmb_main_manuscript.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_en.md`
  - `manuscript/bmb_submission/bmb_main_manuscript.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_en.md`
  - `manuscript/bmb_submission/bmb_supplementary_material.docx`
  - `manuscript/bmb_submission/bmb_title_page_en.md`
  - `manuscript/bmb_submission/bmb_title_page_en.docx`
  - `manuscript/bmb_submission/bmb_cover_letter_en.md`
  - `manuscript/bmb_submission/bmb_cover_letter_en.docx`
  - `manuscript/bmb_submission/bmb_submission_package.zip`
- 同步更新：
  - `manuscript/bmb_submission/README.md`
  - `manuscript/bmb_submission/bmb_submission_checklist.md`
  - `docs/status_report.md`
  - `scripts/render_bmb_markdown_with_references.py`
- 说明：
  - 英文渲染脚本未删除，但已改为在英文源文件缺失时给出 `CN-first workflow` 的明确提示，避免误判为脚本损坏。
  - 历史审计记录与旧日志中保留对英文稿的引用，作为历史上下文，不再追溯清理。

### BMB Chinese prose densification pass

- 时间：
  - `2026-07-05 14:18:38 CST`
- 目的：
  - 在不改变结果、图号和数值的前提下，压缩中文主文和中文补充材料中 mixed-target 相关段落的解释密度，使成稿更接近投稿阅读节奏。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 主文 `4.3` 保留全部核心数字，但压缩了 mixed-target 合理性、human-dominant 证据和 study-noise 比较的冗余句式。
  - 补充材料 `S9-S10` 改写为“主文已有核心结论，补充只保留边界、机制和精确计算口径”的写法，减少与主文的平行复述。
- 验证：
  - 中文主文与中文补充 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`、补图 `S1-S9`、图片引用完整且无正文/补充重复主图。

### BMB Chinese methods/discussion tightening pass

- 时间：
  - `2026-07-05 14:23:17 CST`
- 目的：
  - 继续把中文主文从“技术记录风格”收向“投稿成稿风格”，重点压缩方法前半段和讨论部分的句子密度。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - `3.1-3.3` 改写为更紧凑的方法描述，保留 `MOESM2`/`MOESM3` 分工、PCA 诊断目的和共识标签重建流程，但减少开发式解释。
  - `5. 讨论` 压缩为更清晰的投稿口吻：先讲贡献边界，再讲 prior-enhanced 与 `MOESM3` 控制，再讲 mixed target 与使用边界，最后给出方法学延展方向。
- 验证：
  - 中文主文 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整。

### BMB Chinese ending-section tightening pass

- 时间：
  - `2026-07-05 14:26:57 CST`
- 目的：
  - 继续压缩中文主文结尾部分的句子密度，使“局限性”和“数据与代码可获得性”更接近正式投稿口吻。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - `6. 局限性` 从单段六连句改为三层结构：数据外推边界、任务/标签边界、输入与部署边界。
  - `8. 数据与代码可获得性` 改为更标准的投稿写法，直接交代公开数据来源与版本化仓库入口。
- 验证：
  - 中文主文 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese headline-consistency pass

- 时间：
  - `2026-07-05 14:31:17 CST`
- 目的：
  - 统一中文主文摘要、引言开头与结论的主线表述，明确区分标签诊断、human-only benchmark、prior-enhanced benchmark 与 `0.10` 混合 target 扩展。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - 摘要改为“三条任务线”口径，避免把标签诊断、输入增强和混合 target 结果写成同一类结论。
  - 引言首段中的 benchmark 动机改为更正式的“标签质量先于性能比较”表述，并补强 mixed-target extension 的语义边界。
  - 结论重写为三条边界清楚的结果陈述，对应 human-only 基线、prior-enhanced 设置和 `0.10` mixed target。
- 验证：
  - 中文主文 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese methods boundary-clarification pass

- 时间：
  - `2026-07-05 14:35:26 CST`
- 目的：
  - 继续打磨中文主文 `3.6-3.9` 的方法叙述，进一步区分“输入端加入 mouse prior”和“标签端进行 `0.10` 轻度混合”这两类实验。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 主文 `3.6` 明确写成“改标签，不改输入”，并压缩 mixed target 的构造与验证口径。
  - 主文 `3.7-3.8` 明确对照为“human-only 基线”和“目标不变、只改输入的 prior-enhanced 设定”。
  - 主文 `3.9` 把 permutation、residual analysis 和 no-direct-prior ablation 分别对应到三个独立问题，减少控制实验定义上的歧义。
  - 补充材料同步补入两句边界提示，避免审稿人只看补充时把 mixed target 与 prior-enhanced benchmark 误读为同一类实验。
- 验证：
  - 中文主文与补充 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese results-interpretation tightening pass

- 时间：
  - `2026-07-05 14:39:36 CST`
- 目的：
  - 收紧中文主文 `4.6-4.8` 及对应补充解释的表述，使“prior 不是靠多几列取胜”“prior 很强但不是简单复制”“mixed target 更易预测但不等于 pure human sequence rule 跃迁”三层结论更易读。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 主文 `4.6` 改为先定义 permutation control 排除的疑问，再解释“被打乱的是 prior 数值与 gene identity 的对应关系”。
  - 主文 `4.7` 改写为两阶段 residual decomposition 口径，突出 `prior prediction + human-feature correction` 的解释。
  - 主文 `4.8` 压缩 mixed-target 结果段，明确 `0.855` 属于 mixed-target prior-enhanced 任务，不应表述为 pure human 学习能力跃迁。
  - 补充材料 `S7-S9` 同步统一语气，使其与主文使用相同的边界表述。
- 验证：
  - 中文主文与补充 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese discussion-limitation closing pass

- 时间：
  - `2026-07-05 14:43:19 CST`
- 目的：
  - 收紧中文主文 `5.讨论` 与 `6.局限性` 的结尾口吻，使整篇稿件在末尾更集中地落到“贡献边界、任务边界和使用边界”。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - `5.讨论` 改写为四段式收束：三层问题拆分、prior-enhanced 的方法学解释、`MOESM3` 与 mixed target 的定位、以及实际部署边界与后续方法学延展。
  - `6.局限性` 压缩为三个清楚的限制来源：数据外推、任务/标签边界、模型与输入条件边界。
  - 补充材料 `S15` 删除已不再适用的“中英文稿统一排版”说法，改为更中性的正式投稿版表述。
- 验证：
  - 中文主文与补充 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese terminology-harmonization pass

- 时间：
  - `2026-07-05 14:55:18 CST`
- 目的：
  - 对中文主文与补充做投稿前总清稿，统一标题、摘要、结论、图注和表题中的核心术语，减少代码式字段名和双重口径残留。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.md`
  - `manuscript/bmb_submission/bmb_supplementary_material_cn.docx`
- 关键修改：
  - 统一 `0.10` 路线的英文短语为 `90:10 mixed mammalian stability target`，并把残留的 `score`/`mixed-label extension` 改为 `target`/`mixed-target extension`。
  - 把正文图 7 和相关结果中的 `human_pruned_pc1` 改为更面向读者的 `human no-Gejman PC1`。
  - 将多处泛称 `Saluki PC1` 明确为 `Saluki human PC1` 或 `Saluki 发布 PC1 标签`，并把 `augmented-information` 收口为 `cross-species transfer setting`。
  - 关键词同步改为更贴近正文口径的 `cross-species prior` 与 `prior-enhanced benchmark`。
- 验证：
  - 中文主文与补充 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese reviewer-first-read readability pass

- 时间：
  - `2026-07-05 15:01:20 CST`
- 目的：
  - 按“审稿人第一遍扫读”口径收紧中文主文的可读性，只处理最影响首轮理解的长句和别扭表述，不改实验结论与数值。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - 摘要改为更短的句群结构，减少一段内连续堆叠三层任务、结果和边界的阅读负担。
  - `3. 材料与方法` 开头补出“前两条路线改输入、第三条路线改目标”的对照句，帮助读者快速区分任务边界。
  - 将正文中两处读起来重复的 `mammalian stability target 的 target extension` / `mixed-target extension 或 ... 候选目标` 改为更自然的投稿表述。
- 验证：
  - 中文主文 DOCX 已重建。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无重复主图。

### BMB Chinese ancillary submission-material pass

- 时间：
  - `2026-07-05 15:18:40 CST`
- 目的：
  - 在保持“中文主文与中文补充为 source of truth”的前提下，补齐 BMB 中文内部外围投稿材料，给后续英文 `title page`、`cover letter`、`Author contributions` 和系统表单提供统一口径。
- 更新文件：
  - `manuscript/bmb_submission/bmb_title_page_cn.md`
  - `manuscript/bmb_submission/bmb_title_page_cn.docx`
  - `manuscript/bmb_submission/bmb_cover_letter_cn.md`
  - `manuscript/bmb_submission/bmb_cover_letter_cn.docx`
  - `manuscript/bmb_submission/bmb_editor_highlights_cn.md`
  - `manuscript/bmb_submission/bmb_editor_highlights_cn.docx`
  - `manuscript/bmb_submission/bmb_author_contributions_cn.md`
  - `manuscript/bmb_submission/bmb_author_contributions_cn.docx`
  - `manuscript/bmb_submission/bmb_statements_and_declarations_cn.md`
  - `manuscript/bmb_submission/bmb_statements_and_declarations_cn.docx`
  - `manuscript/bmb_submission/README.md`
  - `manuscript/bmb_submission/bmb_submission_checklist.md`
- 关键修改：
  - 新增中文内部 `title page`，统一题目、作者、单位、通讯作者、funding、关键词和 ORCID 状态。
  - 新增中文内部 `cover letter`，按 BMB `Original research article` 口径概括三条主贡献，并明确 human-only、prior-enhanced 与 mixed-target extension 的 claim boundary。
  - 新增 `editor highlights`，专门整理“编辑最该快速看到的点”和“不建议过度声称的点”，供后续英文附信和系统备注复用。
  - 新增 `author contributions` 草案，给出一个可直接落入主稿的简版声明，以及一个带 `待确认` 标记的 CRediT 细化表。
  - 新增 `Statements and Declarations` 汇总，统一 funding、data/code availability、competing interests、ethics/consent 和 ORCID 口径，并把生成式 AI 披露写成条件性说明而非自动并入主稿。
  - README 与 checklist 同步更新，明确这些外围件均为中文内部版，英文版本仍待中文定稿后统一生成。
- 验证：
  - 五份中文外围材料均已成功导出为 DOCX。
  - 新增 DOCX 文件已在 `manuscript/bmb_submission/` 目录下就位，可直接内部审阅。

### BMB Chinese figure simplification and narrative-tightening pass

- 时间：
  - `2026-07-05 16:30:29 CST`
- 目的：
  - 参考外部审稿意见与 `BMB_mrna_rewrite_cn_0706.docx` 的示意风格，对中文主文做一轮更接近 BMB 审稿阅读习惯的收敛：减少防御式叙述，压缩方法主图复杂度，并把主文图形从“复现说明图”往“主结论图”推进。
- 更新文件：
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.md`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
  - `scripts/redraw_bmb_workflow_figure.py`
  - `scripts/redraw_bmb_study_influence_cn.py`
  - `scripts/plot_bmb_mixed_target_validation.py`
  - `scripts/redraw_bmb_sequence_progression_cn.py`
  - `scripts/redraw_bmb_prior_summary_cn.py`
- 关键修改：
  - 重写摘要、方法 `3.6-3.9`、结果 `4.3/4.6/4.7/4.8`、讨论、局限性和结论中的多处“防误读式”表达，改为更短、更直接的结果-解释-边界结构。
  - 将主图 `Fig 1` 压缩为三段式框架：`Data foundation -> Study-aware label audit -> Benchmark tasks`，并单独给出一条 `Claim boundary`，减少原流程图中过密的小模块。
  - 将 `Fig 4` 重画为横向 `leave-one-study-out` 结果图：左图只保留 `full-label geometry`，右图只保留相对 `Saluki human PC1` 的 `delta Pearson`，突出 `Gejman` 的唯一正向增益。
  - 简化 `Fig 6` 的 mixed-target 对照，只保留 mixed-vs-human、study-study median 和 single-study-vs-human median 三类同尺度比较。
  - 简化 `Fig 7` 的 human-only 主图结构，删除容易引发“融合”误读的额外面板，保留 model progression、feature ablation、label-definition sensitivity 和 subset sensitivity。
  - 重排 `Fig 8` 的 prior-enhanced 概览，并把 prior coverage 改为更紧凑的 partition 展示。
  - 同步更新图 1 与图 4 图注，使其与新版图的结构和信息量一致。
- 验证：
  - 主图 `Fig01/Fig04/Fig06/Fig07/Fig08` 已全部重生成。
  - 中文主文 DOCX 已重新导出为 `bmb_main_manuscript_cn.docx`。
  - `check_bmb_cn_consistency.py` 通过：主文 `图 1-8`，补图 `S1-S9`，中文链路图片引用完整且无主补重复图。

### BMB Fig. 1 office-style redraw with editable PPTX

- 时间：
  - `2026-07-05 16:53:00 CST`
- 目的：
  - 参照外部 `pro` 版本 Fig. 1 的视觉组织方式，把正文主图 1 重画为更接近 WPS/PowerPoint 风格的流程图，同时保留本项目自己的术语、任务边界和数值，并产出可编辑版本。
- 更新文件：
  - `scripts/redraw_bmb_workflow_figure.py`
  - `manuscript/bmb_submission/figures/main/Fig01_workflow.png`
  - `manuscript/bmb_submission/figures/main/Fig01_workflow.tiff`
  - `manuscript/bmb_submission/figures/main/Fig01_workflow.pdf`
  - `manuscript/bmb_submission/figures/main/Fig01_workflow.svg`
  - `manuscript/bmb_submission/figures/main/Fig01_workflow_editable.pptx`
  - `manuscript/bmb_submission/bmb_main_manuscript_cn.docx`
- 关键修改：
  - 用 `PIL` 重画主图 1，采用三段式 office 风格布局：`A Data and label audit`、`B Human-only benchmark`、`C Cross-species extension`。
  - 在 `A` 区显式保留 `MOESM2`、`MOESM3`、`Ensembl orthologs` 与 `study-aware label reconstruction pipeline`，并用单独注释框强调 `Gejman` 是主异常 study。
  - 在 `B` 和 `C` 区按流程图方式分别展示 human-only、prior-enhanced 和 weakly mixed target 路线，并保留 `OOF Pearson ~0.74 / ~0.83 / ~0.855` 作为读图级摘要。
  - 新增基于 `python-pptx` 的 `Fig01_workflow_editable.pptx`，可直接在 WPS 演示或 PowerPoint 中打开继续微调。
- 验证：
  - 新版 `Fig01_workflow.png` 与 `Fig01_workflow_editable.pptx` 已生成。
  - 中文主文 DOCX 已重导，确保新图 1 已嵌入当前主稿。
  - `check_bmb_cn_consistency.py` 再次通过。

### BMB main-docx cleanup

- 时间：
  - `2026-07-05 17:04:13 CST`
- 目的：
  - 简化中文主稿文件管理，只保留一个当前 Word 主文件。
- 更新文件：
  - `manuscript/bmb_submission/README.md`
  - 删除并行 Word 副本
- 关键修改：
  - 当前目录说明改为只维护 `bmb_main_manuscript_cn.docx` 作为中文主文 Word 输出。
  - 移除并行副本口径，当前只维护 `bmb_main_manuscript_cn.docx`。

## 历史里程碑摘要

以下内容整理自 `docs/codex_handoff.md`，用于快速接续，不代表原始执行日期：

- 完成 human/mouse consensus label 重建
- 完成 study influence 分析并识别 human `Gejman` outlier
- 完成 ortholog consistency 分析
- 完成 pure human compact/regulatory baseline
- 完成 Saluki hand-crafted feature package 快扫
- 完成 ortholog 子集 prior benchmark
- 完成 global common-gene-universe prior benchmark
- 首轮 TCBB 投稿包已生成

### BMB final reviewer-risk, reproducibility, and release pass

- 时间：
  - `2026-07-13 CST`
- 目的：
  - 将外部预审意见中有证据依据的部分落实到中英文主文、补充材料、脚本和发布边界，同时保留 `lambda=0.10` 作为独立核心标签构建贡献。
- 关键修改：
  - 主线收敛为三项可独立检验的命题：Saluki-label-independent study influence audit、固定 Saluki human PC1 下的 human-only/cross-species transfer，以及 `0.10 ortholog-informed shrinkage target`。
  - 明确 `Gejman` 由不使用 Saluki label 或下游分数的 PC1-stability 排序恢复；因缺少 size-matched null，正文只称其为 dominant influence，不再称为样本数校正后的 statistical outlier。
  - 明确两类 mouse priors 均为 human CV 之前固定的 mouse-side external covariates；fold-wise permutation 打断 gene-prior identity 后增益消失，结果属于 cross-species transfer 而非 sequence-only prediction。
  - `lambda=0.10` 未降级为附属结果，而是用 label geometry、study-noise comparison 和 human-only cross-target evaluation 建立边界；`0.855` 只用于 regularized-target estimand，不与固定 Saluki human PC1 排名混用。
  - 中英文参考文献改为逐条独立段落；删除活跃目录中的旧 mixed-target 图别名、Pro 样图副本和过时 ortholog bar 脚本。
  - 新增 `scripts/reproduce_bmb_key_results.sh` 分阶段复现入口，README 补入 Saluki Zenodo DOI、18.8-GB datapack 放置结构和任务使用边界。
  - 调整 Git ignore 规则，使 compact result tables/summaries 进入 GitHub 发布范围，同时继续排除大矩阵、逐基因 OOF predictions、历史 smoke/tuning outputs 和 TIFF 投稿资产。
- 验证：
  - BMB package audit：`99 checks, 0 failures, 0 warnings`。
  - Python compileall 与环境 smoke test 全部通过。
  - 英文主文 PDF 视觉核对为 23 页，英文 Online Resource 1 为 18 页；在用主图 Fig. 1-8 与补图 S1-S4 无越界、重叠或不可读文字。
  - 最终上传包 `bmb_upload_candidate_20260713.zip` 通过 `unzip -t`，SHA256 为 `0cee9bc694c7bcc60899ca8d5d9f852db87263e953820a1b7f7c62317c066c1b`。

### BMB Supplementary Figure S1 connector geometry cleanup

- 时间：
  - `2026-07-14 CST`
- 目的：
  - 修复 prediction-use boundary 流程图中 panel c 与 panel d 的框位和箭头几何，使输入关系在正文缩放尺寸下仍清楚、对称。
- 关键修改：
  - panel c 的两个左侧输入框整体下移，两条输入箭头改为等长水平箭头，并围绕 transfer model 对称排布。
  - panel d 的 `Regulatory blocks?` 判断框下移，延长上方纵向 `no` 分支；三个右向分支采用相同水平跨度，下方 `yes/no` 箭头平行排列。
  - 进一步扩大 panel c 两个输入框及 panel d 右侧 compact/base 输出框之间的可见纵向间隔，避免圆角框外扩后产生视觉粘连。
  - `redraw_bmb_flow_figures.py` 新增 `BMB_SKIP_PDF=1` 开关；本轮仅刷新 SVG、PNG 和 TIFF，未更新 PDF。

### BMB Supplementary Figure S4 precision clarification

- 时间：
  - `2026-07-14 CST`
- 关键修改：
  - 确认 tuning-scan 柱长使用结果表中的未舍入 Pearson，而旧图柱端仅显示三位小数，导致近似值看似相等但柱长略有差异。
  - 将柱端标签改为五位小数，并在中英文图注中明确柱长使用未舍入值；本轮未更新 PDF。

### BMB DOCX table-width and pagination cleanup

- 时间：
  - `2026-07-14 CST`
- 关键修改：
  - 在统一 Markdown-to-DOCX 渲染脚本中加入内容感知列宽：叙述列按内容长度加宽，纯数值列收窄并居中，表格固定为版心宽度。
  - 4–8 列表保留竖版；9 列及以上补充表自动使用单页横版，避免长字段挤压和过小字号。
  - 宽表中的下划线标识符和路径可在分隔符后软换行；数值单元格禁止拆行，表头跨页重复，数据行不跨页拆分。
  - 本轮重新生成中英文正文与补充材料 DOCX；未生成 PDF。

### BMB final pre-submission acceptance audit

- 时间：
  - `2026-07-14 CST`
- 审核范围：
  - 中英文正文与补充材料、主图 Fig. 1-8、补图 S1-S4、表格、引文与参考文献、BMB 投稿附件、复现入口、OOF 预测完整性和 GitHub 模拟发布范围。
- 关键处理：
  - 对照 BMB 当前官方指南复核摘要、关键词、作者-年份引文、无编号字母顺序参考文献、连续行号、补充材料 metadata 和 600 dpi 组合图要求。
  - 按原始分辨率逐张检查 12 张活跃图件，未见裁切、重叠、panel label 缺失、数值越界或图注不匹配。
  - 将保留结果快照中的旧 `*_vs_author` / `mouse_author` 列名统一为当前代码使用的 `*_vs_saluki` / `saluki_mouse_prior`，不改变数值。
  - GitHub 模拟发布范围收紧为 277 个文件、24.60 MiB，排除 raw/interim 大数据、大矩阵、逐基因 OOF、DOCX、TIFF、PDF 和 ZIP；未发现绝对路径、密钥或 token 模式。
- 验证：
  - DOCX-only BMB preflight：`136 checks, 0 failures, 0 warnings`。
  - OOF 完整性：9 个 global fixed-target、3 个 cross-target、18 个 ortholog-informed shrinkage 文件以及 11,107-gene residual 子集全部通过。
  - Python compileall、全部 shell 语法、`pip check`、项目 smoke test、MOESM2/MOESM3/Ensembl 115 SHA-256 均通过。
  - 本轮未重跑全部耗时 GPU 训练；校验的是已有 OOF 结果、输入、源码、结果表与复现入口。
- 冻结边界：
  - 最终英文 Online Resource 1 PDF、上传候选目录与 ZIP 继续按作者要求留到 ORCID、Author Contributions、GitHub release 和 license 冻结后一次生成。
  - 完整项目备份路径固定为 `../workspace_archive/mrna-pc1-label_paper_pca_presubmission_20260714.tar.zst`，SHA-256 保存在同名边车文件。

### BMB Figure 4 stability-panel simplification

- 时间：
  - `2026-07-14 CST`
- 关键修改：
  - 删除 Fig. 4a 中从每个相关系数点延伸到 `r = 1` 的水平线，改为排序散点图，仅保留 `r = 1` 的虚线参照。
  - 中英文图注同步说明圆点表示 leave-one-study-out PC1 与完整 human PC1 的相关系数，消除对水平线含义的歧义。
  - 绘图脚本支持 `BMB_SKIP_PDF=1`；本轮只更新 PNG、SVG 和 600 dpi lossless TIFF，并重建中英文 DOCX。
- 验证：
  - 原始分辨率视觉检查通过，panel label、数值标注、坐标轴与图题未见重叠或越界。
  - BMB package audit：`136 checks, 0 failures, 0 warnings`。

### GitHub repository-address migration

- 时间：
  - `2026-07-15 CST`
- 关键修改：
  - 将稿件、补充材料、title page、cover letter、declarations、README、`CITATION.cff`、投稿检查表和自动审计规则中的仓库地址统一改为 `https://github.com/wwzdl/mrna-pc1-label`。
  - 删除活跃文本源中所有原仓库地址，并重建中英文 DOCX；未生成 PDF 或 ZIP。
- 验证：
  - BMB package audit：`136 checks, 0 failures, 0 warnings`。
  - 2026-07-15 在公网访问新地址仍返回 404，因此 Data/Code Availability 继续使用将来时，待仓库公开后再改为当前时态。

### Public GitHub release

- 时间：
  - `2026-07-15 CST`
- 关键处理：
  - 将 `paper_pca` 初始化为独立 Git 仓库，并按 `.gitignore` 发布源码、紧凑结果表、稿件 Markdown、参考文献及 PNG/SVG 图件。
  - 公开仓库为 `https://github.com/wwzdl/mrna-pc1-label`；raw/interim 大数据、大矩阵、逐基因 OOF arrays、feature packages、DOCX、PDF、TIFF 和 ZIP 未进入公开提交。
  - 中英文 Data/Code Availability、title page、cover letter、declarations、补充材料复现说明及 README 均由将来时改为当前时态。
  - `CITATION.cff` 增加 `version: 1.0.0` 与发布日期；稿件审阅状态标记为 `mRNA-PC1-label-v1.0`。
- 验证：
  - 初始公开提交包含 `277` 个文件，范围约 `24.60 MiB`；未发现本机绝对路径、私钥或 GitHub token。
  - 发布后再次执行中英文一致性和 BMB package audit；最终结果记录于本次发布提交。

### Sample-count-aware study audit and reviewer-driven manuscript revision

- 时间：
  - `2026-07-15 CST`
- 目的：
  - 回应预审中最关键的质疑：`Gejman` 有 15/54 个 human samples，其 leave-one-study-out 几何影响是否主要由样本数、动态基因集合或预处理选择造成。
- 新增分析：
  - 500 次 size-matched random removal；每次从 39 个非 Gejman 样本中无放回移除 15 个，并完整重跑 coverage、z-score、PCA imputation、quantile normalization 和 PC1 extraction。
  - Dynamic 与 fixed gene universe 的 leave-one-study-out 对照。
  - Coverage threshold 3/5/10、PCA rank 3/5/10、iterative-PCA 与 median imputation 的 12 组预处理敏感性分析。
  - Equal-study PCA weighting 与 study-mean collapse 两种 study-balanced estimators。
- 结果与解释：
  - Gejman observed PC1 stability 为 `0.9550`；size-matched 零分布中位数 `0.9438`，95% interval `0.9319-0.9561`，低尾 `p=0.964`。因此不能称为 sample-count-adjusted geometric outlier。
  - Gejman 的 Saluki agreement gain 为 `+0.0314`，ortholog gain 为 `+0.0147`；两者在 size-matched 零分布中的高尾经验检验均为 `p=0.002`。
  - 12 个预处理/基因宇宙设定均将 Gejman 排第 1；两种 study-balanced estimators 均将其排第 3。最终主张收紧为“样本数杠杆叠加有方向的外部一致性改善”。
- 稿件修改：
  - Results 按 label audit、human-only fixed target、cross-species transfer、ortholog-informed target shrinkage 重排。
  - 新增 measurement-error shrinkage model、inner-validation early-stopping 说明、扩展 Discussion 与 Limitations。
  - 将 `Saluki-label-independent` 改为 `Saluki-label-independent`，将读者-facing 的 `ortholog-informed target shrinkage` 改为 `ortholog-informed target shrinkage`；代码和结果字段名保持不变。
  - 以新的补充图 S4 和补充表 S4 替换信息价值较低的 XGBoost tuning scan。

### ChatGPT 5.6 Sol review adjudication and final evidence-role cleanup

- 时间：
  - `2026-07-15 CST`
- 采纳的意见：
  - 将 `Gejman` 的 15/54 样本量混杂作为正文核心结果处理，而不是只放在 Limitations；正文和摘要同步报告 size-matched null 与两类 study-balanced estimators。
  - 同时报告 dynamic/fixed gene universe、coverage、PCA rank 和 imputation sensitivity；明确 outer test fold 从未用于 early stopping。
  - Results 按 label audit、fixed-target human-only、cross-species transfer、target shrinkage 重排；摘要定义 PC1 并只保留不同 estimand 下不会混淆的核心数值。
  - 将 `lambda=0.10` 写成 measurement-error working model 下的 shrinkage estimator，并显式给出 variance、species-shift bias 与 MSE 的关系。
  - 将 Saluki agreement 限定为恢复已发布处理选择；将 ortholog concordance 写成作用不同的跨物种比较轴，不再表述为完全独立的外部验证。
  - Fig. 1 中 fixed-target 数值统一为同一 repeated 10-fold 口径：human-only `0.748`、real priors `0.830`、shuffled priors `0.748`；Fig. 2 与补图统一使用 target-shrinkage 和 concordance 术语。
- 未机械采纳的意见：
  - 不把 study-balanced 分析包装成 `Gejman` 仍为第 1；真实结果为第 3，正文据此收紧结论。
  - 不虚构第二个独立 half-life compendium、模拟结果或理论最优 `lambda`；这些内容保留为明确局限和后续验证方向。
  - 不将 `lambda=0.10` 降为无关附属分析；它保留为第三项核心贡献，但与固定 Saluki human PC1 的模型排名完全分开。

### Submission-author metadata revision

- 时间：
  - `2026-07-16 CST`
- 作者与单位：
  - 活跃投稿材料的作者统一为 Wenzhuo Wang 与 Ying Shao；Ying Shao 标记为唯一通讯作者。
  - 两位作者均使用 `School of Science, Dalian Maritime University, Dalian 116026, China`。
  - Ying Shao 的通讯作者邮箱按作者要求暂留空，待最终投稿前补充；两位作者 ORCID 均待确认。
- 声明与贡献：
  - 删除原国家自然科学基金编号，Funding 统一为 `No specific funding was received for this work.`。
  - Author Contributions 改为两作者口径：Wenzhuo Wang 负责分析、软件、数据、验证、图表和初稿；Ying Shao 负责构思、方法、验证、监督、项目管理和审阅修改。正式投稿前仍需两位作者确认其与实际贡献一致。
- 同步范围：
  - 中英文正文、补充材料、title page、cover letter、declarations、中文贡献核对稿、投稿清单、README、`CITATION.cff` 与自动审计规则。
  - 由于不可移动既有 v1.1 标签，两作者稿件与元数据冻结为新标签 `mRNA-PC1-label-v1.2`。

### Ying Shao ORCID confirmation

- 时间：
  - `2026-07-16 CST`
- 更新：
  - Ying Shao 的 ORCID 确认为 `0000-0002-4056-5757`。
  - 已同步中英文正文、补充材料、title page、declarations、中文 cover letter、投稿清单、仓库说明和 `CITATION.cff`。
  - 自动投稿审计新增该 ORCID 在英文主稿、补充材料、title page、declarations 与 `CITATION.cff` 中的一致性检查。
  - Wenzhuo Wang 的 ORCID 与 Ying Shao 的通讯作者邮箱仍待最终确认。

### Ying Shao corresponding-author email confirmation

- 时间：
  - `2026-07-16 CST`
- 更新：
  - Ying Shao 的通讯作者邮箱确认为 `yshao@dlmu.edu.cn`，ORCID 保持为 `0000-0002-4056-5757`。
  - 已同步中英文正文、补充材料、title page、cover letter、declarations、投稿清单和 `CITATION.cff`。
  - 自动投稿审计改为要求英文活动稿件与仓库引用元数据中的邮箱精确匹配，不再允许通讯邮箱留空。
  - 作者元数据目前仅剩 Wenzhuo Wang 的 ORCID 待确认。

### Native Word equation conversion

- 时间：
  - `2026-07-16 CST`
- 稿件更新：
  - 将 study influence、经验 p 值、0.10 target shrinkage、measurement-error decomposition 和 residual variance 的核心定义整理为正文公式（1）-（9）。
  - 中英文补充材料中的 iterative-PCA convergence criterion 与 remaining-variance definition 分别编为补充式（S1）和（S2）；短数学表达使用 `$...$` 行内标记，代码字段继续保留等宽字体。
- 渲染与验证：
  - `render_markdown_to_docx.py` 新增 LaTeX -> MathML -> OMML 转换，生成 Word/WPS 可编辑公式对象及右侧编号。
  - `render_markdown_to_pdf.py` 使用同一 LaTeX 源生成矢量 SVG 公式，避免 WeasyPrint 对 MathML 分数、上下标的降级。
  - 新增公式渲染单元测试和自动审计：英文/中文正文分别检测到 44/46 个 OMML math objects，其中各 9 个编号公式；英文/中文补充分别检测到 22/24 个 OMML math objects，其中各 2 个编号公式。
  - DOCX-only preflight 为 `158 checks, 0 failures, 1 warning`；warning 仅来自打开文档产生的 WPS 临时锁文件。
