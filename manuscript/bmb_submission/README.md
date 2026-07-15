# BMB Submission Package README

本目录是 `paper1` 面向 *Bulletin of Mathematical Biology* 的投稿材料工作目录。

## 当前定位

- 目标期刊：*Bulletin of Mathematical Biology*（Springer, hybrid）
- 当前稿件口径：`Original research article`
- 发表模式默认假设：`subscription` 路线，不强制 APC
- 当前维护策略：中文主文和中文补充材料保留为内部源稿；英文主文、英文补充材料、英文 title page、英文 cover letter 和英文 declarations 已形成当前投稿工作稿

## 目录内容

- `bmb_main_manuscript_cn.md`
  - 中文主文源文件；当前主线稿
- `bmb_main_manuscript_cn.docx`
  - 中文主文 Word 导出稿
- `bmb_main_manuscript_en.md`
  - 英文主文源文件；按当前中文主文同步生成
- `bmb_main_manuscript_en.docx`
  - 英文主文 Word 导出稿
- `bmb_supplementary_material_cn.md`
  - 中文补充材料源文件；当前主线稿
- `bmb_supplementary_material_cn.docx`
  - 中文补充材料 Word 导出稿
- `bmb_supplementary_material_en.md`
  - 英文补充材料源文件；按当前中文补充同步生成
- `bmb_supplementary_material_en.docx`
  - 英文补充材料 Word 导出稿
- `bmb_title_page_cn.md`
  - 中文内部标题页源稿
- `bmb_title_page_cn.docx`
  - 中文内部标题页 Word 导出稿
- `bmb_cover_letter_cn.md`
  - 中文内部附信源稿
- `bmb_cover_letter_cn.docx`
  - 中文内部附信 Word 导出稿
- `bmb_title_page_en.md`
  - 英文 title page 源稿
- `bmb_title_page_en.docx`
  - 英文 title page Word 导出稿
- `bmb_cover_letter_en.md`
  - 英文 cover letter 源稿
- `bmb_cover_letter_en.docx`
  - 英文 cover letter Word 导出稿
- `bmb_editor_highlights_cn.md`
  - 中文内部投稿要点摘要；用于后续英文 cover letter 和系统备注统一口径
- `bmb_editor_highlights_cn.docx`
  - 中文内部投稿要点摘要 Word 导出稿
- `bmb_author_contributions_cn.md`
  - 中文作者贡献草案；需作者最终确认
- `bmb_author_contributions_cn.docx`
  - 中文作者贡献草案 Word 导出稿
- `bmb_statements_and_declarations_cn.md`
  - 中文声明汇总；用于后续英文 `Statements and Declarations` 重建
- `bmb_statements_and_declarations_cn.docx`
  - 中文声明汇总 Word 导出稿
- `bmb_statements_and_declarations_en.md`
  - 英文 `Statements and Declarations` 源稿
- `bmb_statements_and_declarations_en.docx`
  - 英文 `Statements and Declarations` Word 导出稿
- `bmb_submission_requirements.md`
  - 官方投稿要求摘要与链接
- `bmb_submission_checklist.md`
  - 当前执行清单
- `upload_candidate_manifest.md`
  - 最终上传候选目录的文件清单模板
- `upload_candidate/` 与 `bmb_upload_candidate_*.zip`
  - 仅在稿件、补充 PDF 与作者信息全部冻结后生成；当前编辑阶段不保留活动副本
- `references/`
  - 与正文/补充实际引用一致的 BibTeX、RIS 参考文献库及审计说明
- `figures/`
  - 当前 PNG/SVG 图件及本地 600 dpi TIFF 投稿图件；表格以 Markdown 和结果 TSV 为可编辑源

当前目录已经保留英文主文、英文补充材料、英文 title page、英文 cover letter 和英文 declarations。最终补充材料上传名将使用 `Online_Resource_1_supplementary_material.pdf`；当前 PDF、上传候选目录和 ZIP 按作者要求推迟到版本冻结后一次生成。生成式 AI 使用说明已按 BMB 要求写入方法部分。GitHub 公开版本已完成；最终系统上传前仍需作者侧确认两位作者的 ORCID、Author Contributions 与软件许可。

## 最新内容更新

- 中文主文 Word 稿现只保留一个文件：`bmb_main_manuscript_cn.docx`。
- 英文主文、英文补充材料、英文 title page、英文 cover letter 和英文 declarations 已同步写入并导出为 `docx`。
- 主文 Figure 1 与补充 Figure S1 均由 `redraw_bmb_flow_figures.py` 使用 matplotlib 原生矢量对象生成；活跃投稿链路仅引用脚本生成的图件。
- 活跃图件当前保留 SVG、600 dpi PNG 和 600 dpi lossless TIFF；PDF 由同一脚本在最终冻结阶段重建，图内数值与术语受脚本和结果表统一控制。
- 中文补充材料已去掉与正文完全重复的 human PCA、human leave-one-study-out，以及重复度过高的 Saluki feature quick scan / compact ablation / OOF blend / subset scan 图；当前中文补图编号收紧为连续的 `S1-S4`，并把 human 主证据保留在正文、mouse 扩展与必要控制保留在补充中。
- 一致性检查脚本 `scripts/check_bmb_cn_consistency.py` 已支持 `--lang cn|en|both`，可同时检查中英文主文/补充图号连续性、图片存在性，以及 active submission markdown 中是否残留退役标签术语。
- `scripts/render_bmb_submission_bundle.sh` 可一次性重建中英文正文、补充材料和外围投稿材料的 `docx` 导出物，减少手工漏导。
- `scripts/audit_bmb_submission_package.py` 已加入 preflight，可检查英文摘要词数、关键词数、声明、补充材料 metadata、DOCX 行号/页码、引用残留、退役术语、Online Resource 口径、参考文献引用一致性和 600 dpi 图件。
- `scripts/assemble_bmb_upload_candidate.sh` 只重建 `upload_candidate/`，避免每次预检都重复生成 zip；压缩包只在最终发布时生成一次。
- `BMB_SKIP_PDF=1 bash scripts/preflight_bmb_release.sh` 可在编辑阶段重建全部 DOCX 并完成审计，同时跳过 PDF 与上传候选目录，避免混合新旧版本。
- `scripts/audit_oof_integrity.py` 检查关键 OOF 文件的规范基因数、唯一 ID、共同 universe、fold 覆盖和预测缺失。
- 已补齐中文内部外围材料：标题页、附信、投稿要点摘要、作者贡献草案和声明汇总，供后续英文 `title page`、`cover letter` 和系统表单统一转换。
- 主文已纳入 10-fold global prior robustness、ortholog-regularized target 的 label geometry、cross-target evaluation、prior ablation 与 no-direct-prior 解释边界。
- 补充材料保留 ortholog-regularized target 的相关补充表与必要控制结果。
- 旧 TCBB 打包脚本和历史备份目录不再属于当前活跃 BMB 投稿链路；当前项目目录已作为公开 `mrna-pc1-label` 仓库根目录。

## 已对齐的关键官方要求

- 单盲评审
- abstract `150-250` 词
- `4-6` 个 keywords
- 文中参考文献改为作者-年份制
- 英文主文参考文献均在正文中被引用，且按第一作者字母顺序排列
- 主稿与补充材料分开
- figures 保留 `PNG/SVG/TIFF` 高分辨率版本；最终 PDF 尚未冻结
- data availability statement 已写入主稿
- 英文主稿 DOCX 已加入连续行号和页码；补充材料 DOCX 已加入页码
- subscription 路线可投，不强制 APC

## 建议命令

检查中英文主稿/补充材料一致性：

```bash
python3 scripts/check_bmb_cn_consistency.py --lang both
```

重建中英文主稿、补充材料与外围投稿材料 `docx`：

```bash
bash scripts/render_bmb_submission_bundle.sh
```

一键执行提交前自检与重导出：

```bash
bash scripts/preflight_bmb_release.sh
```

编辑阶段只重建 DOCX 并跳过 PDF/上传候选目录：

```bash
BMB_SKIP_PDF=1 bash scripts/preflight_bmb_release.sh
```

单独重建英文上传候选包：

```bash
bash scripts/assemble_bmb_upload_candidate.sh
```

## 仍需作者侧最终确认的项目

- Wenzhuo Wang 和 Yuebin Zhang 的 ORCID
- 最终 `Author contributions` 表述
- 选择并加入软件许可证；当前尚未替作者决定法律许可条款
- 冻结后生成并逐页检查英文 `Online Resource 1` PDF，再组装上传目录和最终 ZIP
- 是否坚持以 `Original research article` 身份提交，还是改走 `Methods`

## 当前判断

按期刊 scope 和现有证据链，`paper1` 以 `Original research article` 口径更稳。它的主线不是单独推出一个算法组件，而是把 study-aware label diagnosis、ortholog validation 和 cross-species prior-enhanced prediction 组织成一个可审计的生物学 benchmark。
