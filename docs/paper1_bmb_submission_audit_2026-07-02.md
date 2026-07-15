# paper1 BMB submission audit (2026-07-02)

## 目的

将 `paper1` 从此前的 `TCBB` 投稿工作目录切换到 *Bulletin of Mathematical Biology*（BMB）投稿口径，并依据 Springer 官方投稿要求建立新的 submission package。

## 本轮完成内容

### 1. 读取并整理 BMB 官方要求

已核对官方页面：

- `Submission guidelines`
- `Aims and scope`
- `How to publish with us`
- `Journal home / metrics`

本轮明确采用的关键规则：

- `single-blind peer review`
- abstract `150-250` 词
- `4-6` 个 keywords
- 文内引用改为作者-年份制
- original research 必须含 `Data Availability Statement`
- supplementary information 单独提交
- subscription 路线可投，`no APC charges apply`

### 2. 新建 BMB 投稿目录

新增目录：

- `paper_pca/manuscript/bmb_submission/`

主要文件：

- `bmb_main_manuscript.md`
- `bmb_main_manuscript_en.md`
- `bmb_main_manuscript.docx`
- `bmb_supplementary_material.md`
- `bmb_supplementary_material_en.md`
- `bmb_supplementary_material.docx`
- `bmb_title_page_en.md`
- `bmb_title_page_en.docx`
- `bmb_cover_letter_en.md`
- `bmb_cover_letter_en.docx`
- `bmb_submission_requirements.md`
- `bmb_submission_checklist.md`
- `README.md`
- `bmb_submission_package.zip`

### 3. 主稿口径切换

主稿已完成以下调整：

- 标题统一为  
  `Study-Aware Label Diagnosis and Cross-Species Prior-Enhanced Prediction of Mammalian mRNA Half-Life`
- 文章类型改为 `Original research article`
- `Index Terms` 改为 `Keywords`
- 关键词压缩到 `6` 个
- `TCBB` 残留表述移除
- 文末改为 `Statements and Declarations`
- 保留 current strongest claim boundary，不把 prior-enhanced 结果写成 sequence-only claim

### 4. 补充材料口径切换

补充材料已完成以下调整：

- 文件头加入 journal / title / author / corresponding author 信息
- 移除 `IEEE/TCBB` 相关补充说明
- 保留图表与补充表结构
- 保持与主稿相同的作者-年份引文风格

### 5. 引文与导出工具

新增脚本：

- `paper_pca/scripts/render_bmb_markdown_with_references.py`

作用：

- 将 Pandoc 风格 `[@key]` 引文渲染成 BMB 所需作者-年份格式
- 生成按第一作者排序的参考文献表

同时修补：

- `paper_pca/scripts/render_markdown_to_docx.py`

修补内容：

- 自动剥离 YAML front matter，避免导出的 `docx` 把元数据误当正文

## 当前核对结果

- 主稿 abstract 约 `223` 词，符合 BMB `150-250` 词要求
- 主稿关键词数：`6`
- `bmb_submission/` 中未再检出 `TCBB` / `IEEE` 残留短语
- `render_bmb_markdown_with_references.py` 与 `render_markdown_to_docx.py` 已 `py_compile` 通过
- `bmb_submission_package.zip` 已生成

## 仍待作者侧确认

- Wenzhuo Wang ORCID
- Yuebin Zhang ORCID
- `Author contributions` 最终表述
- 连续行号与页码的最终提交版
- supplementary 最终 PDF 命名与 `Online Resource` 口径

## 当前建议

若继续沿 BMB 路线推进，下一步最值得做的是：

1. 手工审一遍 `bmb_main_manuscript.docx` 的标题层级、图题和参考文献可读性  
2. 补齐 `Author contributions`  
3. 准备最终带行号的提交版 PDF 或 LaTeX/Word 正式版
