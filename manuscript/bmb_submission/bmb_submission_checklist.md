# BMB Submission Checklist

更新时间：2026-07-16

当前说明：中文主文与中文补充材料保留为内部核对稿；英文主稿、英文补充、英文 title page、英文 cover letter 和英文 declarations 已同步并可重现导出。最终 zip 只在全部预检通过后生成一次。

## Journal fit

- [x] 目标期刊改为 *Bulletin of Mathematical Biology*
- [x] scope 已核对，当前稿件可按 `Original research article` 口径组织
- [x] 投稿模式确认可走 `subscription`，不强制 APC

## Main manuscript

- [x] 主稿已切到 BMB 标题与投稿口径
- [x] abstract 控制在 `150-250` 词范围内（当前自动审计为 `233` 词）
- [x] keywords 收敛到 `6` 个
- [x] 英文主稿和英文补充材料文内引文已改为作者-年份制
- [x] `Statements and Declarations` 已补入主稿结构
- [x] `Data Availability Statement` 已写入
- [x] MOESM2/MOESM3、Ensembl 与 Saluki Zenodo datapack（DOI `10.5281/zenodo.6326409`）来源已明确
- [x] Methods 已加入生成式 AI 使用披露，并明确分析数据、定量结果和最终图件均非 AI 生成
- [x] 英文主稿 DOCX 已加入连续行号
- [x] 英文主稿与英文补充 DOCX 已加入页码
- [x] 中英文正文公式（1）-（9）及补充式（S1）-（S2）均为 Word 原生可编辑 OMML，编号连续且右对齐
- [ ] 作者贡献声明最终确认

## Title page / author metadata

- [x] 作者顺序、单位和通讯作者身份已写入
- [x] Ying Shao 通讯作者邮箱已核对：`yshao@dlmu.edu.cn`
- [x] no-funding statement 已写入，旧基金编号已删除
- [x] competing interests 已写入
- [ ] Wenzhuo Wang ORCID 最终核对
- [x] Ying Shao ORCID 已核对：`0000-0002-4056-5757`

## References

- [x] BibTeX 库已复制到 `references/bmb_references.bib`
- [x] RIS 库已复制到 `references/bmb_references_for_endnote.ris`
- [x] 参考文献改为 BMB 所需作者-年份风格
- [x] 英文主文参考文献已通过自动检查：无未引用条目，且按第一作者字母顺序排列
- [x] 主文 46 条、补充 11 条参考文献的唯一并集为 47 条；45 条 DOI 均通过 Crossref 解析，其中 44 条自动匹配题名/年份/卷，XGBoost 按 ACM 正式记录人工核对
- [x] 无 DOI 的 LightGBM 与 scikit-learn 条目已分别用 NeurIPS 和 JMLR 官方记录核对
- [x] BibTeX/RIS 已清除未进入正文或补充材料的游离条目，当前均为 47 条

## Figures and tables

- [x] 主图与补充图均保留高分辨率源文件
- [x] `PNG/TIFF` 双版本可用
- [x] figure captions 在正文/补充正文中
- [x] 表格保留可编辑源格式
- [x] 12 张当前 PNG 图件已完成视觉检查，未见越界、重叠、遗漏 panel label 或不可读问题
- [x] 当前组合图均提供 600 dpi lossless TIFF，并保留 SVG 矢量源
- [ ] 稿件冻结后重建最终 PDF 图件并复核，避免使用编辑阶段的旧 PDF

## Supplementary information

- [x] supplement 已单独成稿
- [x] supplement 文件头已加入 journal / title / author / corresponding author 信息
- [ ] 稿件冻结后生成英文 supplementary PDF，并复核页码、metadata、图表分页和文字可读性
- [x] 主文已按 Springer 口径把整份补充称为 `Online Resource 1`
- [ ] 若系统要求分拆上传，补充文件命名改成 `ESM_*`

## Submission package

- [x] 中文内部 title page 草案已整理
- [x] 中文内部 cover letter 草案已整理
- [x] 中文内部投稿要点摘要已整理（非系统必传）
- [x] 中文内部 author contributions 已按两位作者重写（待作者确认实际分工）
- [x] 中文内部 declarations 汇总已整理
- [x] 英文 title page 已生成
- [x] 英文 cover letter 已生成
- [x] 英文 `Statements and Declarations` 已生成
- [x] 投稿包自动审计脚本已加入 preflight，当前 DOCX-only 审计为 `0 failure / 0 warning / 165 checks`
- [x] 关键 OOF 完整性审计通过：9 个 global、3 个 cross-target、18 个 target-shrinkage 结果文件均覆盖规范基因宇宙且无缺失预测；residual 11,107 genes 完整
- [x] requirements note 已整理
- [x] README 已整理
- [ ] 英文 Online Resource 1 PDF 冻结后重建 `upload_candidate/`
- [ ] 最终英文上传候选 ZIP 仅在全部作者信息、PDF 和 GitHub 版本冻结后生成一次，并执行 `unzip -t` 与 SHA-256 校验
- [x] GitHub 仓库 `wwzdl/mrna-pc1-label` 已公开，Data/Code Availability 已改为当前时态；本轮严格科学与投稿编辑审计版本标记为 `mRNA-PC1-label-v1.4.3`
- [ ] 选择并加入适合软件公开发布的许可证
- [ ] 上传系统前由两位作者确认作者顺序、Author Contributions、no-funding statement 和英文单位名
