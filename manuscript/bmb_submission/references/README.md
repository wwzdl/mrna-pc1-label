# Reference Library Audit

The active reference-manager library is deliberately limited to references that appear in the English main manuscript or Supplementary Information.

- Main manuscript reference list: 44 entries
- Supplementary reference list: 10 entries
- Unique union across both files: 45 entries
- Entries with DOI: 43
- Entries without DOI: 2 (`ke2017lightgbm`, `pedregosa2011sklearn`)

On 2026-07-14, all 43 DOI records resolved through Crossref and matched the recorded title, year, and volume. The two records without DOI were checked against the official NeurIPS proceedings page and the Journal of Machine Learning Research article page. `bmb_references.bib` is the source of truth; regenerate the EndNote-compatible RIS file with:

```bash
python scripts/export_bib_to_ris.py \
  manuscript/bmb_submission/references/bmb_references.bib \
  manuscript/bmb_submission/references/bmb_references_for_endnote.ris
```
