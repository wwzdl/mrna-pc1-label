# Metadata Format

The paper 1 workflow expects a tab-separated sample metadata file with at least the
following columns:

- `sample_id`: sample name, must match matrix columns exactly
- `study`: study or lab identifier
- `method`: half-life measurement method
- `cell_type`: biological source or cell line
- `species`: species label such as `human` or `mouse`

Optional columns can be added freely and will be preserved.

For the released Saluki supplementary data, this project auto-generates:

- `metadata/real_human_samples.tsv`
- `metadata/real_mouse_samples.tsv`

These files are parsed from the supplement sample names and include an additional
`replicate` column.

The main expression matrix should be a TSV/CSV file with:

- first column: gene identifier
- remaining columns: sample values

One row per gene and one column per sample.
