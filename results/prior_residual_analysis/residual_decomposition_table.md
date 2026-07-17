# Prior residual decomposition result table

| model | genes | prior_features | human_features | pearson | spearman | r2 | delta_r2_vs_prior | remaining_variance_explained |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| prior_only_linear | 11107 | 2 | 0 | 0.792 | 0.782 | 0.627 | 0.000 | 0.00% |
| compact_all_only | 11107 | 0 | 1802 | 0.739 | 0.731 | 0.544 | -0.084 |  |
| prior_plus_compact_residual | 11107 | 2 | 1802 | 0.826 | 0.817 | 0.675 | 0.048 | 12.85% |
