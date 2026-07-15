from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class BiasAwareConfig:
    max_iter: int = 200
    tol: float = 1e-8


@dataclass
class BiasAwareResult:
    gene_consensus: pd.Series
    study_effects: pd.Series
    method_effects: pd.Series
    gene_uncertainty: pd.Series
    fitted_matrix: pd.DataFrame
    residual_matrix: pd.DataFrame
    iterations: int


def _center_series(values: pd.Series) -> pd.Series:
    return values - values.mean()


def fit_bias_aware_consensus(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    config: BiasAwareConfig | None = None,
) -> BiasAwareResult:
    if config is None:
        config = BiasAwareConfig()

    if list(matrix.columns) != metadata["sample_id"].tolist():
        raise ValueError("Matrix columns must match metadata sample order.")

    samples = list(matrix.columns)
    studies = sorted(metadata["study"].unique().tolist())
    methods = sorted(metadata["method"].unique().tolist())

    values = matrix.to_numpy(dtype=float)
    sample_to_study = metadata["study"].astype(str).tolist()
    sample_to_method = metadata["method"].astype(str).tolist()
    study_index = np.array([studies.index(study) for study in sample_to_study], dtype=int)
    method_index = np.array([methods.index(method) for method in sample_to_method], dtype=int)

    alpha_values = np.zeros(len(studies), dtype=float)
    beta_values = np.zeros(len(methods), dtype=float)
    theta_values = np.nanmean(values, axis=1)
    previous = theta_values.copy()

    for iteration in range(1, config.max_iter + 1):
        sample_offsets = alpha_values[study_index] + beta_values[method_index]
        theta_values = np.nanmean(values - sample_offsets[None, :], axis=1)

        residual_without_study = values - theta_values[:, None] - beta_values[method_index][None, :]
        for study_pos in range(len(studies)):
            study_values = residual_without_study[:, study_index == study_pos]
            alpha_values[study_pos] = float(np.nanmean(study_values)) if study_values.size else 0.0
        alpha_values = alpha_values - np.mean(alpha_values)

        residual_without_method = values - theta_values[:, None] - alpha_values[study_index][None, :]
        for method_pos in range(len(methods)):
            method_values = residual_without_method[:, method_index == method_pos]
            beta_values[method_pos] = float(np.nanmean(method_values)) if method_values.size else 0.0
        beta_values = beta_values - np.mean(beta_values)

        max_change = float(np.max(np.abs(theta_values - previous)))
        previous = theta_values.copy()
        if max_change < config.tol:
            break

    theta = pd.Series(theta_values, index=matrix.index, name="bias_aware_consensus")
    alpha = pd.Series(alpha_values, index=studies, name="study_effect")
    beta = pd.Series(beta_values, index=methods, name="method_effect")

    fitted_values = theta_values[:, None] + alpha_values[study_index][None, :] + beta_values[method_index][None, :]
    fitted = pd.DataFrame(fitted_values, index=matrix.index, columns=samples)
    residual = matrix - fitted

    gene_uncertainty = residual.std(axis=1, skipna=True, ddof=0).fillna(0.0).rename("gene_uncertainty")

    return BiasAwareResult(
        gene_consensus=theta,
        study_effects=alpha,
        method_effects=beta,
        gene_uncertainty=gene_uncertainty,
        fitted_matrix=fitted,
        residual_matrix=residual,
        iterations=iteration,
    )
