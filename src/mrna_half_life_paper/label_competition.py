from __future__ import annotations

import argparse
import json
import os
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.benchmark import label_correlation
from mrna_half_life_paper.config import MANUSCRIPT_DIR, METADATA_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.global_mouse_prior import _build_global_settings, _load_compact_all_features
from mrna_half_life_paper.io import align_matrix_and_metadata, load_matrix, load_metadata
from mrna_half_life_paper.ortholog_analysis import ONE2ONE_ORTHOLOG_PATH, extract_mouse_human_one2one
from mrna_half_life_paper.pc1_consensus import compute_pc1_consensus, compute_study_weighted_pc1
from mrna_half_life_paper.preprocess import PreprocessConfig, preprocess_matrix
from mrna_half_life_paper.saluki_naming import saluki_label_column, saluki_label_filename

try:
    from xgboost import XGBRegressor
    from xgboost.core import XGBoostError
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None
    XGBoostError = RuntimeError


@dataclass(frozen=True)
class RecipeSpec:
    species: str
    matrix_source: str
    metadata_variant: str
    min_observed_per_gene: int | None
    preprocess: bool
    weighted: bool


@dataclass(frozen=True)
class CandidateSpec:
    name: str
    human_key: str
    mouse_key: str
    stability_key: str
    description: str


HUMAN_RECIPE_SPECS = {
    "saluki_human_pc1": None,
    "saluki_proxy_pc1": RecipeSpec(
        species="human",
        matrix_source="moesm3",
        metadata_variant="no_gejman",
        min_observed_per_gene=None,
        preprocess=False,
        weighted=False,
    ),
    "saluki_like_full_pc1": RecipeSpec(
        species="human",
        matrix_source="moesm2",
        metadata_variant="full",
        min_observed_per_gene=10,
        preprocess=True,
        weighted=False,
    ),
    "saluki_like_no_gejman_pc1": RecipeSpec(
        species="human",
        matrix_source="moesm2",
        metadata_variant="no_gejman",
        min_observed_per_gene=10,
        preprocess=True,
        weighted=False,
    ),
    "study_weighted_full_pc1": RecipeSpec(
        species="human",
        matrix_source="moesm2",
        metadata_variant="full",
        min_observed_per_gene=10,
        preprocess=True,
        weighted=True,
    ),
    "study_weighted_no_gejman_pc1": RecipeSpec(
        species="human",
        matrix_source="moesm2",
        metadata_variant="no_gejman",
        min_observed_per_gene=10,
        preprocess=True,
        weighted=True,
    ),
}

MOUSE_RECIPE_SPECS = {
    "saluki_human_pc1": None,
    "saluki_like_pc1": RecipeSpec(
        species="mouse",
        matrix_source="moesm2",
        metadata_variant="full",
        min_observed_per_gene=5,
        preprocess=True,
        weighted=False,
    ),
    "study_weighted_pc1": RecipeSpec(
        species="mouse",
        matrix_source="moesm2",
        metadata_variant="full",
        min_observed_per_gene=5,
        preprocess=True,
        weighted=True,
    ),
}

CANDIDATES = [
    CandidateSpec(
        name="saluki_human_pc1",
        human_key="saluki_human_pc1",
        mouse_key="saluki_human_pc1",
        stability_key="saluki_proxy_pc1",
        description="Released Saluki PC1 reference target.",
    ),
    CandidateSpec(
        name="saluki_like_full_pc1",
        human_key="saluki_like_full_pc1",
        mouse_key="saluki_like_pc1",
        stability_key="saluki_like_full_pc1",
        description="MOESM2 reconstruction with Saluki-style coverage thresholds.",
    ),
    CandidateSpec(
        name="saluki_like_no_gejman_pc1",
        human_key="saluki_like_no_gejman_pc1",
        mouse_key="saluki_like_pc1",
        stability_key="saluki_like_no_gejman_pc1",
        description="Saluki-style thresholding plus Gejman removal in human.",
    ),
    CandidateSpec(
        name="study_weighted_full_pc1",
        human_key="study_weighted_full_pc1",
        mouse_key="study_weighted_pc1",
        stability_key="study_weighted_full_pc1",
        description="Study-balanced PCA without manual pruning.",
    ),
    CandidateSpec(
        name="study_weighted_no_gejman_pc1",
        human_key="study_weighted_no_gejman_pc1",
        mouse_key="study_weighted_pc1",
        stability_key="study_weighted_no_gejman_pc1",
        description="Study-balanced PCA plus Gejman removal in human.",
    ),
]


DISPLAY_CANDIDATE_NAMES = {
    "saluki_human_pc1": "Saluki human PC1",
}


def _display_candidate_name(name: str) -> str:
    return DISPLAY_CANDIDATE_NAMES.get(name, name)


def _matrix_path(species: str, matrix_source: str) -> Path:
    file_name = {
        "moesm2": "moesm2_raw_sparse_matrix.tsv",
        "moesm3": "moesm3_saluki_processed_matrix.tsv",
    }[matrix_source]
    return PROCESSED_DIR / "real_data" / species / file_name


def _metadata_path(species: str, metadata_variant: str) -> Path:
    if species == "human" and metadata_variant == "no_gejman":
        return METADATA_DIR / "real_human_samples_no_Gejman.tsv"
    return METADATA_DIR / f"real_{species}_samples.tsv"


def _saluki_path(species: str) -> Path:
    return PROCESSED_DIR / "real_data" / species / saluki_label_filename(species)


def _load_recipe_inputs(spec: RecipeSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = load_matrix(_matrix_path(spec.species, spec.matrix_source))
    metadata = load_metadata(_metadata_path(spec.species, spec.metadata_variant))
    return align_matrix_and_metadata(matrix, metadata)


def _compute_recipe_label(
    spec: RecipeSpec,
    *,
    matrix: pd.DataFrame | None = None,
    metadata: pd.DataFrame | None = None,
) -> pd.Series:
    if matrix is None or metadata is None:
        matrix, metadata = _load_recipe_inputs(spec)

    working = matrix.copy()
    if spec.preprocess:
        if spec.min_observed_per_gene is None:
            raise ValueError("Preprocessed recipes require min_observed_per_gene.")
        working = preprocess_matrix(
            working,
            PreprocessConfig(
                apply_log=False,
                min_observed_per_gene=spec.min_observed_per_gene,
            ),
        )
    if spec.weighted:
        return compute_study_weighted_pc1(working, metadata)
    return compute_pc1_consensus(working)


def _load_saluki_label(species: str) -> pd.Series:
    frame = pd.read_csv(_saluki_path(species), sep="\t")
    column = saluki_label_column(species)
    return frame.set_index("gene_id")[column].astype(float).rename(column)


def _persist_series(series: pd.Series, path: Path, column_name: str = "pc1") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    series.rename(column_name).to_frame().reset_index().rename(columns={"index": "gene_id"}).to_csv(
        path,
        sep="\t",
        index=False,
    )
    return path


def _summary_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    return {
        "pearson": float(np.corrcoef(y_true, y_pred)[0, 1]),
        "spearman": float(spearmanr(y_true, y_pred).statistic),
        "r2": float(r2_score(y_true, y_pred)),
    }


def _xgb_oof_predictions(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    random_state: int = 42,
    cv_splits: int = 5,
    device: str | None = None,
) -> tuple[np.ndarray, list[int], str]:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for label competition benchmarks.")

    chosen_device = device if device is not None else os.environ.get("MRNA_HALF_LIFE_XGB_DEVICE", "cuda")
    x = df[feature_cols].to_numpy(dtype=np.float32)
    y = df[target_col].to_numpy(dtype=np.float32)
    predictions = np.full(len(y), np.nan, dtype=np.float32)
    best_iterations: list[int] = []
    kfold = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)

    def fit_one_fold(fit_idx: np.ndarray, val_idx: np.ndarray, test_idx: np.ndarray, model_device: str) -> int:
        model = XGBRegressor(
            n_estimators=3000,
            early_stopping_rounds=100,
            objective="reg:squarederror",
            tree_method="hist",
            device=model_device,
            max_bin=256,
            max_depth=6,
            learning_rate=0.02,
            min_child_weight=4,
            subsample=0.9,
            colsample_bytree=0.75,
            reg_alpha=0.0,
            reg_lambda=1.0,
            random_state=random_state,
        )
        model.fit(
            x[fit_idx],
            y[fit_idx],
            eval_set=[(x[val_idx], y[val_idx])],
            verbose=False,
        )
        best_iteration = int(getattr(model, "best_iteration", 2999))
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                category=UserWarning,
            )
            predictions[test_idx] = model.predict(x[test_idx], iteration_range=(0, best_iteration + 1))
        return best_iteration

    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(np.arange(len(y))), start=1):
        rng = np.random.RandomState(random_state + fold_idx)
        perm = rng.permutation(train_idx.size)
        val_n = max(512, int(0.1 * train_idx.size))
        val_rel = perm[:val_n]
        fit_rel = perm[val_n:]
        fit_idx = train_idx[fit_rel]
        val_idx = train_idx[val_rel]
        try:
            best_iteration = fit_one_fold(fit_idx, val_idx, test_idx, chosen_device)
        except XGBoostError:
            if chosen_device == "cpu":
                raise
            chosen_device = "cpu"
            best_iteration = fit_one_fold(fit_idx, val_idx, test_idx, chosen_device)
        best_iterations.append(best_iteration)

    return predictions, best_iterations, chosen_device


def _load_prior_feature_table(feature_path: Path) -> tuple[pd.DataFrame, list[str], dict[str, list[str]]]:
    features = _load_compact_all_features(feature_path)
    common_path = RESULTS_DIR / "ortholog_analysis" / "common_ortholog_pairs.tsv"
    saluki_path = RESULTS_DIR / "ortholog_analysis" / "saluki_ortholog_pairs.tsv"
    if not common_path.exists() or not saluki_path.exists():
        raise FileNotFoundError(
            "Missing ortholog prior tables. Run mrna_half_life_paper.ortholog_analysis first."
        )

    mouse_common = pd.read_csv(common_path, sep="\t")[["human_gene_id", "mouse_pc1", "is_high_confidence"]].rename(
        columns={"human_gene_id": "gene_id"}
    )
    saluki_mouse = pd.read_csv(saluki_path, sep="\t")[
        ["human_gene_id", "saluki_mouse_prior", "is_high_confidence"]
    ].rename(columns={"human_gene_id": "gene_id", "is_high_confidence": "saluki_mouse_high_confidence"})

    df = features.merge(mouse_common, on="gene_id", how="left").merge(saluki_mouse, on="gene_id", how="left")
    df["has_mouse_pc1"] = df["mouse_pc1"].notna().astype(np.float32)
    df["mouse_pc1_highconf"] = df["is_high_confidence"].fillna(0).astype(np.float32)
    df["mouse_pc1"] = df["mouse_pc1"].fillna(0.0)
    df["has_saluki_mouse_prior"] = df["saluki_mouse_prior"].notna().astype(np.float32)
    df["saluki_mouse_high_confidence"] = df["saluki_mouse_high_confidence"].fillna(0).astype(np.float32)
    df["saluki_mouse_prior"] = df["saluki_mouse_prior"].fillna(0.0)
    df = df.drop(columns=["is_high_confidence"]).sort_values("gene_id").reset_index(drop=True)

    base_cols = [col for col in features.columns if col != "gene_id"]
    return df, base_cols, _build_global_settings(base_cols)


def _build_shared_ortholog_table(
    candidate_labels: dict[str, tuple[pd.Series, pd.Series]],
) -> pd.DataFrame:
    extract_mouse_human_one2one()
    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")
    for human_label, mouse_label in candidate_labels.values():
        mapping = mapping[
            mapping["human_gene_id"].isin(human_label.index)
            & mapping["mouse_gene_id"].isin(mouse_label.index)
        ].copy()
    return mapping.sort_values(["human_gene_id", "mouse_gene_id"]).reset_index(drop=True)


def _benchmark_stability(
    candidate: CandidateSpec,
    reference_series: pd.Series,
    recipe: RecipeSpec,
) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    matrix, metadata = _load_recipe_inputs(recipe)
    rows: list[dict[str, float | int | str]] = []
    for study in sorted(metadata["study"].unique().tolist()):
        keep_meta = metadata.loc[metadata["study"] != study].copy()
        keep_samples = keep_meta["sample_id"].tolist()
        subset_matrix = matrix.loc[:, keep_samples].copy()
        subset_label = _compute_recipe_label(recipe, matrix=subset_matrix, metadata=keep_meta)
        corr = label_correlation(reference_series, subset_label)
        rows.append(
            {
                "candidate": candidate.name,
                "left_out": study,
                "n_samples_removed": int(metadata.loc[metadata["study"] == study].shape[0]),
                "n_genes_overlap": int(corr["n"]),
                "pearson": corr["pearson"],
                "spearman": corr["spearman"],
            }
        )

    df = pd.DataFrame(rows).sort_values("pearson", ascending=True).reset_index(drop=True)
    summary = {
        "candidate": candidate.name,
        "n_studies": int(df.shape[0]),
        "stability_mean_pearson": float(df["pearson"].mean()),
        "stability_min_pearson": float(df["pearson"].min()),
        "stability_mean_spearman": float(df["spearman"].mean()),
        "stability_worst_left_out": str(df.iloc[0]["left_out"]),
    }
    return df, summary


def _benchmark_orthologs(
    candidate_labels: dict[str, tuple[pd.Series, pd.Series]],
    shared_mapping: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, float | int | str]] = []
    pair_tables: list[pd.DataFrame] = []
    for candidate_name, (human_label, mouse_label) in candidate_labels.items():
        pairs = shared_mapping.copy()
        pairs["human_pc1"] = pairs["human_gene_id"].map(human_label)
        pairs["mouse_pc1"] = pairs["mouse_gene_id"].map(mouse_label)
        corr = label_correlation(
            pairs.set_index("human_gene_id")["human_pc1"],
            pairs.set_index("human_gene_id")["mouse_pc1"],
        )
        rows.append(
            {
                "candidate": candidate_name,
                "n_pairs": int(corr["n"]),
                "pearson": corr["pearson"],
                "spearman": corr["spearman"],
            }
        )
        pair_tables.append(
            pairs.assign(candidate=candidate_name).loc[
                :, ["candidate", "human_gene_id", "mouse_gene_id", "is_high_confidence", "human_pc1", "mouse_pc1"]
            ]
        )
    return pd.DataFrame(rows).sort_values("pearson", ascending=False), pd.concat(pair_tables, ignore_index=True)


def _write_markdown_note(overall: pd.DataFrame, outpath: Path) -> Path:
    leader = overall.sort_values(
        ["ortholog_pearson", "stability_mean_pearson", "pure_human_pearson_vs_target"],
        ascending=False,
    ).iloc[0]
    best_pure = overall.sort_values(["pure_human_pearson_vs_target"], ascending=False).iloc[0]
    best_prior = overall.sort_values(["prior_enhanced_pearson_vs_target"], ascending=False).iloc[0]
    best_ortholog = overall.sort_values(["ortholog_pearson"], ascending=False).iloc[0]
    best_stable = overall.sort_values(["stability_mean_pearson"], ascending=False).iloc[0]

    lines = [
        "# Candidate Label Competition",
        "",
        "## Goal",
        "",
        "在统一评估口径下比较 Saluki human PC1、Saluki 对齐重建标签，以及 study-balanced 标签，",
        "看我们自己的标签能否在外部一致性、study 稳定性和下游可学习性上拿到更扎实的优势。",
        "",
        "## Headline",
        "",
        f"- 综合最靠前的候选是 `{_display_candidate_name(leader['candidate'])}`。",
        f"- ortholog 一致性最好的是 `{_display_candidate_name(best_ortholog['candidate'])}`，Pearson={best_ortholog['ortholog_pearson']:.3f}。",
        f"- leave-one-study-out 平均稳定性最好的是 `{_display_candidate_name(best_stable['candidate'])}`，Mean Pearson={best_stable['stability_mean_pearson']:.3f}。",
        f"- pure human OOF 最好学的是 `{_display_candidate_name(best_pure['candidate'])}`，Pearson(target)={best_pure['pure_human_pearson_vs_target']:.3f}。",
        f"- prior-enhanced OOF 最好学的是 `{_display_candidate_name(best_prior['candidate'])}`，Pearson(target)={best_prior['prior_enhanced_pearson_vs_target']:.3f}。",
        "",
        "## How To Read",
        "",
        "- 如果某个标签只是在“与 Saluki human PC1 的一致性”上更高，它更像是在逼近 Saluki 定义。",
        "- 如果某个标签同时提升了 ortholog 一致性、study 稳定性和 pure/prior OOF，可更有底气说它是“更好的标签”。",
        "",
    ]
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text("\n".join(lines), encoding="utf-8")
    return outpath


def run_label_competition(
    feature_path: Path | None = None,
    results_root: Path | None = None,
    device: str | None = None,
) -> dict[str, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    results_root = results_root if results_root is not None else RESULTS_DIR / "label_competition"
    results_root.mkdir(parents=True, exist_ok=True)
    target_dir = results_root / "targets"
    target_dir.mkdir(parents=True, exist_ok=True)
    oof_dir = results_root / "oof_predictions"
    oof_dir.mkdir(parents=True, exist_ok=True)

    human_labels: dict[str, pd.Series] = {"saluki_human_pc1": _load_saluki_label("human")}
    for key, spec in HUMAN_RECIPE_SPECS.items():
        if spec is None:
            continue
        human_labels[key] = _compute_recipe_label(spec)
    mouse_labels: dict[str, pd.Series] = {"saluki_human_pc1": _load_saluki_label("mouse")}
    for key, spec in MOUSE_RECIPE_SPECS.items():
        if spec is None:
            continue
        mouse_labels[key] = _compute_recipe_label(spec)

    for key, series in human_labels.items():
        column = "saluki_human_pc1" if key == "saluki_human_pc1" else "pc1"
        _persist_series(series, target_dir / "human" / f"{key}.tsv", column_name=column)
    for key, series in mouse_labels.items():
        column = "saluki_human_pc1" if key == "saluki_human_pc1" else "pc1"
        _persist_series(series, target_dir / "mouse" / f"{key}.tsv", column_name=column)

    candidate_pairs = {
        candidate.name: (human_labels[candidate.human_key], mouse_labels[candidate.mouse_key]) for candidate in CANDIDATES
    }
    candidate_info = pd.DataFrame(
        [
            {
                "candidate": candidate.name,
                "description": candidate.description,
                "human_source": candidate.human_key,
                "mouse_source": candidate.mouse_key,
                "n_human_genes": int(human_labels[candidate.human_key].shape[0]),
                "n_mouse_genes": int(mouse_labels[candidate.mouse_key].shape[0]),
            }
            for candidate in CANDIDATES
        ]
    )
    candidate_info.to_csv(results_root / "candidate_info.tsv", sep="\t", index=False)

    saluki_reference = human_labels["saluki_human_pc1"]
    saluki_alignment_rows = []
    for candidate in CANDIDATES:
        corr = label_correlation(human_labels[candidate.human_key], saluki_reference)
        saluki_alignment_rows.append(
            {
                "candidate": candidate.name,
                "n_overlap": int(corr["n"]),
                "human_vs_saluki_pearson": corr["pearson"],
                "human_vs_saluki_spearman": corr["spearman"],
            }
        )
    saluki_alignment = pd.DataFrame(saluki_alignment_rows).sort_values("human_vs_saluki_pearson", ascending=False)
    saluki_alignment.to_csv(results_root / "human_saluki_alignment.tsv", sep="\t", index=False)

    stability_rows = []
    stability_summary_rows = []
    for candidate in CANDIDATES:
        recipe = HUMAN_RECIPE_SPECS[candidate.stability_key]
        if recipe is None:
            raise RuntimeError(f"Missing stability recipe for {candidate.name}")
        reference = human_labels[candidate.human_key]
        df, summary = _benchmark_stability(candidate, reference, recipe)
        stability_rows.append(df)
        stability_summary_rows.append(summary)
    stability_detail = pd.concat(stability_rows, ignore_index=True)
    stability_summary = pd.DataFrame(stability_summary_rows).sort_values(
        "stability_mean_pearson", ascending=False
    )
    stability_detail.to_csv(results_root / "stability_detail.tsv", sep="\t", index=False)
    stability_summary.to_csv(results_root / "stability_summary.tsv", sep="\t", index=False)

    shared_mapping = _build_shared_ortholog_table(candidate_pairs)
    ortholog_summary, ortholog_pairs = _benchmark_orthologs(candidate_pairs, shared_mapping)
    ortholog_summary.to_csv(results_root / "ortholog_summary.tsv", sep="\t", index=False)
    ortholog_pairs.to_csv(results_root / "ortholog_pairs.tsv", sep="\t", index=False)

    compact_all_features = _load_compact_all_features(feature_path).sort_values("gene_id").reset_index(drop=True)
    prior_features, base_cols, global_settings = _load_prior_feature_table(feature_path)
    prior_setting_cols = global_settings["compact_all_plus_both_mouse_priors_global"]
    pure_setting_cols = [col for col in compact_all_features.columns if col != "gene_id"]

    shared_genes = set(compact_all_features["gene_id"].astype(str))
    shared_genes &= set(prior_features["gene_id"].astype(str))
    shared_genes &= set(saluki_reference.index.astype(str))
    for candidate in CANDIDATES:
        shared_genes &= set(human_labels[candidate.human_key].index.astype(str))
    shared_gene_list = sorted(shared_genes)
    pd.DataFrame({"gene_id": shared_gene_list}).to_csv(results_root / "shared_gene_universe.tsv", sep="\t", index=False)

    pure_base = compact_all_features.loc[compact_all_features["gene_id"].isin(shared_genes)].copy()
    pure_base = pure_base.sort_values("gene_id").reset_index(drop=True)
    prior_base = prior_features.loc[prior_features["gene_id"].isin(shared_genes)].copy()
    prior_base = prior_base.sort_values("gene_id").reset_index(drop=True)

    pure_rows = []
    prior_rows = []
    devices_used: set[str] = set()
    saluki_ref_df = saluki_reference.rename("saluki_reference_pc1").to_frame().reset_index().rename(
        columns={"index": "gene_id"}
    )
    for candidate in CANDIDATES:
        target_df = (
            human_labels[candidate.human_key]
            .rename(candidate.name)
            .to_frame()
            .reset_index()
            .rename(columns={"index": "gene_id"})
        )
        merged_pure = (
            pure_base.merge(target_df, on="gene_id", how="inner")
            .merge(saluki_ref_df, on="gene_id", how="inner")
            .sort_values("gene_id")
            .reset_index(drop=True)
        )
        pred, best_iter, device_used = _xgb_oof_predictions(
            merged_pure,
            pure_setting_cols,
            candidate.name,
            device=device,
        )
        devices_used.add(device_used)
        target_values = merged_pure[candidate.name].to_numpy(dtype=np.float32)
        saluki_values = merged_pure["saluki_reference_pc1"].to_numpy(dtype=np.float32)
        target_corr = label_correlation(
            merged_pure.set_index("gene_id")[candidate.name],
            merged_pure.set_index("gene_id")["saluki_reference_pc1"],
        )
        pure_rows.append(
            {
                "candidate": candidate.name,
                "n": int(len(merged_pure)),
                "n_features": int(len(pure_setting_cols)),
                "target_vs_saluki_pearson": target_corr["pearson"],
                "target_vs_saluki_spearman": target_corr["spearman"],
                "pearson_vs_target": _summary_metrics(target_values, pred)["pearson"],
                "spearman_vs_target": _summary_metrics(target_values, pred)["spearman"],
                "r2_vs_target": _summary_metrics(target_values, pred)["r2"],
                "pearson_vs_saluki_reference": _summary_metrics(saluki_values, pred)["pearson"],
                "spearman_vs_saluki_reference": _summary_metrics(saluki_values, pred)["spearman"],
                "r2_vs_saluki_reference": _summary_metrics(saluki_values, pred)["r2"],
                "mean_best_iter": float(np.mean(best_iter)),
                "device_used": device_used,
            }
        )
        oof_pure = merged_pure[["gene_id", candidate.name, "saluki_reference_pc1"]].copy()
        oof_pure["prediction"] = pred
        oof_pure.to_csv(oof_dir / f"{candidate.name}__pure_human.tsv", sep="\t", index=False)

        merged_prior = (
            prior_base.merge(target_df, on="gene_id", how="inner")
            .merge(saluki_ref_df, on="gene_id", how="inner")
            .sort_values("gene_id")
            .reset_index(drop=True)
        )
        pred, best_iter, device_used = _xgb_oof_predictions(
            merged_prior,
            prior_setting_cols,
            candidate.name,
            device=device,
        )
        devices_used.add(device_used)
        target_values = merged_prior[candidate.name].to_numpy(dtype=np.float32)
        saluki_values = merged_prior["saluki_reference_pc1"].to_numpy(dtype=np.float32)
        prior_rows.append(
            {
                "candidate": candidate.name,
                "n": int(len(merged_prior)),
                "n_features": int(len(prior_setting_cols)),
                "n_with_mouse_pc1": int(merged_prior["has_mouse_pc1"].sum()),
                "n_with_saluki_mouse_prior": int(merged_prior["has_saluki_mouse_prior"].sum()),
                "pearson_vs_target": _summary_metrics(target_values, pred)["pearson"],
                "spearman_vs_target": _summary_metrics(target_values, pred)["spearman"],
                "r2_vs_target": _summary_metrics(target_values, pred)["r2"],
                "pearson_vs_saluki_reference": _summary_metrics(saluki_values, pred)["pearson"],
                "spearman_vs_saluki_reference": _summary_metrics(saluki_values, pred)["spearman"],
                "r2_vs_saluki_reference": _summary_metrics(saluki_values, pred)["r2"],
                "mean_best_iter": float(np.mean(best_iter)),
                "device_used": device_used,
            }
        )
        oof_prior = merged_prior[["gene_id", candidate.name, "saluki_reference_pc1"]].copy()
        oof_prior["prediction"] = pred
        oof_prior.to_csv(oof_dir / f"{candidate.name}__prior_enhanced.tsv", sep="\t", index=False)

    pure_summary = pd.DataFrame(pure_rows).sort_values("pearson_vs_target", ascending=False)
    prior_summary = pd.DataFrame(prior_rows).sort_values("pearson_vs_target", ascending=False)
    pure_summary.to_csv(results_root / "pure_human_summary.tsv", sep="\t", index=False)
    prior_summary.to_csv(results_root / "prior_enhanced_summary.tsv", sep="\t", index=False)

    overall = (
        candidate_info.merge(saluki_alignment, on="candidate", how="left")
        .merge(stability_summary, on="candidate", how="left")
        .merge(
            ortholog_summary.rename(columns={"pearson": "ortholog_pearson", "spearman": "ortholog_spearman"}),
            on="candidate",
            how="left",
        )
        .merge(
            pure_summary.rename(
                columns={
                    "pearson_vs_target": "pure_human_pearson_vs_target",
                    "spearman_vs_target": "pure_human_spearman_vs_target",
                    "r2_vs_target": "pure_human_r2_vs_target",
                    "pearson_vs_saluki_reference": "pure_human_pearson_vs_saluki",
                    "spearman_vs_saluki_reference": "pure_human_spearman_vs_saluki",
                    "r2_vs_saluki_reference": "pure_human_r2_vs_saluki",
                    "device_used": "pure_human_device",
                }
            )[
                [
                    "candidate",
                    "pure_human_pearson_vs_target",
                    "pure_human_spearman_vs_target",
                    "pure_human_r2_vs_target",
                    "pure_human_pearson_vs_saluki",
                    "pure_human_spearman_vs_saluki",
                    "pure_human_r2_vs_saluki",
                    "pure_human_device",
                ]
            ],
            on="candidate",
            how="left",
        )
        .merge(
            prior_summary.rename(
                columns={
                    "pearson_vs_target": "prior_enhanced_pearson_vs_target",
                    "spearman_vs_target": "prior_enhanced_spearman_vs_target",
                    "r2_vs_target": "prior_enhanced_r2_vs_target",
                    "pearson_vs_saluki_reference": "prior_enhanced_pearson_vs_saluki",
                    "spearman_vs_saluki_reference": "prior_enhanced_spearman_vs_saluki",
                    "r2_vs_saluki_reference": "prior_enhanced_r2_vs_saluki",
                    "device_used": "prior_enhanced_device",
                }
            )[
                [
                    "candidate",
                    "prior_enhanced_pearson_vs_target",
                    "prior_enhanced_spearman_vs_target",
                    "prior_enhanced_r2_vs_target",
                    "prior_enhanced_pearson_vs_saluki",
                    "prior_enhanced_spearman_vs_saluki",
                    "prior_enhanced_r2_vs_saluki",
                    "prior_enhanced_device",
                ]
            ],
            on="candidate",
            how="left",
        )
    )
    overall["device_used"] = ",".join(sorted(devices_used))
    overall_path = results_root / "overall_summary.tsv"
    overall.sort_values(
        ["ortholog_pearson", "stability_mean_pearson", "prior_enhanced_pearson_vs_target"],
        ascending=False,
    ).to_csv(overall_path, sep="\t", index=False)
    with open(results_root / "overall_summary.json", "w", encoding="utf-8") as handle:
        json.dump(overall.to_dict(orient="records"), handle, indent=2)

    note_path = _write_markdown_note(overall, MANUSCRIPT_DIR.parent / "docs" / "label_competition.md")
    return {
        "overall_summary": overall_path,
        "saluki_alignment": results_root / "human_saluki_alignment.tsv",
        "stability_summary": results_root / "stability_summary.tsv",
        "ortholog_summary": results_root / "ortholog_summary.tsv",
        "pure_human_summary": results_root / "pure_human_summary.tsv",
        "prior_enhanced_summary": results_root / "prior_enhanced_summary.tsv",
        "note": note_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare candidate human half-life labels under shared benchmarks.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_DIR / "label_competition",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Optional XGBoost device override, for example 'cuda' or 'cpu'.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_label_competition(
        feature_path=args.feature_path,
        results_root=args.results_root,
        device=args.device,
    )
    for path in outputs.values():
        print(path)


if __name__ == "__main__":
    main()
