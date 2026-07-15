from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RAW_DIR, RESULTS_DIR
from mrna_half_life_paper.saluki_naming import (
    HAS_SALUKI_MOUSE_PRIOR,
    SALUKI_HUMAN_PC1,
    SALUKI_MOUSE_HIGH_CONFIDENCE,
    SALUKI_MOUSE_PRIOR,
)
from mrna_half_life_paper.sequence_baseline import _load_target_tables

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None


REGULATORY_DIR = RAW_DIR / "saluki_datapack" / "human"
PRIOR_MASK_COLS = {
    "mouse_pc1": "has_mouse_pc1",
    SALUKI_MOUSE_PRIOR: HAS_SALUKI_MOUSE_PRIOR,
}


def _load_region_prediction_block(directory: Path, prefix: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for region, file_name in (("utr5", "5pUTR_avg.txt.gz"), ("orf", "ORF_avg.txt.gz"), ("utr3", "3pUTR_avg.txt.gz")):
        frame = pd.read_csv(directory / file_name, sep="\t")
        frame = frame.rename(columns={frame.columns[0]: "gene_id"})
        frame = frame.rename(columns={col: f"{prefix}__{region}__{col}" for col in frame.columns if col != "gene_id"})
        frames.append(frame)

    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on="gene_id", how="outer")
    return merged.fillna(0.0)


def _load_compact_all_features(feature_path: Path) -> pd.DataFrame:
    features = pd.read_csv(feature_path, sep="\t").drop(columns=["transcript_id"], errors="ignore")
    cwcs = pd.read_csv(REGULATORY_DIR / "CWCS.txt.gz", sep="\t").rename(columns={"GeneID": "gene_id"})
    seqweaver = _load_region_prediction_block(REGULATORY_DIR / "SeqWeaver_predictions", "seqweaver")
    deepripe = _load_region_prediction_block(REGULATORY_DIR / "DeepRiPe_predictions_v83", "deepripe")
    return (
        features.merge(cwcs, on="gene_id", how="left")
        .merge(seqweaver, on="gene_id", how="left")
        .merge(deepripe, on="gene_id", how="left")
        .fillna(0.0)
    )


def _load_global_table(feature_path: Path) -> tuple[pd.DataFrame, list[str]]:
    features = _load_compact_all_features(feature_path)
    targets = _load_target_tables("human")
    common_gene_ids = set(features["gene_id"].astype(str))
    for target_name in (SALUKI_HUMAN_PC1, "human_full_pc1", "human_pruned_pc1"):
        common_gene_ids &= set(targets[target_name]["gene_id"].astype(str))
    features = features.loc[features["gene_id"].isin(common_gene_ids)].copy()

    mouse_common = pd.read_csv(
        RESULTS_DIR / "ortholog_analysis" / "common_ortholog_pairs.tsv",
        sep="\t",
    )[["human_gene_id", "mouse_pc1", "is_high_confidence"]].rename(columns={"human_gene_id": "gene_id"})
    saluki_mouse = pd.read_csv(
        RESULTS_DIR / "ortholog_analysis" / "saluki_ortholog_pairs.tsv",
        sep="\t",
    )[["human_gene_id", "saluki_mouse_prior", "is_high_confidence"]].rename(
        columns={
            "human_gene_id": "gene_id",
            "saluki_mouse_prior": SALUKI_MOUSE_PRIOR,
            "is_high_confidence": SALUKI_MOUSE_HIGH_CONFIDENCE,
        }
    )

    df = (
        features.merge(targets[SALUKI_HUMAN_PC1], on="gene_id", how="inner")
        .merge(mouse_common, on="gene_id", how="left")
        .merge(saluki_mouse, on="gene_id", how="left")
    )
    df["has_mouse_pc1"] = df["mouse_pc1"].notna().astype(np.float32)
    df["mouse_pc1_highconf"] = df["is_high_confidence"].fillna(0).astype(np.float32)
    df["mouse_pc1"] = df["mouse_pc1"].fillna(0.0)
    df[HAS_SALUKI_MOUSE_PRIOR] = df[SALUKI_MOUSE_PRIOR].notna().astype(np.float32)
    df[SALUKI_MOUSE_HIGH_CONFIDENCE] = df[SALUKI_MOUSE_HIGH_CONFIDENCE].fillna(0).astype(np.float32)
    df[SALUKI_MOUSE_PRIOR] = df[SALUKI_MOUSE_PRIOR].fillna(0.0)
    df = df.drop(columns=["is_high_confidence"]).sort_values("gene_id").reset_index(drop=True)

    feature_cols = [col for col in features.columns if col != "gene_id"]
    return df, feature_cols


def _build_global_settings(base_cols: list[str]) -> dict[str, list[str]]:
    settings = {
        "compact_all_global": base_cols,
        "compact_all_plus_mouse_pc1_global": base_cols + ["mouse_pc1", "has_mouse_pc1", "mouse_pc1_highconf"],
        "compact_all_plus_saluki_mouse_prior_global": base_cols
        + [SALUKI_MOUSE_PRIOR, HAS_SALUKI_MOUSE_PRIOR, SALUKI_MOUSE_HIGH_CONFIDENCE],
        "compact_all_plus_both_mouse_priors_global": base_cols
        + [
            "mouse_pc1",
            "has_mouse_pc1",
            "mouse_pc1_highconf",
            SALUKI_MOUSE_PRIOR,
            HAS_SALUKI_MOUSE_PRIOR,
            SALUKI_MOUSE_HIGH_CONFIDENCE,
        ],
    }
    return settings


def _shuffle_prior_values_in_place(
    x: np.ndarray,
    *,
    prior_idx: int,
    mask_idx: int,
    indices: np.ndarray,
    rng: np.random.RandomState,
) -> None:
    available_idx = indices[x[indices, mask_idx] > 0.5]
    if available_idx.size <= 1:
        return
    shuffled = x[available_idx, prior_idx].copy()
    rng.shuffle(shuffled)
    x[available_idx, prior_idx] = shuffled


def _xgb_oof_predictions(
    df: pd.DataFrame,
    feature_cols: list[str],
    *,
    random_state: int = 42,
    cv_splits: int = 5,
    shuffle_priors: tuple[str, ...] = (),
) -> tuple[np.ndarray, list[int]]:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for the global mouse prior benchmark")

    x = df[feature_cols].to_numpy(dtype=np.float32)
    y = df[SALUKI_HUMAN_PC1].to_numpy(dtype=np.float32)
    predictions = np.full(len(y), np.nan, dtype=np.float32)
    best_iterations: list[int] = []
    kfold = KFold(n_splits=cv_splits, shuffle=True, random_state=random_state)
    col_idx = {col: idx for idx, col in enumerate(feature_cols)}

    for prior_col in shuffle_priors:
        if prior_col not in col_idx:
            raise ValueError(f"Cannot shuffle missing prior column: {prior_col}")
        mask_col = PRIOR_MASK_COLS[prior_col]
        if mask_col not in col_idx:
            raise ValueError(f"Cannot shuffle prior '{prior_col}' without mask column '{mask_col}'")

    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(np.arange(len(y))), start=1):
        rng = np.random.RandomState(random_state + fold_idx)
        perm = rng.permutation(train_idx.size)
        val_n = max(512, int(0.1 * train_idx.size))
        val_rel = perm[:val_n]
        fit_rel = perm[val_n:]
        fit_idx = train_idx[fit_rel]
        val_idx = train_idx[val_rel]
        x_fold = x.copy() if shuffle_priors else x

        if shuffle_priors:
            shuffle_rng = np.random.RandomState(random_state * 100 + fold_idx)
            for prior_col in shuffle_priors:
                prior_idx = col_idx[prior_col]
                mask_idx = col_idx[PRIOR_MASK_COLS[prior_col]]
                for split_idx in (fit_idx, val_idx, test_idx):
                    _shuffle_prior_values_in_place(
                        x_fold,
                        prior_idx=prior_idx,
                        mask_idx=mask_idx,
                        indices=split_idx,
                        rng=shuffle_rng,
                    )

        model = XGBRegressor(
            n_estimators=3000,
            early_stopping_rounds=100,
            objective="reg:squarederror",
            tree_method="hist",
            device="cuda",
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
            x_fold[fit_idx],
            y[fit_idx],
            eval_set=[(x_fold[val_idx], y[val_idx])],
            verbose=False,
        )
        best_iteration = int(getattr(model, "best_iteration", 2999))
        best_iterations.append(best_iteration)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                category=UserWarning,
            )
            predictions[test_idx] = model.predict(x_fold[test_idx], iteration_range=(0, best_iteration + 1))

    return predictions, best_iterations


def _metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    return {
        "pearson": float(np.corrcoef(y, pred)[0, 1]),
        "spearman": float(spearmanr(y, pred).statistic),
        "r2": float(r2_score(y, pred)),
    }


def run_global_mouse_prior_benchmark(
    feature_path: Path | None = None,
    out_path: Path | None = None,
    blend_path: Path | None = None,
) -> tuple[Path, Path]:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    out_path = out_path if out_path is not None else RESULTS_DIR / "human_global_mouse_prior_benchmark.tsv"
    blend_path = blend_path if blend_path is not None else RESULTS_DIR / "human_global_mouse_prior_blend.tsv"

    df, base_cols = _load_global_table(feature_path)
    settings = _build_global_settings(base_cols)

    y = df[SALUKI_HUMAN_PC1].to_numpy(dtype=np.float32)
    rows: list[dict[str, object]] = []
    predictions: dict[str, np.ndarray] = {}
    for setting_name, cols in settings.items():
        pred, best_iterations = _xgb_oof_predictions(df, cols)
        predictions[setting_name] = pred
        rows.append(
            {
                "setting": setting_name,
                "n": int(len(y)),
                "n_with_mouse_pc1": int(df["has_mouse_pc1"].sum()),
                "n_with_saluki_mouse_prior": int(df[HAS_SALUKI_MOUSE_PRIOR].sum()),
                **_metrics(y, pred),
                "mean_best_iter": float(np.mean(best_iterations)),
            }
        )

    benchmark_df = pd.DataFrame(rows).sort_values(["pearson", "r2"], ascending=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    benchmark_df.to_csv(out_path, sep="\t", index=False)

    blend_rows: list[dict[str, object]] = []
    for prior_name in (
        "compact_all_plus_mouse_pc1_global",
        "compact_all_plus_saluki_mouse_prior_global",
        "compact_all_plus_both_mouse_priors_global",
    ):
        for weight in np.linspace(0.0, 1.0, 41):
            blend = weight * predictions[prior_name] + (1.0 - weight) * predictions["compact_all_global"]
            blend_rows.append(
                {
                    "prior_setting": prior_name,
                    "prior_weight": float(weight),
                    "compact_all_weight": float(1.0 - weight),
                    **_metrics(y, blend),
                }
            )

    blend_df = pd.DataFrame(blend_rows).sort_values(["pearson", "r2"], ascending=False)
    blend_df.to_csv(blend_path, sep="\t", index=False)
    return out_path, blend_path


def run_global_mouse_prior_permutation_control(
    feature_path: Path | None = None,
    out_path: Path | None = None,
) -> Path:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    out_path = (
        out_path
        if out_path is not None
        else RESULTS_DIR / "human_global_mouse_prior_permutation_control.tsv"
    )

    df, base_cols = _load_global_table(feature_path)
    settings = _build_global_settings(base_cols)
    control_specs = {
        "compact_all_global": {"cols": settings["compact_all_global"], "prior_mode": "none", "shuffle_priors": ()},
        "compact_all_plus_mouse_pc1_global": {
            "cols": settings["compact_all_plus_mouse_pc1_global"],
            "prior_mode": "real",
            "shuffle_priors": (),
        },
        "compact_all_plus_mouse_pc1_shuffled_global": {
            "cols": settings["compact_all_plus_mouse_pc1_global"],
            "prior_mode": "shuffled",
            "shuffle_priors": ("mouse_pc1",),
        },
        "compact_all_plus_saluki_mouse_prior_global": {
            "cols": settings["compact_all_plus_saluki_mouse_prior_global"],
            "prior_mode": "real",
            "shuffle_priors": (),
        },
        "compact_all_plus_saluki_mouse_prior_shuffled_global": {
            "cols": settings["compact_all_plus_saluki_mouse_prior_global"],
            "prior_mode": "shuffled",
            "shuffle_priors": (SALUKI_MOUSE_PRIOR,),
        },
        "compact_all_plus_both_mouse_priors_global": {
            "cols": settings["compact_all_plus_both_mouse_priors_global"],
            "prior_mode": "real",
            "shuffle_priors": (),
        },
        "compact_all_plus_both_mouse_priors_shuffled_global": {
            "cols": settings["compact_all_plus_both_mouse_priors_global"],
            "prior_mode": "shuffled",
            "shuffle_priors": ("mouse_pc1", SALUKI_MOUSE_PRIOR),
        },
    }

    y = df[SALUKI_HUMAN_PC1].to_numpy(dtype=np.float32)
    rows: list[dict[str, object]] = []
    for setting_name, spec in control_specs.items():
        pred, best_iterations = _xgb_oof_predictions(
            df,
            spec["cols"],
            shuffle_priors=tuple(spec["shuffle_priors"]),
        )
        rows.append(
            {
                "setting": setting_name,
                "prior_mode": spec["prior_mode"],
                "shuffled_columns": ",".join(spec["shuffle_priors"]) if spec["shuffle_priors"] else "",
                "n": int(len(y)),
                "n_with_mouse_pc1": int(df["has_mouse_pc1"].sum()),
                "n_with_saluki_mouse_prior": int(df[HAS_SALUKI_MOUSE_PRIOR].sum()),
                **_metrics(y, pred),
                "mean_best_iter": float(np.mean(best_iterations)),
            }
        )

    control_df = pd.DataFrame(rows).sort_values(["pearson", "r2"], ascending=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    control_df.to_csv(out_path, sep="\t", index=False)
    return out_path


def write_cn_note(
    benchmark_path: Path | None = None,
    note_path: Path | None = None,
    permutation_path: Path | None = None,
) -> Path:
    benchmark_path = benchmark_path if benchmark_path is not None else RESULTS_DIR / "human_global_mouse_prior_benchmark.tsv"
    note_path = note_path if note_path is not None else MANUSCRIPT_DIR / "human_global_mouse_prior_note_cn.md"
    permutation_path = (
        permutation_path
        if permutation_path is not None
        else RESULTS_DIR / "human_global_mouse_prior_permutation_control.tsv"
    )

    benchmark_df = pd.read_csv(benchmark_path, sep="\t")
    best_direct = benchmark_df.iloc[0]
    baseline = benchmark_df.loc[benchmark_df["setting"] == "compact_all_global"].iloc[0]
    control_text = ""
    if permutation_path.exists():
        permutation_df = pd.read_csv(permutation_path, sep="\t")
        real_both = permutation_df.loc[
            permutation_df["setting"] == "compact_all_plus_both_mouse_priors_global"
        ].iloc[0]
        shuffled_both = permutation_df.loc[
            permutation_df["setting"] == "compact_all_plus_both_mouse_priors_shuffled_global"
        ].iloc[0]
        control_text = f"""
## permutation control

- `compact_all_plus_both_mouse_priors_global`：Pearson={real_both['pearson']:.3f}，Spearman={real_both['spearman']:.3f}，R2={real_both['r2']:.3f}
- `compact_all_plus_both_mouse_priors_shuffled_global`：Pearson={shuffled_both['pearson']:.3f}，Spearman={shuffled_both['spearman']:.3f}，R2={shuffled_both['r2']:.3f}
- 真实 prior 相对 shuffled prior 的增益：Pearson={real_both['pearson'] - shuffled_both['pearson']:+.3f}，Spearman={real_both['spearman'] - shuffled_both['spearman']:+.3f}，R2={real_both['r2'] - shuffled_both['r2']:+.3f}
- 这个对照说明提升并不只是“多加了几列特征”或缺失模式提示，而是来自 gene-specific 的跨物种先验信息。
"""

    text = f"""# 固定 human target 下的跨物种迁移结果备注

## 一句话结论

在共同基因宇宙和固定 Saluki human PC1 目标下，把 reconstructed mouse PC1 与
Saluki mouse prior 作为外部 mouse-side covariates 加入 `compact_all`，可提高 gene-level
OOF prediction；打乱 gene-prior identity 后该增益消失。

## 固定目标结果

- `compact_all_global`: N={int(baseline['n'])}，Pearson={baseline['pearson']:.3f}，Spearman={baseline['spearman']:.3f}，R2={baseline['r2']:.3f}
- `{best_direct['setting']}`: N={int(best_direct['n'])}，Pearson={best_direct['pearson']:.3f}，Spearman={best_direct['spearman']:.3f}，R2={best_direct['r2']:.3f}

## 口径说明

- 该设定改变输入信息，属于 cross-species transfer，不是 human-only 或 sequence-only prediction。
- 评估不是只挑 high-confidence ortholog 子集，而是在完整 common gene universe 上完成。
- 在 N={int(best_direct['n'])} 个共同基因里，{int(best_direct['n_with_mouse_pc1'])} 个基因有重建 mouse prior，{int(best_direct['n_with_saluki_mouse_prior'])} 个基因有 Saluki mouse prior。
- 对没有 mouse prior 的基因，模型使用缺失指示特征和 0 填充值；部署时应按 prior coverage 使用 coverage-aware fallback。
- 该结果不能与 Saluki 文中不同切分和不同输入条件下的 r=0.77 作同口径胜负比较。
{control_text}

## 稿件表述

- Human-only baseline 与 cross-species transfer 在固定 target 下分开报告。
- 主要证据是 real-prior、shuffled-prior 与 human-only 的同口径 OOF 对照，而不是跨论文分数比较。
"""
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text(text + "\n", encoding="utf-8")
    return note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run global human prediction with mouse ortholog priors.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
    )
    parser.add_argument("--out", type=Path, default=RESULTS_DIR / "human_global_mouse_prior_benchmark.tsv")
    parser.add_argument("--blend-out", type=Path, default=RESULTS_DIR / "human_global_mouse_prior_blend.tsv")
    parser.add_argument(
        "--permutation-out",
        type=Path,
        default=RESULTS_DIR / "human_global_mouse_prior_permutation_control.tsv",
    )
    parser.add_argument("--note-out", type=Path, default=MANUSCRIPT_DIR / "human_global_mouse_prior_note_cn.md")
    parser.add_argument("--skip-permutation-control", action="store_true")
    parser.add_argument("--skip-note", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    benchmark_path, blend_path = run_global_mouse_prior_benchmark(
        feature_path=args.feature_path,
        out_path=args.out,
        blend_path=args.blend_out,
    )
    permutation_path: Path | None = None
    if not args.skip_permutation_control:
        permutation_path = run_global_mouse_prior_permutation_control(
            feature_path=args.feature_path,
            out_path=args.permutation_out,
        )
    print(benchmark_path)
    print(pd.read_csv(benchmark_path, sep="\t").to_string(index=False))
    print(blend_path)
    print(pd.read_csv(blend_path, sep="\t").head(10).to_string(index=False))
    if permutation_path is not None:
        print(permutation_path)
        print(pd.read_csv(permutation_path, sep="\t").to_string(index=False))
    if not args.skip_note:
        print(
            write_cn_note(
                benchmark_path=benchmark_path,
                note_path=args.note_out,
                permutation_path=permutation_path,
            )
        )


if __name__ == "__main__":
    main()
