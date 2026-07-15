from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.config import MANUSCRIPT_DIR, PROCESSED_DIR, RESULTS_DIR
from mrna_half_life_paper.saluki_naming import saluki_label_filename

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None


def _load_inputs(feature_path: Path) -> pd.DataFrame:
    feature_df = pd.read_csv(feature_path, sep="\t")
    saluki_human = pd.read_csv(
        PROCESSED_DIR / "real_data" / "human" / saluki_label_filename("human"),
        sep="\t",
    )[["gene_id", "saluki_human_pc1"]]
    ortholog_pairs = pd.read_csv(
        RESULTS_DIR / "ortholog_analysis" / "common_ortholog_pairs.tsv",
        sep="\t",
    )[["human_gene_id", "mouse_pc1", "is_high_confidence"]].rename(columns={"human_gene_id": "gene_id"})
    merged = feature_df.merge(saluki_human, on="gene_id", how="inner").merge(ortholog_pairs, on="gene_id", how="inner")
    return merged


def _run_one_setting(
    df: pd.DataFrame,
    feature_cols: list[str],
    target_col: str,
    *,
    use_mouse_prior: bool,
    random_state: int = 42,
) -> dict[str, float | int | str]:
    if XGBRegressor is None:
        raise ImportError("xgboost is not installed. Please install xgboost to run the GPU ortholog-prior benchmark.")

    cols = feature_cols.copy()
    if use_mouse_prior:
        cols.append("mouse_pc1")

    x = df[cols].to_numpy(dtype=np.float32)
    y = df[target_col].to_numpy(dtype=np.float32)
    kfold = KFold(n_splits=5, shuffle=True, random_state=random_state)

    predictions = np.zeros(len(y), dtype=np.float32)
    best_iterations: list[int] = []
    for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(x), start=1):
        x_train = x[train_idx]
        y_train = y[train_idx]
        x_test = x[test_idx]

        rng = np.random.RandomState(random_state + fold_idx)
        perm = rng.permutation(len(x_train))
        val_n = max(256, int(0.1 * len(x_train)))
        val_idx = perm[:val_n]
        fit_idx = perm[val_n:]

        model = XGBRegressor(
            n_estimators=5000,
            early_stopping_rounds=120,
            objective="reg:squarederror",
            tree_method="hist",
            device="cuda",
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
            x_train[fit_idx],
            y_train[fit_idx],
            eval_set=[(x_train[val_idx], y_train[val_idx])],
            verbose=False,
        )
        best_iteration = int(getattr(model, "best_iteration", 4999))
        best_iterations.append(best_iteration)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                category=UserWarning,
            )
            predictions[test_idx] = model.predict(x_test, iteration_range=(0, best_iteration + 1))

    return {
        "n": int(len(y)),
        "pearson": float(np.corrcoef(y, predictions)[0, 1]),
        "spearman": float(spearmanr(y, predictions).statistic),
        "r2": float(r2_score(y, predictions)),
        "mean_best_iter": float(np.mean(best_iterations)),
    }


def run_cross_species_prior_benchmark(
    feature_path: Path | None = None,
    out_path: Path | None = None,
) -> Path:
    feature_path = feature_path if feature_path is not None else RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz"
    out_path = out_path if out_path is not None else RESULTS_DIR / "human_cross_species_prior_benchmark.tsv"

    merged = _load_inputs(feature_path)
    feature_cols = [
        col
        for col in merged.columns
        if col not in {"gene_id", "transcript_id", "saluki_human_pc1", "mouse_pc1", "is_high_confidence"}
    ]

    settings = {
        "seq_only__ortholog_all": merged.copy(),
        "seq_plus_mouse_pc1__ortholog_all": merged.copy(),
        "seq_only__ortholog_highconf": merged.loc[merged["is_high_confidence"] == 1].copy(),
        "seq_plus_mouse_pc1__ortholog_highconf": merged.loc[merged["is_high_confidence"] == 1].copy(),
    }

    rows: list[dict[str, float | int | str]] = []
    for setting_name, setting_df in settings.items():
        use_mouse_prior = "plus_mouse_pc1" in setting_name
        metrics = _run_one_setting(
            setting_df,
            feature_cols,
            "saluki_human_pc1",
            use_mouse_prior=use_mouse_prior,
        )
        rows.append({"setting": setting_name, **metrics})

    result_df = pd.DataFrame(rows).sort_values("r2", ascending=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(out_path, sep="\t", index=False)
    return out_path


def write_cn_note(
    benchmark_path: Path | None = None,
    note_path: Path | None = None,
) -> Path:
    benchmark_path = benchmark_path if benchmark_path is not None else RESULTS_DIR / "human_cross_species_prior_benchmark.tsv"
    note_path = note_path if note_path is not None else MANUSCRIPT_DIR / "human_cross_species_prior_note_cn.md"

    benchmark_df = pd.read_csv(benchmark_path, sep="\t")
    best = benchmark_df.sort_values("r2", ascending=False).iloc[0]
    baseline = benchmark_df.loc[benchmark_df["setting"] == "seq_only__ortholog_highconf"].iloc[0]
    enhanced = benchmark_df.loc[benchmark_df["setting"] == "seq_plus_mouse_pc1__ortholog_highconf"].iloc[0]

    text = f"""# 跨物种先验增强的人类半衰期预测备注

## 一句话结论

对于存在一对一 ortholog 的人类基因，如果在序列特征之外加入 `mouse_pc1` 作为跨物种先验，
则在 high-confidence ortholog 子集上可把人类 `saluki_human_pc1` 的预测性能提升到
Pearson={enhanced['pearson']:.3f}，Spearman={enhanced['spearman']:.3f}，R2={enhanced['r2']:.3f}。

## 关键对照

- `seq_only__ortholog_highconf`：N={int(baseline['n'])}，Pearson={baseline['pearson']:.3f}，Spearman={baseline['spearman']:.3f}，R2={baseline['r2']:.3f}
- `seq_plus_mouse_pc1__ortholog_highconf`：N={int(enhanced['n'])}，Pearson={enhanced['pearson']:.3f}，Spearman={enhanced['spearman']:.3f}，R2={enhanced['r2']:.3f}

## 解释边界

- `mouse_pc1` 是由 mouse measurements 经 ortholog mapping 得到的外部协变量；human target 不参与其构建。
- 加入该协变量改变了输入信息，因此结果属于 cross-species transfer，不是 human-only baseline。
- 该子集分析用于检查方向一致性；主文的全局结论应以完整 common-gene universe 和 fold-wise permutation control 为准。

## 稿件表述

- 固定 Saluki human PC1 作为目标，分别报告 human-only 与 cross-species-transfer 输入。
- 不把 high-confidence ortholog 子集结果外推为所有 human genes 的全局性能。
- 不与不同切分、不同输入条件下的文献分数作同口径胜负比较。

## 当前子集结果

- 该 benchmark 中数值最高的 setting 为 `{best['setting']}`，N={int(best['n'])}，Pearson={best['pearson']:.3f}，Spearman={best['spearman']:.3f}，R2={best['r2']:.3f}；其解释限定在该 ortholog 子集。
"""
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text(text + "\n", encoding="utf-8")
    return note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark a GPU XGBoost model with an optional mouse ortholog prior.")
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=RESULTS_DIR / "human_sequence_features_regions3mer4mer.tsv.gz",
        help="Human feature table used as the sequence backbone.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=RESULTS_DIR / "human_cross_species_prior_benchmark.tsv",
        help="Where to write the benchmark TSV.",
    )
    parser.add_argument(
        "--note-out",
        type=Path,
        default=MANUSCRIPT_DIR / "human_cross_species_prior_note_cn.md",
        help="Where to write the Chinese manuscript note.",
    )
    parser.add_argument("--skip-note", action="store_true", help="Only write the TSV benchmark.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    benchmark_path = run_cross_species_prior_benchmark(
        feature_path=args.feature_path,
        out_path=args.out,
    )
    benchmark_df = pd.read_csv(benchmark_path, sep="\t")
    print(benchmark_path)
    print(benchmark_df.to_string(index=False))
    if not args.skip_note:
        note_path = write_cn_note(benchmark_path=benchmark_path, note_path=args.note_out)
        print(note_path)


if __name__ == "__main__":
    main()
