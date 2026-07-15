from __future__ import annotations

import argparse
import csv
import gzip
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack, load_npz, save_npz
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold

from mrna_half_life_paper.config import INTERIM_DIR, MANUSCRIPT_DIR, RAW_DIR, RESULTS_DIR
from mrna_half_life_paper.saluki_naming import SALUKI_HUMAN_PC1
from mrna_half_life_paper.sequence_baseline import BenchmarkConfig, XGBRegressor, _load_target_tables

try:
    from xgboost.core import XGBoostError
except ImportError:  # pragma: no cover - xgboost is an explicit runtime dependency here
    XGBoostError = RuntimeError


SALUKI_DIR = RAW_DIR / "saluki_datapack" / "human"
BASE_FEATURE_FILE = SALUKI_DIR / "seqFeatWithKmerFreqs.txt.gz"
CWCS_FILE = SALUKI_DIR / "CWCS.txt.gz"
SEQWEAVER_DIR = SALUKI_DIR / "SeqWeaver_predictions"
DEEPRIPE_DIR = SALUKI_DIR / "DeepRiPe_predictions_v83"
CACHE_DIR = INTERIM_DIR / "saluki_datapack" / "human"


@dataclass(frozen=True)
class FeatureBlock:
    name: str
    matrix: csr_matrix
    gene_ids: np.ndarray
    feature_names: list[str]


@dataclass(frozen=True)
class FeatureSetting:
    name: str
    include_cwcs: bool = False
    include_seqweaver: bool = False
    include_deepripe: bool = False
    min_base_nnz: int = 1


DEFAULT_SETTINGS = (
    FeatureSetting(name="base"),
    FeatureSetting(name="base_nnz5", min_base_nnz=5),
    FeatureSetting(name="base_nnz20", min_base_nnz=20),
    FeatureSetting(name="base_cwcs", include_cwcs=True),
    FeatureSetting(name="base_seqweaver", include_seqweaver=True),
    FeatureSetting(name="base_deepripe", include_deepripe=True),
    FeatureSetting(name="base_cwcs_seqweaver", include_cwcs=True, include_seqweaver=True),
    FeatureSetting(name="base_cwcs_deepripe", include_cwcs=True, include_deepripe=True),
    FeatureSetting(name="base_seqweaver_deepripe", include_seqweaver=True, include_deepripe=True),
    FeatureSetting(
        name="base_cwcs_seqweaver_deepripe",
        include_cwcs=True,
        include_seqweaver=True,
        include_deepripe=True,
    ),
    FeatureSetting(
        name="base_cwcs_seqweaver_deepripe_nnz20",
        include_cwcs=True,
        include_seqweaver=True,
        include_deepripe=True,
        min_base_nnz=20,
    ),
)


def _default_benchmark_path() -> Path:
    return RESULTS_DIR / "human_saluki_feature_benchmark.tsv"


def _default_note_path() -> Path:
    return MANUSCRIPT_DIR / "human_saluki_feature_note_cn.md"


def _cache_paths(block_name: str) -> tuple[Path, Path, Path]:
    matrix_path = CACHE_DIR / f"{block_name}.npz"
    gene_path = CACHE_DIR / f"{block_name}.gene_ids.npy"
    feature_path = CACHE_DIR / f"{block_name}.feature_names.txt.gz"
    return matrix_path, gene_path, feature_path


def _write_feature_names(path: Path, feature_names: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for name in feature_names:
            handle.write(f"{name}\n")


def _read_feature_names(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def _save_block_to_cache(block: FeatureBlock) -> None:
    matrix_path, gene_path, feature_path = _cache_paths(block.name)
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    save_npz(matrix_path, block.matrix)
    np.save(gene_path, block.gene_ids)
    _write_feature_names(feature_path, block.feature_names)


def _load_block_from_cache(block_name: str) -> FeatureBlock | None:
    matrix_path, gene_path, feature_path = _cache_paths(block_name)
    if not matrix_path.exists() or not gene_path.exists() or not feature_path.exists():
        return None
    return FeatureBlock(
        name=block_name,
        matrix=load_npz(matrix_path).tocsr(),
        gene_ids=np.load(gene_path, allow_pickle=False),
        feature_names=_read_feature_names(feature_path),
    )


def _parse_dense_table_as_block(
    path: Path,
    block_name: str,
    gene_col: str,
    prefix: str,
    force_rebuild: bool = False,
) -> FeatureBlock:
    cached = None if force_rebuild else _load_block_from_cache(block_name)
    if cached is not None:
        return cached

    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={gene_col: "gene_id"})
    feature_names = [f"{prefix}{col}" for col in df.columns if col != "gene_id"]
    values = df.drop(columns=["gene_id"]).to_numpy(dtype=np.float32, copy=False)
    matrix = csr_matrix(values)
    block = FeatureBlock(
        name=block_name,
        matrix=matrix,
        gene_ids=df["gene_id"].astype(str).to_numpy(dtype=str),
        feature_names=feature_names,
    )
    _save_block_to_cache(block)
    return block


def _parse_region_prediction_block(
    directory: Path,
    block_name: str,
    prefix: str,
    force_rebuild: bool = False,
) -> FeatureBlock:
    cached = None if force_rebuild else _load_block_from_cache(block_name)
    if cached is not None:
        return cached

    region_frames = []
    for region_name, file_name in (("utr5", "5pUTR_avg.txt.gz"), ("orf", "ORF_avg.txt.gz"), ("utr3", "3pUTR_avg.txt.gz")):
        frame = pd.read_csv(directory / file_name, sep="\t")
        frame = frame.rename(columns={frame.columns[0]: "gene_id"})
        renamed = {
            col: f"{prefix}{region_name}__{col}"
            for col in frame.columns
            if col != "gene_id"
        }
        region_frames.append(frame.rename(columns=renamed))

    merged = region_frames[0]
    for frame in region_frames[1:]:
        merged = merged.merge(frame, on="gene_id", how="outer")
    merged = merged.fillna(0.0)
    feature_names = [col for col in merged.columns if col != "gene_id"]
    matrix = csr_matrix(merged[feature_names].to_numpy(dtype=np.float32, copy=False))
    block = FeatureBlock(
        name=block_name,
        matrix=matrix,
        gene_ids=merged["gene_id"].astype(str).to_numpy(dtype=str),
        feature_names=feature_names,
    )
    _save_block_to_cache(block)
    return block


def _parse_base_block(force_rebuild: bool = False) -> tuple[FeatureBlock, pd.DataFrame]:
    cached = None if force_rebuild else _load_block_from_cache("base")
    label_path = CACHE_DIR / "base_halflife.tsv"
    if cached is not None and label_path.exists():
        labels = pd.read_csv(label_path, sep="\t")
        return cached, labels

    gene_ids: list[str] = []
    halflife: list[float] = []
    data_chunks: list[np.ndarray] = []
    index_chunks: list[np.ndarray] = []
    indptr = [0]
    feature_names: list[str] = []

    with gzip.open(BASE_FEATURE_FILE, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        feature_names = header[2:]
        for row in reader:
            gene_ids.append(row[0])
            halflife.append(float(row[1]))
            values = np.asarray(row[2:], dtype=np.float32)
            nonzero = np.flatnonzero(values)
            data_chunks.append(values[nonzero])
            index_chunks.append(nonzero.astype(np.int32, copy=False))
            indptr.append(indptr[-1] + int(nonzero.size))

    data = np.concatenate(data_chunks) if data_chunks else np.empty(0, dtype=np.float32)
    indices = np.concatenate(index_chunks) if index_chunks else np.empty(0, dtype=np.int32)
    matrix = csr_matrix(
        (data, indices, np.asarray(indptr, dtype=np.int64)),
        shape=(len(gene_ids), len(feature_names)),
    )
    block = FeatureBlock(
        name="base",
        matrix=matrix,
        gene_ids=np.asarray(gene_ids, dtype=str),
        feature_names=feature_names,
    )
    _save_block_to_cache(block)

    labels = pd.DataFrame({"gene_id": gene_ids, "saluki_halflife": halflife})
    label_path.parent.mkdir(parents=True, exist_ok=True)
    labels.to_csv(label_path, sep="\t", index=False)
    return block, labels


def _align_block_to_base(block: FeatureBlock, base_gene_ids: np.ndarray) -> FeatureBlock:
    if np.array_equal(block.gene_ids, base_gene_ids):
        return block

    gene_to_row = {gene_id: idx for idx, gene_id in enumerate(block.gene_ids)}
    target_rows: list[int] = []
    source_rows: list[int] = []
    for target_idx, gene_id in enumerate(base_gene_ids):
        source_idx = gene_to_row.get(gene_id)
        if source_idx is None:
            continue
        target_rows.append(target_idx)
        source_rows.append(source_idx)

    if not source_rows:
        aligned_matrix = csr_matrix((len(base_gene_ids), block.matrix.shape[1]), dtype=block.matrix.dtype)
    else:
        subset = block.matrix[source_rows].tocoo()
        target_rows_arr = np.asarray(target_rows, dtype=np.int32)
        aligned_matrix = csr_matrix(
            (subset.data, (target_rows_arr[subset.row], subset.col)),
            shape=(len(base_gene_ids), block.matrix.shape[1]),
        )
    return FeatureBlock(
        name=block.name,
        matrix=aligned_matrix.tocsr(),
        gene_ids=np.asarray(base_gene_ids, dtype=str),
        feature_names=block.feature_names,
    )


def load_feature_blocks(force_rebuild: bool = False) -> tuple[FeatureBlock, FeatureBlock, FeatureBlock, FeatureBlock, pd.DataFrame]:
    base_block, saluki_labels = _parse_base_block(force_rebuild=force_rebuild)
    cwcs_block = _align_block_to_base(
        _parse_dense_table_as_block(CWCS_FILE, "cwcs", "GeneID", "cwcs__", force_rebuild=force_rebuild),
        base_block.gene_ids,
    )
    seqweaver_block = _align_block_to_base(
        _parse_region_prediction_block(SEQWEAVER_DIR, "seqweaver", "seqweaver__", force_rebuild=force_rebuild),
        base_block.gene_ids,
    )
    deepripe_block = _align_block_to_base(
        _parse_region_prediction_block(DEEPRIPE_DIR, "deepripe", "deepripe__", force_rebuild=force_rebuild),
        base_block.gene_ids,
    )
    return base_block, cwcs_block, seqweaver_block, deepripe_block, saluki_labels


def _filter_base_block(base_block: FeatureBlock, min_nnz: int) -> FeatureBlock:
    if min_nnz <= 1:
        return base_block
    nnz_per_feature = np.diff(base_block.matrix.tocsc().indptr)
    keep = nnz_per_feature >= min_nnz
    filtered = base_block.matrix[:, keep]
    feature_names = [name for name, flag in zip(base_block.feature_names, keep, strict=False) if flag]
    return FeatureBlock(
        name=f"{base_block.name}_nnz{min_nnz}",
        matrix=filtered.tocsr(),
        gene_ids=base_block.gene_ids,
        feature_names=feature_names,
    )


def _assemble_setting_matrix(
    setting: FeatureSetting,
    base_block: FeatureBlock,
    cwcs_block: FeatureBlock,
    seqweaver_block: FeatureBlock,
    deepripe_block: FeatureBlock,
) -> tuple[csr_matrix, list[str]]:
    blocks = [_filter_base_block(base_block, setting.min_base_nnz)]
    if setting.include_cwcs:
        blocks.append(cwcs_block)
    if setting.include_seqweaver:
        blocks.append(seqweaver_block)
    if setting.include_deepripe:
        blocks.append(deepripe_block)
    matrix = hstack([block.matrix for block in blocks], format="csr")
    feature_names: list[str] = []
    for block in blocks:
        feature_names.extend(block.feature_names)
    return matrix, feature_names


def _setting_lookup(settings: list[str] | None) -> list[FeatureSetting]:
    if settings is None:
        return list(DEFAULT_SETTINGS)
    lookup = {setting.name: setting for setting in DEFAULT_SETTINGS}
    missing = [name for name in settings if name not in lookup]
    if missing:
        raise ValueError(f"Unknown feature settings: {', '.join(missing)}")
    return [lookup[name] for name in settings]


def run_saluki_feature_benchmark(
    benchmark_path: Path | None = None,
    note_path: Path | None = None,
    targets: list[str] | None = None,
    settings: list[str] | None = None,
    common_only: bool = False,
    force_rebuild: bool = False,
    benchmark_config: BenchmarkConfig = BenchmarkConfig(),
) -> tuple[Path, Path]:
    if XGBRegressor is None:
        raise ImportError("xgboost is required for saluki_feature_benchmark")

    benchmark_path = Path(benchmark_path) if benchmark_path is not None else _default_benchmark_path()
    note_path = Path(note_path) if note_path is not None else _default_note_path()

    base_block, cwcs_block, seqweaver_block, deepripe_block, saluki_labels = load_feature_blocks(
        force_rebuild=force_rebuild
    )
    target_tables = _load_target_tables("human")
    target_tables["saluki_halflife"] = saluki_labels

    if targets is None:
        targets = [SALUKI_HUMAN_PC1]

    selected_settings = _setting_lookup(settings)
    kfold = KFold(
        n_splits=benchmark_config.cv_splits,
        shuffle=True,
        random_state=benchmark_config.random_state,
    )

    rows: list[dict[str, object]] = []
    benchmark_path.parent.mkdir(parents=True, exist_ok=True)
    for setting in selected_settings:
        x_full, feature_names = _assemble_setting_matrix(
            setting=setting,
            base_block=base_block,
            cwcs_block=cwcs_block,
            seqweaver_block=seqweaver_block,
            deepripe_block=deepripe_block,
        )
        feature_gene_df = pd.DataFrame({"gene_id": base_block.gene_ids, "row_idx": np.arange(x_full.shape[0], dtype=np.int32)})

        common_gene_ids: set[str] | None = None
        if common_only:
            common_gene_ids = set(feature_gene_df["gene_id"].astype(str))
            for target_name in targets:
                common_gene_ids &= set(target_tables[target_name]["gene_id"].astype(str))

        for target_name in targets:
            merged = feature_gene_df.merge(target_tables[target_name], on="gene_id", how="inner")
            if common_gene_ids is not None:
                merged = merged.loc[merged["gene_id"].isin(common_gene_ids)].copy()
            row_idx = merged["row_idx"].to_numpy(dtype=np.int32)
            x = x_full[row_idx]
            y = merged[target_name].to_numpy(dtype=np.float32)

            pred = np.full(len(y), np.nan, dtype=np.float32)
            best_iterations: list[int] = []
            try:
                for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(np.arange(len(y))), start=1):
                    rng = np.random.RandomState(benchmark_config.random_state + fold_idx)
                    perm = rng.permutation(train_idx.size)
                    val_n = max(512, int(0.1 * train_idx.size))
                    val_rel_idx = perm[:val_n]
                    fit_rel_idx = perm[val_n:]

                    train_rows = train_idx
                    fit_rows = train_rows[fit_rel_idx]
                    val_rows = train_rows[val_rel_idx]

                    model = XGBRegressor(
                        n_estimators=benchmark_config.xgb_n_estimators,
                        early_stopping_rounds=benchmark_config.xgb_early_stopping_rounds,
                        objective="reg:squarederror",
                        tree_method="hist",
                        device=benchmark_config.xgb_device,
                        max_bin=benchmark_config.xgb_max_bin,
                        max_depth=benchmark_config.xgb_max_depth,
                        learning_rate=benchmark_config.xgb_learning_rate,
                        min_child_weight=benchmark_config.xgb_min_child_weight,
                        subsample=benchmark_config.xgb_subsample,
                        colsample_bytree=benchmark_config.xgb_colsample_bytree,
                        reg_alpha=benchmark_config.xgb_reg_alpha,
                        reg_lambda=benchmark_config.xgb_reg_lambda,
                        random_state=benchmark_config.random_state,
                    )
                    model.fit(
                        x[fit_rows],
                        y[fit_rows],
                        eval_set=[(x[val_rows], y[val_rows])],
                        verbose=False,
                    )
                    best_iteration = int(getattr(model, "best_iteration", benchmark_config.xgb_n_estimators - 1))
                    best_iterations.append(best_iteration)
                    with warnings.catch_warnings():
                        warnings.filterwarnings(
                            "ignore",
                            message=".*Falling back to prediction using DMatrix due to mismatched devices.*",
                            category=UserWarning,
                        )
                        pred[test_idx] = model.predict(x[test_idx], iteration_range=(0, best_iteration + 1))

                rows.append(
                    {
                        "target": target_name,
                        "feature_setting": setting.name,
                        "gene_universe": "common" if common_only else "target_specific",
                        "n": int(len(y)),
                        "n_features": int(len(feature_names)),
                        "pearson": float(np.corrcoef(y, pred)[0, 1]),
                        "spearman": float(spearmanr(y, pred).statistic),
                        "r2": float(r2_score(y, pred)),
                        "mean_best_iter": float(np.mean(best_iterations)),
                        "status": "ok",
                    }
                )
                pd.DataFrame(rows).sort_values(["target", "r2"], ascending=[True, False]).to_csv(
                    benchmark_path,
                    sep="\t",
                    index=False,
                )
                latest = rows[-1]
                print(
                    f"[saluki_feature_benchmark] {target_name} {setting.name} "
                    f"Pearson={latest['pearson']:.4f} R2={latest['r2']:.4f} "
                    f"features={latest['n_features']}",
                    flush=True,
                )
            except XGBoostError as exc:
                rows.append(
                    {
                        "target": target_name,
                        "feature_setting": setting.name,
                        "gene_universe": "common" if common_only else "target_specific",
                        "n": int(len(y)),
                        "n_features": int(len(feature_names)),
                        "status": "xgboost_error",
                        "error": str(exc).splitlines()[0],
                    }
                )
                pd.DataFrame(rows).to_csv(benchmark_path, sep="\t", index=False)
                print(
                    f"[saluki_feature_benchmark] {target_name} {setting.name} failed: {rows[-1]['error']}",
                    flush=True,
                )

    benchmark_df = pd.DataFrame(rows).sort_values(["target", "r2"], ascending=[True, False])
    benchmark_df.to_csv(benchmark_path, sep="\t", index=False)
    write_saluki_feature_note(benchmark_df, note_path)
    return benchmark_path, note_path


def write_saluki_feature_note(benchmark_df: pd.DataFrame, note_path: Path | None = None) -> Path:
    note_path = Path(note_path) if note_path is not None else _default_note_path()
    best = benchmark_df.sort_values("r2", ascending=False).groupby("target", as_index=False).first()
    rows = {row["target"]: row for _, row in best.iterrows()}

    lines = [
        "# 原始 Saluki 特征包复现实验补充",
        "",
        "## 一句话结论",
        "",
        "把原始 `seqFeatWithKmerFreqs` 与 `CWCS`、`SeqWeaver`、`DeepRiPe` 接回 GPU XGBoost 后，",
        "可以直接检验哪些原文特征块真的还能继续抬升人类半衰期预测上限。",
        "",
        "## 当前最好结果",
        "",
    ]

    for target in best["target"].tolist():
        row = rows[target]
        lines.append(
            f"- `{target}`: 最佳设置 `{row['feature_setting']}`，N={int(row['n'])}，"
            f"features={int(row['n_features'])}，Pearson={row['pearson']:.3f}，"
            f"Spearman={row['spearman']:.3f}，R2={row['r2']:.3f}"
        )

    lines.extend(
        [
            "",
            "## 解读",
            "",
            "- 这组实验更接近原始 Saluki 的手工特征路线，因此比我们前面的轻量区域特征更适合做上限追踪。",
            "- 如果 `CWCS/SeqWeaver/DeepRiPe` 中某一块单独带来稳定增益，它就值得进入正式论文消融。",
            "- 如果最佳结果仍然主要来自 `base` 本身，说明后续更值得继续发力的是跨物种先验，而不是盲目堆更多手工特征。",
        ]
    )

    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return note_path


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark the original Saluki feature pack with GPU XGBoost.",
    )
    parser.add_argument(
        "--benchmark-path",
        type=Path,
        default=_default_benchmark_path(),
        help="Where to write the benchmark table.",
    )
    parser.add_argument(
        "--note-path",
        type=Path,
        default=_default_note_path(),
        help="Where to write the Chinese summary note.",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=[SALUKI_HUMAN_PC1],
        help=f"Human targets to benchmark. Defaults to {SALUKI_HUMAN_PC1} only.",
    )
    parser.add_argument(
        "--settings",
        nargs="+",
        default=None,
        help="Optional feature-setting names to run.",
    )
    parser.add_argument(
        "--common-only",
        action="store_true",
        help="Restrict to the common gene universe shared by all requested targets.",
    )
    parser.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Rebuild cached sparse feature blocks from the raw Saluki datapack.",
    )
    parser.add_argument(
        "--xgb-device",
        default="cuda",
        help="XGBoost device string, e.g. cuda or cpu.",
    )
    parser.add_argument(
        "--cv-splits",
        type=int,
        default=5,
        help="Number of cross-validation folds.",
    )
    parser.add_argument(
        "--xgb-n-estimators",
        type=int,
        default=5000,
        help="Maximum number of XGBoost trees.",
    )
    parser.add_argument(
        "--xgb-early-stopping-rounds",
        type=int,
        default=120,
        help="Early stopping patience for XGBoost.",
    )
    parser.add_argument(
        "--xgb-max-depth",
        type=int,
        default=6,
        help="Maximum tree depth.",
    )
    parser.add_argument(
        "--xgb-learning-rate",
        type=float,
        default=0.02,
        help="XGBoost learning rate.",
    )
    parser.add_argument(
        "--xgb-min-child-weight",
        type=int,
        default=4,
        help="XGBoost min_child_weight.",
    )
    parser.add_argument(
        "--xgb-subsample",
        type=float,
        default=0.9,
        help="XGBoost subsample.",
    )
    parser.add_argument(
        "--xgb-colsample-bytree",
        type=float,
        default=0.75,
        help="XGBoost colsample_bytree.",
    )
    parser.add_argument(
        "--xgb-reg-alpha",
        type=float,
        default=0.0,
        help="XGBoost L1 penalty.",
    )
    parser.add_argument(
        "--xgb-reg-lambda",
        type=float,
        default=1.0,
        help="XGBoost L2 penalty.",
    )
    parser.add_argument(
        "--xgb-max-bin",
        type=int,
        default=64,
        help="Number of histogram bins for GPU XGBoost.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_argparser()
    args = parser.parse_args(argv)
    benchmark_config = BenchmarkConfig(
        cv_splits=args.cv_splits,
        include_xgboost=True,
        xgb_device=args.xgb_device,
        xgb_n_estimators=args.xgb_n_estimators,
        xgb_early_stopping_rounds=args.xgb_early_stopping_rounds,
        xgb_max_depth=args.xgb_max_depth,
        xgb_learning_rate=args.xgb_learning_rate,
        xgb_min_child_weight=args.xgb_min_child_weight,
        xgb_subsample=args.xgb_subsample,
        xgb_colsample_bytree=args.xgb_colsample_bytree,
        xgb_reg_alpha=args.xgb_reg_alpha,
        xgb_reg_lambda=args.xgb_reg_lambda,
        xgb_max_bin=args.xgb_max_bin,
    )
    run_saluki_feature_benchmark(
        benchmark_path=args.benchmark_path,
        note_path=args.note_path,
        targets=list(args.targets),
        settings=args.settings,
        common_only=args.common_only,
        force_rebuild=args.force_rebuild,
        benchmark_config=benchmark_config,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
