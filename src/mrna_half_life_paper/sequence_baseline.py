from __future__ import annotations

import argparse
import gzip
import itertools
import math
import re
import urllib.request
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBRegressor
except ImportError:  # pragma: no cover - optional dependency
    XGBRegressor = None

from mrna_half_life_paper.config import MANUSCRIPT_DIR, PROCESSED_DIR, RAW_DIR, RESULTS_DIR
from mrna_half_life_paper.real_data import extract_real_data
from mrna_half_life_paper.saluki_naming import (
    SALUKI_HUMAN_PC1,
    saluki_label_column,
    saluki_label_filename,
)


ENSEMBL_RELEASE = 115
ATTRIBUTE_RE = re.compile(r'(\S+) "([^"]+)"')
TAG_RE = re.compile(r'tag "([^"]+)"')
BASES = "ACGT"
NONSTOP_CODONS = [codon for codon in map("".join, itertools.product(BASES, repeat=3)) if codon not in {"TAA", "TAG", "TGA"}]
TRIMERS = ["".join(kmer) for kmer in itertools.product(BASES, repeat=3)]
TETRAMERS = ["".join(kmer) for kmer in itertools.product(BASES, repeat=4)]

ENSEMBL_REFERENCE_URLS = {
    "human": {
        "gtf": (
            "https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.115.gtf.gz"
        ),
        "cdna": (
            "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/"
            "Homo_sapiens.GRCh38.cdna.all.fa.gz"
        ),
        "cds": (
            "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cds/"
            "Homo_sapiens.GRCh38.cds.all.fa.gz"
        ),
    },
    "mouse": {
        "gtf": (
            "https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/"
            "Mus_musculus.GRCm39.115.gtf.gz"
        ),
        "cdna": (
            "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cdna/"
            "Mus_musculus.GRCm39.cdna.all.fa.gz"
        ),
        "cds": (
            "https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/cds/"
            "Mus_musculus.GRCm39.cds.all.fa.gz"
        ),
    },
}


@dataclass(frozen=True)
class FeatureConfig:
    min_region_length: int = 10
    include_log_lengths: bool = True
    include_region_trimers: bool = True
    include_utr3_tetramers: bool = True


@dataclass(frozen=True)
class BenchmarkConfig:
    cv_splits: int = 5
    random_state: int = 42
    ridge_alpha_grid: tuple[float, ...] = tuple(np.logspace(-3, 3, 15))
    hgb_learning_rate: float = 0.05
    hgb_max_depth: int = 6
    hgb_max_iter: int = 300
    hgb_min_samples_leaf: int = 20
    xt_estimators: int = 400
    xt_min_samples_leaf: int = 3
    ensemble_hgb_weight: float = 0.8
    include_xgboost: bool = True
    xgb_device: str = "cuda"
    xgb_n_estimators: int = 5000
    xgb_early_stopping_rounds: int = 120
    xgb_max_depth: int = 6
    xgb_learning_rate: float = 0.02
    xgb_min_child_weight: int = 4
    xgb_subsample: float = 0.9
    xgb_colsample_bytree: float = 0.75
    xgb_reg_alpha: float = 0.0
    xgb_reg_lambda: float = 1.0
    xgb_max_bin: int = 256


@dataclass(frozen=True)
class ReferencePaths:
    gtf: Path
    cdna: Path
    cds: Path


def _default_reference_dir(species: str, release: int = ENSEMBL_RELEASE) -> Path:
    return RAW_DIR / "ensembl" / f"{species}_release{release}"


def _default_feature_path(species: str) -> Path:
    return RESULTS_DIR / f"{species}_sequence_features_regions3mer4mer.tsv.gz"


def _default_benchmark_path(species: str) -> Path:
    return RESULTS_DIR / f"{species}_sequence_prediction_benchmark.tsv"


def _default_note_path(species: str) -> Path:
    return MANUSCRIPT_DIR / f"{species}_sequence_prediction_note_cn.md"


def _download_file(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return
    urllib.request.urlretrieve(url, dest)


def ensure_reference_files(
    species: str,
    release: int = ENSEMBL_RELEASE,
    reference_dir: Path | None = None,
    force_download: bool = False,
) -> ReferencePaths:
    if species not in ENSEMBL_REFERENCE_URLS:
        raise ValueError(f"Unsupported species: {species}")

    reference_dir = reference_dir if reference_dir is not None else _default_reference_dir(species, release)
    reference_dir.mkdir(parents=True, exist_ok=True)

    file_names = {
        "human": {
            "gtf": "Homo_sapiens.GRCh38.115.gtf.gz",
            "cdna": "Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "cds": "Homo_sapiens.GRCh38.cds.all.fa.gz",
        },
        "mouse": {
            "gtf": "Mus_musculus.GRCm39.115.gtf.gz",
            "cdna": "Mus_musculus.GRCm39.cdna.all.fa.gz",
            "cds": "Mus_musculus.GRCm39.cds.all.fa.gz",
        },
    }[species]

    gtf_path = reference_dir / file_names["gtf"]
    cdna_path = reference_dir / file_names["cdna"]
    cds_path = reference_dir / file_names["cds"]

    if force_download:
        for path in (gtf_path, cdna_path, cds_path):
            if path.exists():
                path.unlink()

    urls = ENSEMBL_REFERENCE_URLS[species]
    if not gtf_path.exists():
        _download_file(urls["gtf"], gtf_path)
    if not cdna_path.exists():
        _download_file(urls["cdna"], cdna_path)
    if not cds_path.exists():
        _download_file(urls["cds"], cds_path)

    return ReferencePaths(gtf=gtf_path, cdna=cdna_path, cds=cds_path)


def _extract_attrs(attr_field: str) -> dict[str, str]:
    return dict(ATTRIBUTE_RE.findall(attr_field))


def _read_fasta_subset(path: Path, wanted_ids: set[str]) -> dict[str, str]:
    seqs: dict[str, str] = {}
    with gzip.open(path, "rt") as handle:
        current_id = ""
        keep = False
        chunks: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if keep and current_id:
                    seqs[current_id] = "".join(chunks).upper()
                current_id = line[1:].split()[0].split(".")[0]
                keep = current_id in wanted_ids
                chunks = []
                continue
            if keep:
                chunks.append(line.strip())
        if keep and current_id:
            seqs[current_id] = "".join(chunks).upper()
    return seqs


def _choose_representative_transcripts(gtf_path: Path, gene_ids: set[str]) -> dict[str, dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 9 or parts[2] != "transcript":
                continue
            attrs = _extract_attrs(parts[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            transcript_id = attrs.get("transcript_id", "").split(".")[0]
            if gene_id not in gene_ids or not transcript_id:
                continue

            tags = set(TAG_RE.findall(parts[8]))
            biotype = attrs.get("transcript_biotype", attrs.get("transcript_type", ""))
            tsl_raw = attrs.get("transcript_support_level", "999").split()[0]
            try:
                tsl = int(tsl_raw)
            except ValueError:
                tsl = 999

            transcript_length = int(parts[4]) - int(parts[3]) + 1
            rank = (
                "MANE_Select" not in tags,
                "Ensembl_canonical" not in tags,
                "gencode_primary" not in tags,
                "gencode_basic" not in tags,
                biotype != "protein_coding",
                tsl,
                -transcript_length,
                transcript_id,
            )
            previous = transcripts.get(gene_id)
            if previous is None or rank < previous["rank"]:
                transcripts[gene_id] = {
                    "rank": rank,
                    "transcript_id": transcript_id,
                    "transcript_length": transcript_length,
                    "biotype": biotype,
                }
    return transcripts


def _count_cds_exons(gtf_path: Path, transcript_ids: set[str]) -> dict[str, int]:
    exon_counts = {transcript_id: 0 for transcript_id in transcript_ids}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = _extract_attrs(parts[8])
            transcript_id = attrs.get("transcript_id", "").split(".")[0]
            if transcript_id in exon_counts:
                exon_counts[transcript_id] += 1
    return exon_counts


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _base_fraction(seq: str, base: str) -> float:
    if not seq:
        return 0.0
    return seq.count(base) / len(seq)


def _count_kmers(seq: str, k: int) -> tuple[Counter[str], int]:
    counts: Counter[str] = Counter()
    valid = 0
    for idx in range(len(seq) - k + 1):
        kmer = seq[idx : idx + k]
        if set(kmer) <= set(BASES):
            counts[kmer] += 1
            valid += 1
    return counts, valid


def _load_target_tables(species: str) -> dict[str, pd.DataFrame]:
    extract_real_data()
    if species == "human":
        full_path = RESULTS_DIR / "real_human_moesm2" / "consensus_labels.tsv"
        pruned_path = RESULTS_DIR / "real_human_moesm2_no_gejman" / "consensus_labels.tsv"
        saluki_path = PROCESSED_DIR / "real_data" / "human" / saluki_label_filename("human")

        if not full_path.exists():
            raise FileNotFoundError(f"Missing {full_path}. Run scripts/run_real_data_workflow.sh first.")
        if not pruned_path.exists():
            raise FileNotFoundError(
                f"Missing {pruned_path}. Run scripts/run_study_influence_analysis.sh --run-pruned-human first."
            )

        return {
            "human_full_pc1": pd.read_csv(full_path, sep="\t")[["gene_id", "pc1"]].rename(columns={"pc1": "human_full_pc1"}),
            "human_pruned_pc1": pd.read_csv(pruned_path, sep="\t")[["gene_id", "pc1"]].rename(columns={"pc1": "human_pruned_pc1"}),
            SALUKI_HUMAN_PC1: pd.read_csv(saluki_path, sep="\t")[
                ["gene_id", saluki_label_column("human")]
            ].rename(
                columns={saluki_label_column("human"): SALUKI_HUMAN_PC1}
            ),
        }

    if species == "mouse":
        full_path = RESULTS_DIR / "real_mouse_moesm2" / "consensus_labels.tsv"
        saluki_path = PROCESSED_DIR / "real_data" / "mouse" / saluki_label_filename("mouse")
        if not full_path.exists():
            raise FileNotFoundError(f"Missing {full_path}. Run scripts/run_real_data_workflow.sh first.")
        return {
            "mouse_full_pc1": pd.read_csv(full_path, sep="\t")[["gene_id", "pc1"]].rename(columns={"pc1": "mouse_full_pc1"}),
            "saluki_mouse_prior": pd.read_csv(saluki_path, sep="\t")[
                ["gene_id", saluki_label_column("mouse")]
            ].rename(
                columns={saluki_label_column("mouse"): "saluki_mouse_prior"}
            ),
        }

    raise ValueError(f"Unsupported species: {species}")


def build_sequence_feature_table(
    species: str = "human",
    feature_path: Path | None = None,
    release: int = ENSEMBL_RELEASE,
    reference_dir: Path | None = None,
    feature_config: FeatureConfig = FeatureConfig(),
    force_download: bool = False,
) -> Path:
    feature_path = feature_path if feature_path is not None else _default_feature_path(species)
    feature_path.parent.mkdir(parents=True, exist_ok=True)

    target_tables = _load_target_tables(species)
    gene_ids = set()
    for table in target_tables.values():
        gene_ids.update(table["gene_id"].astype(str))

    references = ensure_reference_files(
        species=species,
        release=release,
        reference_dir=reference_dir,
        force_download=force_download,
    )

    transcript_meta = _choose_representative_transcripts(references.gtf, gene_ids)
    transcript_ids = {meta["transcript_id"] for meta in transcript_meta.values()}
    exon_counts = _count_cds_exons(references.gtf, transcript_ids)
    cdna_sequences = _read_fasta_subset(references.cdna, transcript_ids)
    cds_sequences = _read_fasta_subset(references.cds, transcript_ids)

    rows: list[dict[str, object]] = []
    for gene_id, meta in transcript_meta.items():
        transcript_id = str(meta["transcript_id"])
        transcript_seq = cdna_sequences.get(transcript_id)
        cds_seq = cds_sequences.get(transcript_id)
        if not transcript_seq or not cds_seq:
            continue

        cds_start = transcript_seq.find(cds_seq)
        if cds_start < 0:
            continue

        utr5_seq = transcript_seq[:cds_start]
        utr3_seq = transcript_seq[cds_start + len(cds_seq) :]
        if min(len(utr5_seq), len(utr3_seq), len(cds_seq)) < feature_config.min_region_length:
            continue

        intron_length = max(int(meta["transcript_length"]) - len(transcript_seq), 0)
        orf_exon_density = exon_counts.get(transcript_id, 0) * 1000.0 / len(cds_seq)
        row: dict[str, object] = {
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "UTR5LEN": len(utr5_seq),
            "CDSLEN": len(cds_seq),
            "INTRONLEN": intron_length,
            "UTR3LEN": len(utr3_seq),
            "UTR5GC": _gc_fraction(utr5_seq),
            "CDSGC": _gc_fraction(cds_seq),
            "UTR3GC": _gc_fraction(utr3_seq),
            "ORFEXONDENSITY": orf_exon_density,
            "UTR3_A_frac": _base_fraction(utr3_seq, "A"),
            "UTR3_U_frac": _base_fraction(utr3_seq, "T"),
            "UTR3_ARE_AUUUA": utr3_seq.count("ATTTA"),
            "UTR3_polyAATAAA": utr3_seq.count("AATAAA"),
        }
        if feature_config.include_log_lengths:
            row.update(
                {
                    "LOG_UTR5LEN": math.log10(len(utr5_seq) + 0.1),
                    "LOG_CDSLEN": math.log10(len(cds_seq) + 0.1),
                    "LOG_INTRONLEN": math.log10(intron_length + 0.1),
                    "LOG_UTR3LEN": math.log10(len(utr3_seq) + 0.1),
                    "LOG_ORFEXONDENSITY": math.log10(orf_exon_density + 0.1),
                }
            )

        trimmed_orf = cds_seq[: (len(cds_seq) // 3) * 3]
        codon_counts = Counter(
            trimmed_orf[idx : idx + 3]
            for idx in range(0, len(trimmed_orf), 3)
            if set(trimmed_orf[idx : idx + 3]) <= set(BASES)
        )
        codon_total = max(sum(codon_counts.values()), 1)
        for codon in NONSTOP_CODONS:
            row[f"CODON_{codon}"] = codon_counts[codon] / codon_total

        if feature_config.include_region_trimers:
            for prefix, seq in (("UTR5_3MER", utr5_seq), ("ORF_3MER", cds_seq), ("UTR3_3MER", utr3_seq)):
                counts, total = _count_kmers(seq, 3)
                total = max(total, 1)
                for kmer in TRIMERS:
                    row[f"{prefix}_{kmer}"] = counts[kmer] / total

        if feature_config.include_utr3_tetramers:
            counts, total = _count_kmers(utr3_seq, 4)
            total = max(total, 1)
            for kmer in TETRAMERS:
                row[f"UTR3_4MER_{kmer}"] = counts[kmer] / total

        rows.append(row)

    feature_df = pd.DataFrame(rows)
    feature_df.to_csv(feature_path, sep="\t", index=False, compression="gzip")
    return feature_path


def run_prediction_benchmark(
    species: str = "human",
    feature_path: Path | None = None,
    benchmark_path: Path | None = None,
    targets: list[str] | None = None,
    common_only: bool = False,
    benchmark_config: BenchmarkConfig = BenchmarkConfig(),
) -> Path:
    feature_path = feature_path if feature_path is not None else _default_feature_path(species)
    benchmark_path = benchmark_path if benchmark_path is not None else _default_benchmark_path(species)

    feature_df = pd.read_csv(feature_path, sep="\t")
    target_tables = _load_target_tables(species)
    if targets is None:
        targets = list(target_tables.keys())

    common_gene_ids: set[str] | None = None
    if common_only:
        common_gene_ids = set(feature_df["gene_id"].astype(str))
        for target_name in targets:
            common_gene_ids &= set(target_tables[target_name]["gene_id"].astype(str))

    kfold = KFold(n_splits=benchmark_config.cv_splits, shuffle=True, random_state=benchmark_config.random_state)
    ridge = make_pipeline(StandardScaler(), RidgeCV(alphas=np.asarray(benchmark_config.ridge_alpha_grid, dtype=float)))
    extra_trees = ExtraTreesRegressor(
        n_estimators=benchmark_config.xt_estimators,
        min_samples_leaf=benchmark_config.xt_min_samples_leaf,
        random_state=benchmark_config.random_state,
        n_jobs=-1,
    )

    rows: list[dict[str, object]] = []
    for target_name in targets:
        if target_name not in target_tables:
            raise ValueError(f"Unknown target {target_name!r} for species {species}")
        merged = feature_df.merge(target_tables[target_name], on="gene_id", how="inner")
        if common_gene_ids is not None:
            merged = merged.loc[merged["gene_id"].isin(common_gene_ids)].copy()
        feature_cols = [col for col in merged.columns if col not in {"gene_id", "transcript_id", target_name}]
        x = merged[feature_cols].to_numpy(dtype=float)
        y = merged[target_name].to_numpy(dtype=float)

        predictions: dict[str, np.ndarray] = {}
        for model_name, model in (
            ("ridge", ridge),
            (
                "hgb",
                HistGradientBoostingRegressor(
                    loss="squared_error",
                    learning_rate=benchmark_config.hgb_learning_rate,
                    max_depth=benchmark_config.hgb_max_depth,
                    max_iter=benchmark_config.hgb_max_iter,
                    min_samples_leaf=benchmark_config.hgb_min_samples_leaf,
                    random_state=benchmark_config.random_state,
                ),
            ),
            ("xt", extra_trees),
        ):
            pred = np.full(len(y), np.nan)
            for train_idx, test_idx in kfold.split(x):
                model.fit(x[train_idx], y[train_idx])
                pred[test_idx] = model.predict(x[test_idx])
            predictions[model_name] = pred
            rows.append(
                {
                    "target": target_name,
                    "model": model_name,
                    "gene_universe": "common" if common_only else "target_specific",
                    "n": int(len(y)),
                    "pearson": float(np.corrcoef(y, pred)[0, 1]),
                    "spearman": float(spearmanr(y, pred).statistic),
                    "r2": float(r2_score(y, pred)),
                }
            )

        if benchmark_config.include_xgboost and XGBRegressor is not None:
            pred = np.full(len(y), np.nan)
            best_iterations: list[int] = []
            for fold_idx, (train_idx, test_idx) in enumerate(kfold.split(x), start=1):
                x_train = x[train_idx]
                y_train = y[train_idx]
                x_test = x[test_idx]

                rng = np.random.RandomState(benchmark_config.random_state + fold_idx)
                perm = rng.permutation(len(x_train))
                val_n = max(512, int(0.1 * len(x_train)))
                val_idx = perm[:val_n]
                fit_idx = perm[val_n:]

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
                    x_train[fit_idx],
                    y_train[fit_idx],
                    eval_set=[(x_train[val_idx], y_train[val_idx])],
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
                    pred[test_idx] = model.predict(x_test, iteration_range=(0, best_iteration + 1))

            predictions["xgb_cuda"] = pred
            rows.append(
                {
                    "target": target_name,
                    "model": "xgb_cuda" if benchmark_config.xgb_device == "cuda" else f"xgb_{benchmark_config.xgb_device}",
                    "gene_universe": "common" if common_only else "target_specific",
                    "n": int(len(y)),
                    "pearson": float(np.corrcoef(y, pred)[0, 1]),
                    "spearman": float(spearmanr(y, pred).statistic),
                    "r2": float(r2_score(y, pred)),
                    "mean_best_iter": float(np.mean(best_iterations)),
                }
            )

        ensemble = (
            benchmark_config.ensemble_hgb_weight * predictions["hgb"]
            + (1.0 - benchmark_config.ensemble_hgb_weight) * predictions["ridge"]
        )
        rows.append(
            {
                "target": target_name,
                "model": f"hgb_ridge_{benchmark_config.ensemble_hgb_weight:.1f}_{1.0 - benchmark_config.ensemble_hgb_weight:.1f}",
                "gene_universe": "common" if common_only else "target_specific",
                "n": int(len(y)),
                "pearson": float(np.corrcoef(y, ensemble)[0, 1]),
                "spearman": float(spearmanr(y, ensemble).statistic),
                "r2": float(r2_score(y, ensemble)),
            }
        )

    benchmark_df = pd.DataFrame(rows).sort_values(["target", "r2"], ascending=[True, False])
    benchmark_path.parent.mkdir(parents=True, exist_ok=True)
    benchmark_df.to_csv(benchmark_path, sep="\t", index=False)
    return benchmark_path


def write_cn_prediction_note(
    species: str = "human",
    benchmark_path: Path | None = None,
    note_path: Path | None = None,
) -> Path:
    benchmark_path = Path(benchmark_path) if benchmark_path is not None else _default_benchmark_path(species)
    note_path = Path(note_path) if note_path is not None else _default_note_path(species)

    benchmark_df = pd.read_csv(benchmark_path, sep="\t")
    best = benchmark_df.sort_values("r2", ascending=False).groupby("target", as_index=False).first()
    by_target = {row["target"]: row for _, row in best.iterrows()}

    lines = [
        "# 论文思路补充：标签清洗后的人类 mRNA 半衰期序列预测基线",
        "",
        "## 一句话结论",
        "",
        "在代表转录本、区域长度/GC、codon 频率、5'UTR/ORF/3'UTR 3-mer 和 3'UTR 4-mer 特征的基础上，",
        "人类半衰期标签的可预测性会随着标签质量提升而同步提高。",
        "",
        "## 当前最好结果",
        "",
    ]

    for target in best["target"].tolist():
        row = by_target[target]
        universe = row.get("gene_universe", "")
        universe_text = "" if not universe else f"，gene universe={universe}"
        lines.append(
            f"- `{target}`: 最佳模型 `{row['model']}`，N={int(row['n'])}，Pearson={row['pearson']:.3f}，"
            f"Spearman={row['spearman']:.3f}，R2={row['r2']:.3f}{universe_text}"
        )

    if {"human_full_pc1", "human_pruned_pc1"} <= set(by_target):
        full_row = by_target["human_full_pc1"]
        pruned_row = by_target["human_pruned_pc1"]
        lines.extend(
            [
                "",
                "## 这组结果说明了什么",
                "",
                f"- 从全量 human 标签切换到去除 `Gejman` 后的标签，最佳基线的 Pearson 从 {full_row['pearson']:.3f} "
                f"提升到 {pruned_row['pearson']:.3f}，R2 从 {full_row['r2']:.3f} 提升到 {pruned_row['r2']:.3f}。",
                "- 该结果表明，study-level 标签处理会改变 matched human-only feature setting 下的可预测性。",
                "- 标签定义敏感性与固定 Saluki human PC1 benchmark 应分开报告，避免把不同 target 的分数放在同一排行榜中。",
            ]
        )

    if any(str(row["model"]).startswith("xgb_") for row in by_target.values()):
        lines.extend(
            [
                "",
                "## GPU 训练备注",
                "",
                "- XGBoost/CUDA 用于统一的 gene-level OOF benchmark；CPU/GPU 设备差异不改变任务定义。",
                "- 公开结果应同时记录固定参数、fold 划分和 random state，而不依赖特定显卡型号解释。",
            ]
        )

    lines.extend(
        [
            "",
            "## 解释边界",
            "",
            "- Study influence、label-definition sensitivity 和 fixed-target human-only prediction 是相关但不同的证据层。",
            "- 主稿只把同一 target、同一 gene universe 和同一 OOF 口径下的模型作为直接比较。",
        ]
    )

    if species == "human":
        tuning_lines: list[str] = []
        hgb_path = RESULTS_DIR / "human_prediction_hgb_tuning.tsv"
        xgb_path = RESULTS_DIR / "human_xgboost_gpu_saluki_refine.tsv"

        if xgb_path.exists():
            xgb_df = pd.read_csv(xgb_path, sep="\t")
            xgb_best = xgb_df.sort_values("r2", ascending=False).iloc[0]
            tuning_lines.extend(
                [
                    f"- 在 `{SALUKI_HUMAN_PC1}` 上，所测 GPU XGBoost 组合中的最高值为 "
                    f"Pearson={xgb_best['pearson']:.3f}，Spearman={xgb_best['spearman']:.3f}，R2={xgb_best['r2']:.3f}。",
                    f"- 这组 GPU 参数对应 `learning_rate={xgb_best['learning_rate']:.2f}`，"
                    f"`max_depth={int(xgb_best['max_depth'])}`，"
                    f"`min_child_weight={int(xgb_best['min_child_weight'])}`，"
                    f"`subsample={xgb_best['subsample']:.2f}`，"
                    f"`colsample_bytree={xgb_best['colsample_bytree']:.2f}`。",
                ]
            )

        ensemble_plus_path = RESULTS_DIR / "human_gpu_ensemble_plus_5mer.tsv"
        if ensemble_plus_path.exists():
            ensemble_df = pd.read_csv(ensemble_plus_path, sep="\t")
            ensemble_best = ensemble_df.sort_values("r2", ascending=False).iloc[0]
            tuning_lines.extend(
                [
                    f"- 探索性 OOF 线性组合（基础 `xgb` + `5mer xgb` + `hgb` + `ridge`）得到 "
                    f"Pearson={ensemble_best['pearson']:.3f}，"
                    f"Spearman={ensemble_best['spearman']:.3f}，R2={ensemble_best['r2']:.3f}。",
                    "- 该线性组合只作为 sensitivity analysis，不作为主稿的单模型结果或跨论文比较依据。",
                ]
            )

        if hgb_path.exists():
            hgb_df = pd.read_csv(hgb_path, sep="\t")
            hgb_best = hgb_df.sort_values("ens_r2", ascending=False).iloc[0]
            tuning_lines.extend(
                [
                    f"- 作为对照，局部调参后的 HGB + Ridge 最高可到 "
                    f"Pearson={hgb_best['ens_pearson']:.3f}，Spearman={hgb_best['ens_spearman']:.3f}，"
                    f"R2={hgb_best['ens_r2']:.3f}。",
                ]
            )

        if tuning_lines:
            lines.extend(["", "## 补充参数与敏感性结果", ""])
            lines.extend(tuning_lines)

    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build sequence features and benchmark lightweight half-life predictors.")
    parser.add_argument("--species", choices=["human", "mouse"], default="human", help="Species to analyze.")
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=None,
        help="Directory used to store Ensembl GTF/cDNA/CDS reference files.",
    )
    parser.add_argument(
        "--feature-path",
        type=Path,
        default=None,
        help="Output path for the compressed feature table. Defaults to results/<species>_sequence_features_regions3mer4mer.tsv.gz.",
    )
    parser.add_argument(
        "--benchmark-out",
        type=Path,
        default=None,
        help="Output path for the benchmark TSV. Defaults to results/<species>_sequence_prediction_benchmark.tsv.",
    )
    parser.add_argument(
        "--note-out",
        type=Path,
        default=None,
        help="Output path for the Chinese prediction note.",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help=f"Targets to benchmark. Human defaults to human_full_pc1 human_pruned_pc1 {SALUKI_HUMAN_PC1}.",
    )
    parser.add_argument("--skip-feature-build", action="store_true", help="Reuse an existing feature table instead of rebuilding it.")
    parser.add_argument("--skip-benchmark", action="store_true", help="Only build the feature table.")
    parser.add_argument("--skip-note", action="store_true", help="Do not write the manuscript note.")
    parser.add_argument(
        "--common-only",
        action="store_true",
        help="Benchmark all requested targets on the shared gene universe instead of each target's own coverage.",
    )
    parser.add_argument(
        "--disable-xgboost",
        action="store_true",
        help="Skip the optional XGBoost benchmark even if xgboost is installed.",
    )
    parser.add_argument(
        "--xgb-device",
        choices=["cuda", "cpu"],
        default="cuda",
        help="Device used by the optional XGBoost benchmark.",
    )
    parser.add_argument("--force-download", action="store_true", help="Re-download Ensembl reference files even if they already exist.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    feature_path = args.feature_path if args.feature_path is not None else _default_feature_path(args.species)
    benchmark_path = args.benchmark_out if args.benchmark_out is not None else _default_benchmark_path(args.species)
    note_path = args.note_out if args.note_out is not None else _default_note_path(args.species)

    if not args.skip_feature_build:
        built_feature_path = build_sequence_feature_table(
            species=args.species,
            feature_path=feature_path,
            reference_dir=args.reference_dir,
            force_download=args.force_download,
        )
        print(built_feature_path)
    else:
        print(feature_path)

    if not args.skip_benchmark:
        built_benchmark_path = run_prediction_benchmark(
            species=args.species,
            feature_path=feature_path,
            benchmark_path=benchmark_path,
            targets=args.targets,
            common_only=args.common_only,
            benchmark_config=BenchmarkConfig(
                include_xgboost=not args.disable_xgboost,
                xgb_device=args.xgb_device,
            ),
        )
        benchmark_df = pd.read_csv(built_benchmark_path, sep="\t")
        print(built_benchmark_path)
        print(benchmark_df.to_string(index=False))

        if not args.skip_note:
            written_note = write_cn_prediction_note(
                species=args.species,
                benchmark_path=built_benchmark_path,
                note_path=note_path,
            )
            print(written_note)


if __name__ == "__main__":
    main()
