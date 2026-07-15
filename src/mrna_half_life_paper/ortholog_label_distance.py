from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mrna_half_life_paper.config import MANUSCRIPT_DIR, RESULTS_DIR
from mrna_half_life_paper.ortholog_analysis import ONE2ONE_ORTHOLOG_PATH, extract_mouse_human_one2one
from mrna_half_life_paper.ortholog_regularized_label import (
    _load_saluki_label,
    _load_human_base_label,
    _load_reconstructed_mouse_pc1,
    _mouse_to_human_map,
    _zscore,
    make_ortholog_regularized_label,
)
from mrna_half_life_paper.saluki_naming import SALUKI_HUMAN_PC1_Z, SALUKI_MOUSE_PRIOR_Z


def _metrics(left: pd.Series, right: pd.Series, *, comparison: str, subset: str) -> dict[str, float | int | str]:
    aligned = pd.concat([left.rename("left"), right.rename("right")], axis=1).dropna()
    diff = aligned["left"] - aligned["right"]
    abs_diff = diff.abs()
    return {
        "comparison": comparison,
        "subset": subset,
        "n": int(aligned.shape[0]),
        "pearson": float(aligned["left"].corr(aligned["right"], method="pearson")),
        "spearman": float(aligned["left"].corr(aligned["right"], method="spearman")),
        "mse": float(np.mean(np.square(diff))),
        "rmse": float(np.sqrt(np.mean(np.square(diff)))),
        "mae": float(np.mean(abs_diff)),
        "median_abs_error": float(np.median(abs_diff)),
        "mean_signed_diff": float(np.mean(diff)),
        "sd_signed_diff": float(np.std(diff, ddof=0)),
        "p90_abs_error": float(np.quantile(abs_diff, 0.90)),
        "p95_abs_error": float(np.quantile(abs_diff, 0.95)),
    }


def _comparison_series(frame: pd.DataFrame, left_col: str, right_col: str) -> tuple[pd.Series, pd.Series]:
    return (
        frame.set_index("human_gene_id")[left_col].astype(float),
        frame.set_index("human_gene_id")[right_col].astype(float),
    )


def _write_note(summary: pd.DataFrame, note_path: Path) -> Path:
    def pick(comparison: str, subset: str = "all_one2one") -> pd.Series:
        return summary.loc[(summary["comparison"] == comparison) & (summary["subset"] == subset)].iloc[0]

    pre = pick("human_no_gejman_pc1_vs_reconstructed_mouse_pc1")
    post = pick("orthoreg_lambda0p10_vs_reconstructed_mouse_pc1")
    shift = pick("orthoreg_lambda0p10_vs_human_no_gejman_pc1")
    saluki = pick("saluki_human_pc1_vs_saluki_mouse_prior")
    rmse_reduction = 1.0 - float(post["rmse"]) / float(pre["rmse"])
    mae_reduction = 1.0 - float(post["mae"]) / float(pre["mae"])
    mse_reduction = 1.0 - float(post["mse"]) / float(pre["mse"])

    lines = [
        "# Ortholog Label Distance",
        "",
        "## Purpose",
        "",
        "量化 ortholog-informed target shrinkage 之前 human no-Gejman PC1 与 mouse reconstructed prior 的距离，",
        "并检查 `lambda=0.10` 正则化后人鼠标签距离缩小多少。",
        "",
        "## Key Results",
        "",
        f"- 正则化前：human no-Gejman PC1 vs mouse reconstructed prior，N={int(pre['n'])}，"
        f"Pearson={pre['pearson']:.4f}，Spearman={pre['spearman']:.4f}，"
        f"MSE={pre['mse']:.4f}，RMSE={pre['rmse']:.4f}，MAE={pre['mae']:.4f}。",
        f"- `lambda=0.10` 后：orthoreg target vs mouse reconstructed prior，"
        f"Pearson={post['pearson']:.4f}，Spearman={post['spearman']:.4f}，"
        f"MSE={post['mse']:.4f}，RMSE={post['rmse']:.4f}，MAE={post['mae']:.4f}。",
        f"- 相对正则化前，MSE 下降 {mse_reduction:.1%}，RMSE 下降 {rmse_reduction:.1%}，"
        f"MAE 下降 {mae_reduction:.1%}。",
        f"- `lambda=0.10` 标签相对 human no-Gejman PC1 的改变量较小："
        f"Pearson={shift['pearson']:.4f}，RMSE={shift['rmse']:.4f}，MAE={shift['mae']:.4f}。",
        f"- 作为参考，Saluki human PC1 vs Saluki mouse PC1："
        f"Pearson={saluki['pearson']:.4f}，RMSE={saluki['rmse']:.4f}，MAE={saluki['mae']:.4f}。",
        "",
        "## Interpretation",
        "",
        "- 所有距离指标都在 z-scored PC1 单位下计算，因此 MSE/MAE/RMSE 可直接解释为标准化标签尺度上的差异。",
        "- `lambda=0.10` 没有把 human 标签大幅改成 mouse 标签，而是让 human 标签向 mouse ortholog 信号轻度移动。",
        "- 该结果适合作为 ortholog-informed shrinkage target 的补充防守：正则化前人鼠标签存在明显但非随机的差距，轻度正则化提高相关性并降低距离，同时仍保持接近 human 标签。",
        "",
    ]
    note_path.parent.mkdir(parents=True, exist_ok=True)
    note_path.write_text("\n".join(lines), encoding="utf-8")
    return note_path


def run_ortholog_label_distance(results_root: Path | None = None, lambda_: float = 0.10) -> dict[str, Path]:
    results_root = results_root if results_root is not None else RESULTS_DIR / "ortholog_label_distance"
    results_root.mkdir(parents=True, exist_ok=True)

    human_base = _load_human_base_label()
    reconstructed_mouse_pc1 = _load_reconstructed_mouse_pc1()
    human_saluki = _load_saluki_label("human")
    saluki_mouse_prior = _load_saluki_label("mouse")
    orthoreg = make_ortholog_regularized_label(
        human_label=human_base,
        mouse_label=reconstructed_mouse_pc1,
        lambda_=lambda_,
        saluki_reference=human_saluki,
    ).rename(f"orthoreg_lambda{lambda_:g}")

    extract_mouse_human_one2one()
    mapping = pd.read_csv(ONE2ONE_ORTHOLOG_PATH, sep="\t")
    mapping = mapping[
        mapping["human_gene_id"].isin(human_base.index)
        & mapping["mouse_gene_id"].isin(reconstructed_mouse_pc1.index)
    ].copy()
    mapping["human_no_gejman_pc1_z"] = mapping["human_gene_id"].map(_zscore(human_base))
    mapping["reconstructed_mouse_pc1_z"] = mapping["mouse_gene_id"].map(_zscore(reconstructed_mouse_pc1))
    mapping["orthoreg_lambda0p10_z"] = mapping["human_gene_id"].map(orthoreg)
    mapping[SALUKI_HUMAN_PC1_Z] = mapping["human_gene_id"].map(_zscore(human_saluki))
    mapping[SALUKI_MOUSE_PRIOR_Z] = mapping["mouse_gene_id"].map(_zscore(saluki_mouse_prior))
    mapping["mouse_to_human_map_z"] = mapping["human_gene_id"].map(
        _mouse_to_human_map(human_base.index, reconstructed_mouse_pc1)
    )
    mapping["regularization_shift"] = mapping["orthoreg_lambda0p10_z"] - mapping["human_no_gejman_pc1_z"]
    mapping["pre_regularization_diff"] = (
        mapping["human_no_gejman_pc1_z"] - mapping["reconstructed_mouse_pc1_z"]
    )
    mapping["post_regularization_diff"] = (
        mapping["orthoreg_lambda0p10_z"] - mapping["reconstructed_mouse_pc1_z"]
    )
    mapping = mapping.dropna(
        subset=[
            "human_no_gejman_pc1_z",
            "reconstructed_mouse_pc1_z",
            "orthoreg_lambda0p10_z",
        ]
    )

    pair_path = results_root / "ortholog_label_values.tsv"
    mapping.to_csv(pair_path, sep="\t", index=False)

    comparisons = [
        (
            "human_no_gejman_pc1_vs_reconstructed_mouse_pc1",
            "human_no_gejman_pc1_z",
            "reconstructed_mouse_pc1_z",
        ),
        (
            "orthoreg_lambda0p10_vs_reconstructed_mouse_pc1",
            "orthoreg_lambda0p10_z",
            "reconstructed_mouse_pc1_z",
        ),
        (
            "orthoreg_lambda0p10_vs_human_no_gejman_pc1",
            "orthoreg_lambda0p10_z",
            "human_no_gejman_pc1_z",
        ),
        (
            "saluki_human_pc1_vs_saluki_mouse_prior",
            SALUKI_HUMAN_PC1_Z,
            SALUKI_MOUSE_PRIOR_Z,
        ),
        (
            "orthoreg_lambda0p10_vs_saluki_human_pc1",
            "orthoreg_lambda0p10_z",
            SALUKI_HUMAN_PC1_Z,
        ),
    ]
    rows: list[dict[str, float | int | str]] = []
    subsets = {
        "all_one2one": mapping,
        "high_confidence": mapping.loc[mapping["is_high_confidence"].astype(int) == 1].copy(),
    }
    for subset_name, subset_frame in subsets.items():
        for comparison, left_col, right_col in comparisons:
            if left_col not in subset_frame.columns or right_col not in subset_frame.columns:
                continue
            left, right = _comparison_series(subset_frame, left_col, right_col)
            rows.append(_metrics(left, right, comparison=comparison, subset=subset_name))

    summary = pd.DataFrame(rows)
    summary_path = results_root / "distance_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    note_path = _write_note(summary, MANUSCRIPT_DIR.parent / "docs" / "ortholog_label_distance.md")
    return {"summary": summary_path, "pairs": pair_path, "note": note_path}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute human-mouse ortholog label distance metrics.")
    parser.add_argument("--results-root", type=Path, default=RESULTS_DIR / "ortholog_label_distance")
    parser.add_argument("--lambda", dest="lambda_", type=float, default=0.10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = run_ortholog_label_distance(results_root=args.results_root, lambda_=args.lambda_)
    for path in outputs.values():
        print(path)
        if path.suffix == ".tsv":
            print(pd.read_csv(path, sep="\t").head(20).to_string(index=False))


if __name__ == "__main__":
    main()
