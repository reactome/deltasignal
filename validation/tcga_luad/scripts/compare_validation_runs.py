#!/usr/bin/env python3
"""Compare two TCGA LUAD DeltaSignal runs for score and rank stability."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


METRIC = "terminal_mean_activity"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-run", type=Path, required=True)
    parser.add_argument("--candidate-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def load_activity(run_dir: Path, label: str) -> pd.DataFrame:
    path = run_dir / "tables" / "sample_pathway_activity.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing {label} activity table: {path}")
    table = pd.read_csv(path, sep="\t")
    required = {"sample_id", "pathway_id", METRIC}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{label} activity table is missing columns: {sorted(missing)}")
    return table[["sample_id", "pathway_id", METRIC]].rename(
        columns={METRIC: f"{METRIC}_{label}"}
    )


def finite_correlation(func, x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(func(x, y).statistic)


def top_fraction_overlap(x: np.ndarray, y: np.ndarray, fraction: float = 0.25) -> float:
    count = max(1, math.ceil(len(x) * fraction))
    x_top = set(np.argsort(x)[-count:])
    y_top = set(np.argsort(y)[-count:])
    return len(x_top & y_top) / len(x_top | y_top)


def quantile_bin_agreement(x: np.ndarray, y: np.ndarray, bins: int = 4) -> float:
    x_rank = pd.Series(x).rank(method="average", pct=True).to_numpy()
    y_rank = pd.Series(y).rank(method="average", pct=True).to_numpy()
    x_bin = np.minimum((x_rank * bins).astype(int), bins - 1)
    y_bin = np.minimum((y_rank * bins).astype(int), bins - 1)
    return float(np.mean(x_bin == y_bin))


def summarize_pathway(pathway_id: str, group: pd.DataFrame) -> dict[str, object]:
    reference = group[f"{METRIC}_reference"].to_numpy(dtype=float)
    candidate = group[f"{METRIC}_candidate"].to_numpy(dtype=float)
    delta = candidate - reference
    saturation_low = float(np.mean(candidate <= 1.0))
    saturation_high = float(np.mean(candidate >= 99.0))
    spearman = finite_correlation(spearmanr, reference, candidate)
    pearson = finite_correlation(pearsonr, reference, candidate)

    warnings: list[str] = []
    if np.std(candidate) < 1e-6:
        warnings.append("candidate_scores_constant")
    if saturation_low + saturation_high >= 0.10:
        warnings.append("candidate_boundary_saturation_ge_10pct")
    if math.isnan(spearman) or spearman < 0.80:
        warnings.append("rank_stability_below_0.80")

    return {
        "pathway_id": pathway_id,
        "paired_samples": len(group),
        "spearman_rank_correlation": spearman,
        "pearson_correlation": pearson,
        "mean_absolute_score_change": float(np.mean(np.abs(delta))),
        "median_score_change": float(np.median(delta)),
        "reference_q05": float(np.quantile(reference, 0.05)),
        "reference_median": float(np.median(reference)),
        "reference_q95": float(np.quantile(reference, 0.95)),
        "candidate_q05": float(np.quantile(candidate, 0.05)),
        "candidate_median": float(np.median(candidate)),
        "candidate_q95": float(np.quantile(candidate, 0.95)),
        "candidate_std": float(np.std(candidate)),
        "candidate_fraction_le_1": saturation_low,
        "candidate_fraction_ge_99": saturation_high,
        "top_quartile_jaccard": top_fraction_overlap(reference, candidate),
        "rank_quartile_agreement": quantile_bin_agreement(reference, candidate),
        "warnings": ";".join(warnings),
    }


def write_report(summary: pd.DataFrame, output_path: Path) -> None:
    lines = [
        "# DeltaSignal Validation Run Comparison",
        "",
        "This report compares identical sample/pathway pairs between two runs. Score shifts and rank shifts are reported separately because a monotone rescaling can preserve patient ranking while changing the absolute activity scale.",
        "",
        "| Pathway | N | Spearman | MAE | Candidate median | <=1 | >=99 | Top-quartile Jaccard | Warnings |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.pathway_id} | {row.paired_samples} | "
            f"{row.spearman_rank_correlation:.3f} | {row.mean_absolute_score_change:.3f} | "
            f"{row.candidate_median:.3f} | {row.candidate_fraction_le_1:.1%} | "
            f"{row.candidate_fraction_ge_99:.1%} | {row.top_quartile_jaccard:.3f} | "
            f"{row.warnings or 'none'} |"
        )
    lines.extend(
        [
            "",
            "Interpretation guardrails:",
            "",
            "- High Spearman with a large absolute shift indicates mostly monotone recalibration.",
            "- Low Spearman indicates that samples changed order; downstream survival associations can therefore change materially.",
            "- Boundary saturation indicates loss of dynamic range even when convergence succeeds numerically.",
            "- This comparison measures stability, not biological correctness. External perturbation ground truth is required for predictive accuracy.",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference = load_activity(args.reference_run, "reference")
    candidate = load_activity(args.candidate_run, "candidate")
    paired = reference.merge(candidate, on=["sample_id", "pathway_id"], how="inner")
    if paired.empty:
        raise ValueError("The two runs have no shared sample/pathway pairs.")

    summary = pd.DataFrame(
        [summarize_pathway(pathway_id, group) for pathway_id, group in paired.groupby("pathway_id")]
    ).sort_values("pathway_id")
    paired.to_csv(args.output_dir / "paired_scores.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "score_rank_stability.tsv", sep="\t", index=False)
    write_report(summary, args.output_dir / "RUN_COMPARISON.md")
    metadata = {
        "reference_run": str(args.reference_run.resolve()),
        "candidate_run": str(args.candidate_run.resolve()),
        "paired_records": len(paired),
        "pathways": summary["pathway_id"].tolist(),
    }
    (args.output_dir / "comparison_manifest.json").write_text(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
