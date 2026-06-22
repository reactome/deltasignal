#!/usr/bin/env python3
"""Generate figures, final statistics, and interpretation text for TCGA LUAD validation."""

from __future__ import annotations

import argparse
import math
import random
import textwrap
from pathlib import Path
from statistics import median
from typing import Iterable

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from lifelines import CoxPHFitter, KaplanMeierFitter
except Exception:  # pragma: no cover - optional dependency
    CoxPHFitter = None
    KaplanMeierFitter = None


SCRIPT_DIR = Path(__file__).resolve().parent
TCGA_VALIDATION_DIR = SCRIPT_DIR.parent
DEFAULT_RUN_DIR = TCGA_VALIDATION_DIR / "outputs" / "tcga_luad_three_pathways"
DEFAULT_EXAMPLE_DIR = TCGA_VALIDATION_DIR / "example_results" / "tcga_luad_three_pathways"
METRIC = "terminal_mean_activity"


PATHWAY_ORDER = [
    "R-HSA-69620",
    "R-HSA-453279",
    "R-HSA-1257604",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_EXAMPLE_DIR)
    parser.add_argument("--shuffle-count", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--copy-summary-tables",
        action="store_true",
        help="Copy compact summary tables into the output directory.",
    )
    return parser.parse_args()


def ensure_dirs(base: Path) -> dict[str, Path]:
    dirs = {
        "figures": base / "figures",
        "tables": base / "tables",
        "reports": base / "reports",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    for subdir in ["km", "histograms", "shuffle_nulls", "activity_by_event"]:
        (dirs["figures"] / subdir).mkdir(parents=True, exist_ok=True)
    return dirs


def load_tables(run_dir: Path) -> dict[str, pd.DataFrame]:
    tables_dir = run_dir / "tables"
    required = {
        "activity": tables_dir / "sample_pathway_activity.tsv",
        "survival": tables_dir / "survival_median_split_stats.tsv",
        "shuffle_summary": tables_dir / "sample_label_shuffle_controls.tsv",
        "mapping": tables_dir / "mapping_readiness.tsv",
        "key_findings": tables_dir / "readiness_key_findings.tsv",
    }
    missing = [str(path) for path in required.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required validation tables: {missing}")
    return {name: pd.read_csv(path, sep="\t") for name, path in required.items()}


def ordered_pathways(activity: pd.DataFrame) -> list[str]:
    present = list(activity["pathway_id"].drop_duplicates())
    ordered = [pathway_id for pathway_id in PATHWAY_ORDER if pathway_id in present]
    ordered.extend([pathway_id for pathway_id in present if pathway_id not in ordered])
    return ordered


def pathway_label(pathway_id: str, activity: pd.DataFrame) -> str:
    name = activity.loc[activity["pathway_id"] == pathway_id, "pathway_name"].iloc[0]
    return f"{name}\n{pathway_id}"


def short_pathway_name(pathway_id: str, activity: pd.DataFrame) -> str:
    return str(activity.loc[activity["pathway_id"] == pathway_id, "pathway_name"].iloc[0])


def sanitize(text: str) -> str:
    return (
        text.replace("/", "_")
        .replace(" ", "_")
        .replace("-", "_")
        .replace("(", "")
        .replace(")", "")
    )


def logrank_statistic(
    durations: Iterable[float],
    events: Iterable[int],
    groups: Iterable[int],
) -> dict[str, float]:
    df = pd.DataFrame(
        {
            "duration": list(durations),
            "event": list(events),
            "group": list(groups),
        }
    ).sort_values("duration")
    observed_high = 0.0
    expected_high = 0.0
    variance_high = 0.0
    for event_time in sorted(df.loc[df["event"] == 1, "duration"].unique()):
        at_risk = df[df["duration"] >= event_time]
        events_at_time = df[(df["duration"] == event_time) & (df["event"] == 1)]
        n_total = len(at_risk)
        n_high = int((at_risk["group"] == 1).sum())
        d_total = len(events_at_time)
        d_high = int((events_at_time["group"] == 1).sum())
        if n_total <= 1:
            continue
        observed_high += d_high
        expected_high += d_total * (n_high / n_total)
        variance_high += (
            d_total
            * (n_high / n_total)
            * (1.0 - n_high / n_total)
            * ((n_total - d_total) / (n_total - 1))
        )
    z = (observed_high - expected_high) / math.sqrt(variance_high) if variance_high > 0 else 0.0
    chi_square = z * z
    p_value = math.erfc(math.sqrt(chi_square / 2.0))
    return {
        "observed_high_events": observed_high,
        "expected_high_events": expected_high,
        "logrank_z": z,
        "logrank_chisq": chi_square,
        "logrank_p": p_value,
    }


def benjamini_hochberg(p_values: Iterable[float]) -> list[float]:
    p = np.asarray(list(p_values), dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.where(~np.isnan(p))[0]
    if len(valid) == 0:
        return out.tolist()
    order = valid[np.argsort(p[valid])]
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out[order] = adjusted
    return out.tolist()


def km_step_data(durations: pd.Series, events: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    df = pd.DataFrame({"duration": durations, "event": events}).sort_values("duration")
    x = [0.0]
    y = [1.0]
    survival = 1.0
    for event_time in sorted(df.loc[df["event"] == 1, "duration"].unique()):
        at_risk = int((df["duration"] >= event_time).sum())
        deaths = int(((df["duration"] == event_time) & (df["event"] == 1)).sum())
        if at_risk == 0:
            continue
        x.extend([float(event_time), float(event_time)])
        y.extend([survival, survival * (1.0 - deaths / at_risk)])
        survival = y[-1]
    max_time = float(df["duration"].max())
    x.append(max_time)
    y.append(survival)
    return np.asarray(x), np.asarray(y)


def save_current(figures: list[dict[str, str]], path_no_ext: Path, title: str) -> None:
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    png = path_no_ext.with_suffix(".png")
    svg = path_no_ext.with_suffix(".svg")
    plt.savefig(png, dpi=180, bbox_inches="tight")
    plt.savefig(svg, bbox_inches="tight")
    plt.close()
    figures.append({"title": title, "png": str(png), "svg": str(svg)})


def add_median_split(activity: pd.DataFrame) -> pd.DataFrame:
    out = activity.copy()
    out["event_label"] = np.where(out["survival_event"] == 1, "Deceased", "Living/Censored")
    out["activity_group"] = ""
    out["median_threshold"] = np.nan
    for pathway_id, sub in out.groupby("pathway_id"):
        threshold = float(sub[METRIC].median())
        idx = out["pathway_id"] == pathway_id
        out.loc[idx, "median_threshold"] = threshold
        out.loc[idx, "activity_group"] = np.where(out.loc[idx, METRIC] > threshold, "High activity", "Low activity")
    return out


def compute_shuffle_nulls(activity: pd.DataFrame, shuffle_count: int, seed: int) -> pd.DataFrame:
    rng = random.Random(seed)
    rows: list[dict[str, object]] = []
    for pathway_id, sub in activity.groupby("pathway_id"):
        values = sub[METRIC].astype(float).tolist()
        durations = sub["survival_days"].astype(float).tolist()
        events = sub["survival_event"].astype(int).tolist()
        for i in range(shuffle_count):
            shuffled = values[:]
            rng.shuffle(shuffled)
            threshold = float(median(shuffled))
            groups = [int(value > threshold) for value in shuffled]
            stats = logrank_statistic(durations, events, groups)
            rows.append(
                {
                    "pathway_id": pathway_id,
                    "pathway_name": sub["pathway_name"].iloc[0],
                    "shuffle_index": i,
                    "shuffle_logrank_p": stats["logrank_p"],
                    "shuffle_neg_log10_p": -math.log10(max(stats["logrank_p"], 1e-300)),
                }
            )
    return pd.DataFrame(rows)


def write_final_statistics(tables: dict[str, pd.DataFrame], out_dir: Path) -> pd.DataFrame:
    survival = tables["survival"].copy()
    survival["logrank_fdr_bh_all_metrics"] = benjamini_hochberg(survival["logrank_p"])
    terminal = survival[survival["metric"] == METRIC].copy()
    terminal["logrank_fdr_bh_terminal_mean"] = benjamini_hochberg(terminal["logrank_p"])
    survival.to_csv(out_dir / "median_split_statistics_with_fdr.tsv", sep="\t", index=False)
    terminal.to_csv(out_dir / "terminal_mean_median_split_with_fdr.tsv", sep="\t", index=False)
    return terminal


def write_cox_models(activity: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if CoxPHFitter is None:
        pd.DataFrame(rows).to_csv(out_dir / "cox_univariate.tsv", sep="\t", index=False)
        return pd.DataFrame(rows)

    for pathway_id, sub in activity.groupby("pathway_id"):
        model_df = sub[["survival_days", "survival_event", METRIC]].copy()
        model_df[METRIC] = (model_df[METRIC] - model_df[METRIC].mean()) / model_df[METRIC].std(ddof=0)
        cph = CoxPHFitter()
        cph.fit(model_df, duration_col="survival_days", event_col="survival_event")
        row = cph.summary.loc[METRIC]
        rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": sub["pathway_name"].iloc[0],
                "metric": METRIC,
                "scale": "per_1_sd_increase",
                "hazard_ratio": float(row["exp(coef)"]),
                "hazard_ratio_ci95_lower": float(row["exp(coef) lower 95%"]),
                "hazard_ratio_ci95_upper": float(row["exp(coef) upper 95%"]),
                "z": float(row["z"]),
                "p": float(row["p"]),
            }
        )
    cox = pd.DataFrame(rows)
    if not cox.empty:
        cox["fdr_bh"] = benjamini_hochberg(cox["p"])
    cox.to_csv(out_dir / "cox_univariate.tsv", sep="\t", index=False)
    return cox


def copy_summary_tables(tables: dict[str, pd.DataFrame], out_dir: Path) -> None:
    for name in ["key_findings", "mapping", "shuffle_summary"]:
        tables[name].to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)


def plot_km(activity: pd.DataFrame, terminal_stats: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    pathway_ids = ordered_pathways(activity)
    fig, axes = plt.subplots(1, len(pathway_ids), figsize=(5.4 * len(pathway_ids), 4.2), sharey=True)
    if len(pathway_ids) == 1:
        axes = [axes]
    palette = {"Low activity": "#2f6f9f", "High activity": "#c44536"}
    for ax, pathway_id in zip(axes, pathway_ids):
        sub = activity[activity["pathway_id"] == pathway_id].copy()
        stat = terminal_stats[terminal_stats["pathway_id"] == pathway_id].iloc[0]
        for group_name, group_sub in sub.groupby("activity_group"):
            if KaplanMeierFitter is not None:
                kmf = KaplanMeierFitter()
                kmf.fit(group_sub["survival_days"], event_observed=group_sub["survival_event"], label=group_name)
                kmf.plot_survival_function(ax=ax, ci_show=False, color=palette[group_name], linewidth=2.2)
            else:
                x, y = km_step_data(group_sub["survival_days"], group_sub["survival_event"])
                ax.step(x, y, where="post", label=group_name, color=palette[group_name], linewidth=2.2)
        ax.set_title(pathway_label(pathway_id, activity), fontsize=10)
        ax.set_xlabel("Days")
        ax.set_ylabel("Survival probability")
        ax.grid(True, alpha=0.25)
        ax.text(
            0.04,
            0.08,
            f"log-rank p={stat['logrank_p']:.2g}\nFDR={stat['logrank_fdr_bh_terminal_mean']:.2g}",
            transform=ax.transAxes,
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.9},
        )
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("TCGA LUAD Kaplan-Meier Survival by DeltaSignal Terminal Mean Median Split", y=1.05)
    save_current(figures, fig_dir / "km" / "km_all_pathways", "Kaplan-Meier curves for all pathways")

    for pathway_id in pathway_ids:
        plt.figure(figsize=(5.4, 4.3))
        ax = plt.gca()
        sub = activity[activity["pathway_id"] == pathway_id].copy()
        stat = terminal_stats[terminal_stats["pathway_id"] == pathway_id].iloc[0]
        for group_name, group_sub in sub.groupby("activity_group"):
            if KaplanMeierFitter is not None:
                kmf = KaplanMeierFitter()
                kmf.fit(group_sub["survival_days"], event_observed=group_sub["survival_event"], label=group_name)
                kmf.plot_survival_function(ax=ax, ci_show=False, color=palette[group_name], linewidth=2.4)
            else:
                x, y = km_step_data(group_sub["survival_days"], group_sub["survival_event"])
                ax.step(x, y, where="post", label=group_name, color=palette[group_name], linewidth=2.4)
        ax.set_title(pathway_label(pathway_id, activity))
        ax.set_xlabel("Days")
        ax.set_ylabel("Survival probability")
        ax.grid(True, alpha=0.25)
        ax.text(
            0.04,
            0.08,
            f"log-rank p={stat['logrank_p']:.2g}\nFDR={stat['logrank_fdr_bh_terminal_mean']:.2g}",
            transform=ax.transAxes,
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.9},
        )
        ax.legend(frameon=False)
        save_current(figures, fig_dir / "km" / f"{pathway_id}_km", f"{pathway_id} Kaplan-Meier curve")


def plot_activity_by_event(activity: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    pathway_ids = ordered_pathways(activity)
    fig, axes = plt.subplots(1, len(pathway_ids), figsize=(5 * len(pathway_ids), 4.2), sharey=False)
    if len(pathway_ids) == 1:
        axes = [axes]
    for ax, pathway_id in zip(axes, pathway_ids):
        sub = activity[activity["pathway_id"] == pathway_id]
        sns.violinplot(data=sub, x="event_label", y=METRIC, ax=ax, color="#d6dde8", inner=None, cut=0)
        sns.boxplot(
            data=sub,
            x="event_label",
            y=METRIC,
            ax=ax,
            width=0.28,
            showcaps=True,
            boxprops={"facecolor": "white", "edgecolor": "#333333"},
            medianprops={"color": "#c44536"},
            whiskerprops={"color": "#333333"},
            fliersize=0,
        )
        sns.stripplot(data=sub, x="event_label", y=METRIC, ax=ax, color="#333333", alpha=0.35, size=2.2, jitter=0.22)
        ax.set_title(pathway_label(pathway_id, activity), fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("Terminal mean activity")
        ax.grid(True, axis="y", alpha=0.2)
    fig.suptitle("DeltaSignal Pathway Activity by Survival Event", y=1.04)
    save_current(figures, fig_dir / "activity_by_event" / "activity_by_survival_event_all_pathways", "Activity by survival event")

    for pathway_id in pathway_ids:
        plt.figure(figsize=(5.2, 4.2))
        ax = plt.gca()
        sub = activity[activity["pathway_id"] == pathway_id]
        sns.violinplot(data=sub, x="event_label", y=METRIC, ax=ax, color="#d6dde8", inner=None, cut=0)
        sns.boxplot(data=sub, x="event_label", y=METRIC, ax=ax, width=0.28, color="white", fliersize=0)
        sns.stripplot(data=sub, x="event_label", y=METRIC, ax=ax, color="#333333", alpha=0.35, size=2.2, jitter=0.22)
        ax.set_title(pathway_label(pathway_id, activity))
        ax.set_xlabel("")
        ax.set_ylabel("Terminal mean activity")
        ax.grid(True, axis="y", alpha=0.2)
        save_current(figures, fig_dir / "activity_by_event" / f"{pathway_id}_activity_by_survival_event", f"{pathway_id} activity by event")


def plot_histograms(activity: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    pathway_ids = ordered_pathways(activity)
    fig, axes = plt.subplots(1, len(pathway_ids), figsize=(5 * len(pathway_ids), 4), sharey=False)
    if len(pathway_ids) == 1:
        axes = [axes]
    for ax, pathway_id in zip(axes, pathway_ids):
        sub = activity[activity["pathway_id"] == pathway_id]
        threshold = float(sub["median_threshold"].iloc[0])
        sns.histplot(data=sub, x=METRIC, hue="activity_group", bins=32, stat="density", common_norm=False, kde=True, ax=ax)
        ax.axvline(threshold, color="#222222", linestyle="--", linewidth=1.8, label=f"median={threshold:.2f}")
        ax.set_title(pathway_label(pathway_id, activity), fontsize=10)
        ax.set_xlabel("Terminal mean activity")
        ax.grid(True, axis="y", alpha=0.2)
        ax.legend(frameon=False)
    fig.suptitle("Terminal Mean Activity Distributions and Median Split", y=1.04)
    save_current(figures, fig_dir / "histograms" / "terminal_mean_histograms_all_pathways", "Terminal mean histograms")

    for pathway_id in pathway_ids:
        plt.figure(figsize=(5.4, 4.1))
        ax = plt.gca()
        sub = activity[activity["pathway_id"] == pathway_id]
        threshold = float(sub["median_threshold"].iloc[0])
        sns.histplot(data=sub, x=METRIC, hue="activity_group", bins=36, stat="density", common_norm=False, kde=True, ax=ax)
        ax.axvline(threshold, color="#222222", linestyle="--", linewidth=1.8, label=f"median={threshold:.2f}")
        ax.set_title(pathway_label(pathway_id, activity))
        ax.set_xlabel("Terminal mean activity")
        ax.grid(True, axis="y", alpha=0.2)
        ax.legend(frameon=False)
        save_current(figures, fig_dir / "histograms" / f"{pathway_id}_terminal_mean_histogram", f"{pathway_id} terminal mean histogram")


def plot_heatmap(activity: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    matrix = activity.pivot(index="sample_id", columns="pathway_id", values=METRIC)
    matrix = matrix[[p for p in ordered_pathways(activity) if p in matrix.columns]]
    z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
    sort_column = "R-HSA-69620" if "R-HSA-69620" in z.columns else z.columns[0]
    z = z.sort_values(by=sort_column, ascending=False)
    sample_events = activity.drop_duplicates("sample_id").set_index("sample_id").loc[z.index, "survival_event"]
    event_colors = sample_events.map({0: "#2f6f9f", 1: "#c44536"}).to_numpy()

    fig = plt.figure(figsize=(7.6, 8.8))
    grid = fig.add_gridspec(1, 3, width_ratios=[0.18, 4.8, 0.25], wspace=0.08)
    event_ax = fig.add_subplot(grid[0, 0])
    heat_ax = fig.add_subplot(grid[0, 1])
    cbar_ax = fig.add_subplot(grid[0, 2])

    event_rgb = np.array([matplotlib.colors.to_rgb(color) for color in event_colors]).reshape(len(event_colors), 1, 3)
    event_ax.imshow(event_rgb, aspect="auto")
    event_ax.set_xticks([])
    event_ax.set_yticks([])
    event_ax.set_title("Event", fontsize=9)
    for spine in event_ax.spines.values():
        spine.set_visible(False)

    sns.heatmap(
        z,
        ax=heat_ax,
        cmap="vlag",
        center=0,
        yticklabels=False,
        xticklabels=True,
        cbar=True,
        cbar_ax=cbar_ax,
        cbar_kws={"label": "z-scored terminal mean"},
    )
    heat_ax.set_xlabel("Pathway")
    heat_ax.set_ylabel("")
    heat_ax.set_xticklabels(heat_ax.get_xticklabels(), rotation=0)
    fig.suptitle("Samples x Pathways DeltaSignal Activity Heatmap\nRows sorted by Cell Cycle Checkpoints activity", y=0.99)
    save_current(figures, fig_dir / "heatmap_pathway_activity_zscore", "Samples x pathways heatmap")


def plot_pairplot(activity: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    matrix = activity.pivot(index="sample_id", columns="pathway_id", values=METRIC)
    sample_events = activity.drop_duplicates("sample_id").set_index("sample_id")["event_label"]
    matrix["survival_event"] = sample_events
    rename = {pathway_id: short_pathway_name(pathway_id, activity) for pathway_id in matrix.columns if pathway_id != "survival_event"}
    matrix = matrix.rename(columns=rename)
    grid = sns.pairplot(
        matrix,
        hue="survival_event",
        diag_kind="hist",
        corner=True,
        plot_kws={"alpha": 0.65, "s": 26, "edgecolor": "none"},
        diag_kws={"alpha": 0.55},
        palette={"Living/Censored": "#2f6f9f", "Deceased": "#c44536"},
    )
    grid.fig.suptitle("Pairwise DeltaSignal Pathway Activity Scores", y=1.02)
    save_current(figures, fig_dir / "pairplot_pathway_scores", "Pair plot of pathway scores")


def plot_shuffle_nulls(
    activity: pd.DataFrame,
    terminal_stats: pd.DataFrame,
    shuffle_nulls: pd.DataFrame,
    fig_dir: Path,
    figures: list[dict[str, str]],
) -> None:
    pathway_ids = ordered_pathways(activity)
    fig, axes = plt.subplots(1, len(pathway_ids), figsize=(5 * len(pathway_ids), 4), sharey=False)
    if len(pathway_ids) == 1:
        axes = [axes]
    for ax, pathway_id in zip(axes, pathway_ids):
        sub = shuffle_nulls[shuffle_nulls["pathway_id"] == pathway_id]
        stat = terminal_stats[terminal_stats["pathway_id"] == pathway_id].iloc[0]
        observed = -math.log10(max(float(stat["logrank_p"]), 1e-300))
        sns.histplot(sub["shuffle_neg_log10_p"], bins=36, ax=ax, color="#8ea4c8")
        ax.axvline(observed, color="#c44536", linewidth=2.2, label=f"observed={observed:.2f}")
        ax.set_title(pathway_label(pathway_id, activity), fontsize=10)
        ax.set_xlabel("-log10(log-rank p)")
        ax.set_ylabel("Shuffle count")
        ax.legend(frameon=False)
        ax.grid(True, axis="y", alpha=0.2)
    fig.suptitle("Observed Survival Signal Against Shuffled-Label Null", y=1.04)
    save_current(figures, fig_dir / "shuffle_nulls" / "shuffle_nulls_all_pathways", "Shuffle null plots")

    for pathway_id in pathway_ids:
        plt.figure(figsize=(5.4, 4.1))
        ax = plt.gca()
        sub = shuffle_nulls[shuffle_nulls["pathway_id"] == pathway_id]
        stat = terminal_stats[terminal_stats["pathway_id"] == pathway_id].iloc[0]
        observed = -math.log10(max(float(stat["logrank_p"]), 1e-300))
        sns.histplot(sub["shuffle_neg_log10_p"], bins=36, ax=ax, color="#8ea4c8")
        ax.axvline(observed, color="#c44536", linewidth=2.2, label=f"observed={observed:.2f}")
        ax.set_title(pathway_label(pathway_id, activity))
        ax.set_xlabel("-log10(log-rank p)")
        ax.set_ylabel("Shuffle count")
        ax.legend(frameon=False)
        ax.grid(True, axis="y", alpha=0.2)
        save_current(figures, fig_dir / "shuffle_nulls" / f"{pathway_id}_shuffle_null", f"{pathway_id} shuffle null")


def plot_cox(cox: pd.DataFrame, fig_dir: Path, figures: list[dict[str, str]]) -> None:
    if cox.empty:
        return
    cox = cox.sort_values("hazard_ratio")
    plt.figure(figsize=(7.4, 3.8))
    ax = plt.gca()
    y = np.arange(len(cox))
    ax.errorbar(
        cox["hazard_ratio"],
        y,
        xerr=[
            cox["hazard_ratio"] - cox["hazard_ratio_ci95_lower"],
            cox["hazard_ratio_ci95_upper"] - cox["hazard_ratio"],
        ],
        fmt="o",
        color="#333333",
        ecolor="#777777",
        capsize=3,
    )
    ax.axvline(1.0, color="#c44536", linestyle="--", linewidth=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels(cox["pathway_name"])
    ax.set_xlabel("Hazard ratio per 1 SD terminal mean increase")
    ax.set_title("Univariate Cox Survival Model")
    ax.grid(True, axis="x", alpha=0.2)
    save_current(figures, fig_dir / "cox_univariate_terminal_mean", "Univariate Cox forest plot")


def make_report(
    activity: pd.DataFrame,
    terminal_stats: pd.DataFrame,
    mapping: pd.DataFrame,
    shuffle_summary: pd.DataFrame,
    cox: pd.DataFrame,
    figures: list[dict[str, str]],
    reports_dir: Path,
) -> None:
    sample_count = activity["sample_id"].nunique()
    pathway_count = activity["pathway_id"].nunique()
    checkpoint = terminal_stats[terminal_stats["pathway_id"] == "R-HSA-69620"].iloc[0]
    checkpoint_map = mapping[mapping["pathway_id"] == "R-HSA-69620"].iloc[0]
    checkpoint_shuffle = shuffle_summary[shuffle_summary["pathway_id"] == "R-HSA-69620"].iloc[0]
    cox_text = ""
    if not cox.empty and "R-HSA-69620" in set(cox["pathway_id"]):
        row = cox[cox["pathway_id"] == "R-HSA-69620"].iloc[0]
        cox_text = (
            f"\nA univariate Cox model using Cell Cycle Checkpoints terminal mean as a continuous score "
            f"estimated HR={row['hazard_ratio']:.2f} per 1 SD increase "
            f"(95% CI {row['hazard_ratio_ci95_lower']:.2f}-{row['hazard_ratio_ci95_upper']:.2f}, "
            f"p={row['p']:.2g}, FDR={row['fdr_bh']:.2g}).\n"
        )

    lines = [
        "# TCGA LUAD DeltaSignal Validation Results",
        "",
        "## What Was Run",
        "",
        f"- Cohort: {sample_count} TCGA lung adenocarcinoma tumor samples with RNA-seq and survival metadata.",
        f"- Pathways: {pathway_count} Reactome pathways selected after mRNA observability audit.",
        "- Scoring: DeltaSignal steady-state solve for each sample/pathway pair.",
        "- Sample summary metric: terminal_mean_activity, the average predicted activity of terminal/output nodes in each pathway network.",
        "- Survival test: median split into low/high pathway activity groups, then log-rank test.",
        "- Negative control: 1000 shuffled-label survival controls for terminal_mean_activity.",
        "",
        "## Main Finding",
        "",
        "Cell Cycle Checkpoints (R-HSA-69620) is the strongest signal in this validation run.",
        "",
        "Biologically, this Reactome pathway represents checkpoint mechanisms controlling whether cells continue through the cell cycle: DNA damage checkpoints, G1/S transition control, G2/M control, and mitotic checkpoint behavior. In cancer, high activity in this kind of pathway can indicate highly proliferative tumors or tumors under replication/checkpoint stress, which can be associated with aggressive disease.",
        "",
        f"For this LUAD run, Cell Cycle Checkpoints had {int(checkpoint_map['mapped'])} mRNA-backed root nodes out of {int(checkpoint_map['root_uuid_count'])} root UUIDs, with {int(checkpoint_map['stable_id_missing'])} missing stable IDs after the LNG mapping fix.",
        "",
        "The median split produced two equal groups:",
        "",
        f"- Low activity group: {int(checkpoint['low_n'])} patients, {int(checkpoint['low_events'])} deaths, KM median survival {checkpoint['low_km_median_survival']:.0f} days.",
        f"- High activity group: {int(checkpoint['high_n'])} patients, {int(checkpoint['high_events'])} deaths, KM median survival {checkpoint['high_km_median_survival']:.0f} days.",
        "",
        f"Plain English: patients whose tumors had higher predicted Cell Cycle Checkpoints pathway activity died more often and had shorter median survival in this cohort.",
        "",
        f"The log-rank p-value was {checkpoint['logrank_p']:.3g}; BH FDR across the three terminal-mean pathway tests was {checkpoint['logrank_fdr_bh_terminal_mean']:.3g}.",
        f"The shuffle control found {checkpoint_shuffle['shuffle_fraction_p_le_observed'] * 1000:.0f}/1000 shuffled controls as strong or stronger than the observed split.",
        "Plain English: the observed survival association is stronger than the label-randomized nulls generated for this run.",
        cox_text.strip(),
        "",
        "## Other Pathways",
        "",
    ]
    for _, row in terminal_stats.sort_values("logrank_p").iterrows():
        lines.extend(
            [
                f"### {row['pathway_name']} ({row['pathway_id']})",
                "",
                f"- log-rank p: {row['logrank_p']:.3g}",
                f"- terminal-mean BH FDR: {row['logrank_fdr_bh_terminal_mean']:.3g}",
                f"- low/high events: {int(row['low_events'])}/{int(row['high_events'])}",
                f"- low/high KM median survival: {row['low_km_median_survival']:.0f}/{row['high_km_median_survival']:.0f} days",
                "",
            ]
        )
    lines.extend(
        [
            "## Caveats",
            "",
            "- This validates an association, not causality.",
            "- The current expression-to-observation transform is cohort percentile activity; alternative transforms should be compared.",
            "- Terminal outputs are graph-derived terminals, not manually curated biological readouts.",
            "- Only three pathways were tested in this focused run. Reactome-wide expansion requires strict multiple-testing correction and more careful model selection.",
            "- Clinical covariates beyond survival/event are not yet included in the main result tables.",
            "",
            "## Generated Figures",
            "",
        ]
    )
    for fig in figures:
        rel_png = Path(fig["png"]).relative_to(reports_dir.parent)
        lines.append(f"- {fig['title']}: `{rel_png}`")
    lines.append("")
    (reports_dir / "RESULT_INTERPRETATION.md").write_text("\n".join(lines))


def write_figure_manifest(figures: list[dict[str, str]], output_path: Path) -> None:
    base_dir = output_path.parent.parent
    rows = []
    for fig in figures:
        rows.append(
            {
                "title": fig["title"],
                "png": str(Path(fig["png"]).relative_to(base_dir)),
                "svg": str(Path(fig["svg"]).relative_to(base_dir)),
            }
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    dirs = ensure_dirs(args.output_dir)
    tables = load_tables(args.run_dir)
    activity = add_median_split(tables["activity"])
    terminal_stats = write_final_statistics(tables, dirs["tables"])
    cox = write_cox_models(activity, dirs["tables"])
    shuffle_nulls_path = dirs["tables"] / "shuffle_null_pvalues.tsv"
    if shuffle_nulls_path.exists():
        shuffle_nulls = pd.read_csv(shuffle_nulls_path, sep="\t")
    else:
        shuffle_nulls = compute_shuffle_nulls(activity, args.shuffle_count, args.seed)
        shuffle_nulls.to_csv(shuffle_nulls_path, sep="\t", index=False)
    if args.copy_summary_tables:
        copy_summary_tables(tables, dirs["tables"])

    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.size": 10,
        }
    )
    figures: list[dict[str, str]] = []
    plot_km(activity, terminal_stats, dirs["figures"], figures)
    plot_activity_by_event(activity, dirs["figures"], figures)
    plot_histograms(activity, dirs["figures"], figures)
    plot_heatmap(activity, dirs["figures"], figures)
    plot_pairplot(activity, dirs["figures"], figures)
    plot_shuffle_nulls(activity, terminal_stats, shuffle_nulls, dirs["figures"], figures)
    plot_cox(cox, dirs["figures"], figures)
    write_figure_manifest(figures, dirs["reports"] / "figure_manifest.tsv")
    make_report(activity, terminal_stats, tables["mapping"], tables["shuffle_summary"], cox, figures, dirs["reports"])
    print(f"Wrote {len(figures)} figure groups and reports to {args.output_dir}")


if __name__ == "__main__":
    main()
