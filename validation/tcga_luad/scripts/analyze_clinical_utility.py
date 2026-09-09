#!/usr/bin/env python3
"""Evaluate whether historical TCGA LUAD DeltaSignal scores add clinical signal.

This analysis is deliberately separate from the perturbation benchmark. It
uses the completed 502-tumor historical run and asks whether a continuous
DeltaSignal pathway score remains associated with overall survival after
clinical adjustment and after comparison with simple expression-only scores.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import chi2, pearsonr, spearmanr


PATHWAY_ORDER = ["R-HSA-1257604", "R-HSA-453279", "R-HSA-69620"]
PATHWAY_SHORT_NAMES = {
    "R-HSA-1257604": "PIP3/AKT",
    "R-HSA-453279": "Mitotic G1/G1-S",
    "R-HSA-69620": "Cell Cycle Checkpoints",
}
ACTIVITY_METRICS = [
    "terminal_mean_activity",
    "terminal_median_activity",
    "terminal_q25_activity",
    "terminal_q75_activity",
]
CLINICAL_FEATURES = ["age_z", "male", "ever_smoker", "stage_II", "stage_III", "stage_IV"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activity-tsv", type=Path, required=True)
    parser.add_argument("--clinical-tsv", type=Path, required=True)
    parser.add_argument("--expression-tsv", type=Path, required=True)
    parser.add_argument("--root-audit-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--cv-repeats", type=int, default=20)
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_revision(path: Path) -> dict[str, object]:
    root = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=False,
    )
    if root.returncode != 0:
        return {"available": False}
    repo = Path(root.stdout.strip())

    def git(*parts: str) -> str:
        result = subprocess.run(
            ["git", "-C", str(repo), *parts],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    return {
        "available": True,
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "dirty": bool(git("status", "--porcelain")),
    }


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    adjusted = np.empty(len(p), dtype=float)
    running = 1.0
    for rank_index in range(len(p) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p[original_index] * len(p) / rank)
        adjusted[original_index] = min(running, 1.0)
    return pd.Series(adjusted, index=values.index)


def stage_group(value: object) -> str | None:
    if pd.isna(value):
        return None
    match = re.match(r"Stage\s+(IV|III|II|I)", str(value).strip(), flags=re.IGNORECASE)
    return match.group(1).upper() if match else None


def load_clinical(path: Path, sample_ids: set[str]) -> pd.DataFrame:
    columns = [
        "sampleID",
        "age_at_initial_pathologic_diagnosis",
        "gender",
        "pathologic_stage",
        "tobacco_smoking_history",
    ]
    clinical = pd.read_csv(path, sep="\t", usecols=columns, low_memory=False)
    clinical = clinical[clinical["sampleID"].isin(sample_ids)].drop_duplicates("sampleID")
    clinical = clinical.rename(columns={"sampleID": "sample_id"}).set_index("sample_id")
    clinical["age"] = pd.to_numeric(
        clinical["age_at_initial_pathologic_diagnosis"], errors="coerce"
    )
    clinical["male"] = clinical["gender"].map({"FEMALE": 0.0, "MALE": 1.0})
    smoking = pd.to_numeric(clinical["tobacco_smoking_history"], errors="coerce")
    clinical["ever_smoker"] = np.where(smoking.notna(), (smoking != 1).astype(float), np.nan)
    clinical["stage_group"] = clinical["pathologic_stage"].map(stage_group)
    return clinical


def load_expression_baselines(
    expression_path: Path,
    root_audit_path: Path,
    sample_ids: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    audit = pd.read_csv(root_audit_path, sep="\t")
    audit = audit[
        audit["pathway_id"].isin(PATHWAY_ORDER) & (audit["mapping_status"] == "mapped")
    ].copy()
    audit["gene_symbol"] = audit["gene_symbol"].astype("string").str.strip()
    audit = audit[audit["gene_symbol"].notna() & (audit["gene_symbol"] != "")]
    genes_by_pathway = {
        pathway_id: sorted(set(sub["gene_symbol"]))
        for pathway_id, sub in audit.groupby("pathway_id")
    }
    requested_genes = sorted(set().union(*genes_by_pathway.values(), {"MKI67"}))

    expression = pd.read_csv(expression_path, sep="\t")
    expression = expression.rename(columns={expression.columns[0]: "gene"})
    missing_samples = sorted(set(sample_ids) - set(expression.columns))
    if missing_samples:
        raise ValueError(f"Expression matrix is missing {len(missing_samples)} scored samples")
    expression = expression[expression["gene"].isin(requested_genes)][["gene", *sample_ids]]
    expression = expression.groupby("gene", as_index=True).mean(numeric_only=True)
    percentiles = expression.rank(axis=1, method="average", pct=True) * 100.0

    rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    for pathway_id in PATHWAY_ORDER:
        requested = genes_by_pathway.get(pathway_id, [])
        found = sorted(set(requested).intersection(percentiles.index))
        if not found:
            raise ValueError(f"No mapped genes found in expression for {pathway_id}")
        mean_percentile = percentiles.loc[found].mean(axis=0)
        for sample_id, value in mean_percentile.items():
            rows.append(
                {
                    "sample_id": sample_id,
                    "pathway_id": pathway_id,
                    "expression_baseline": float(value),
                    "mki67_percentile": (
                        float(percentiles.at["MKI67", sample_id])
                        if "MKI67" in percentiles.index
                        else np.nan
                    ),
                }
            )
        coverage_rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": PATHWAY_SHORT_NAMES[pathway_id],
                "unique_mapped_genes_requested": len(requested),
                "unique_mapped_genes_in_expression": len(found),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(coverage_rows)


def prepare_pathway_frame(
    activity: pd.DataFrame,
    baselines: pd.DataFrame,
    clinical: pd.DataFrame,
    pathway_id: str,
) -> pd.DataFrame:
    sub = activity[activity["pathway_id"] == pathway_id].copy()
    if sub["sample_id"].duplicated().any():
        raise ValueError(f"Duplicate activity rows for {pathway_id}")
    sub = sub.merge(
        baselines[baselines["pathway_id"] == pathway_id],
        on=["sample_id", "pathway_id"],
        how="left",
        validate="one_to_one",
    )
    sub = sub.join(clinical, on="sample_id", validate="many_to_one")
    stage_dummies = pd.get_dummies(sub["stage_group"], prefix="stage", dtype=float)
    for column in ["stage_II", "stage_III", "stage_IV"]:
        sub[column] = stage_dummies[column] if column in stage_dummies else 0.0
    required = [
        "survival_days",
        "survival_event",
        "age",
        "male",
        "ever_smoker",
        "stage_group",
        "expression_baseline",
        "mki67_percentile",
        *ACTIVITY_METRICS,
    ]
    sub = sub.dropna(subset=required)
    sub = sub[sub["stage_group"].isin(["I", "II", "III", "IV"])].copy()
    sub["age_z"] = standardize(sub["age"])
    sub["ds_z"] = standardize(sub["terminal_mean_activity"])
    sub["baseline_z"] = standardize(sub["expression_baseline"])
    sub["mki67_z"] = standardize(sub["mki67_percentile"])
    return sub


def standardize(values: pd.Series) -> pd.Series:
    std = float(values.std(ddof=0))
    if not math.isfinite(std) or std == 0:
        raise ValueError(f"Cannot standardize constant values: {values.name}")
    return (values - values.mean()) / std


def fit_cox(frame: pd.DataFrame, features: list[str]) -> CoxPHFitter:
    model = frame[["survival_days", "survival_event", *features]].copy()
    fitter = CoxPHFitter()
    fitter.fit(model, duration_col="survival_days", event_col="survival_event")
    return fitter


def coefficient_row(
    pathway_id: str,
    model_name: str,
    fitter: CoxPHFitter,
    feature: str,
) -> dict[str, object]:
    row = fitter.summary.loc[feature]
    return {
        "pathway_id": pathway_id,
        "pathway_name": PATHWAY_SHORT_NAMES[pathway_id],
        "model": model_name,
        "feature": feature,
        "n": int(fitter._n_examples),
        "events": int(fitter.event_observed.sum()),
        "hazard_ratio": float(row["exp(coef)"]),
        "hazard_ratio_ci95_lower": float(row["exp(coef) lower 95%"]),
        "hazard_ratio_ci95_upper": float(row["exp(coef) upper 95%"]),
        "z": float(row["z"]),
        "p": float(row["p"]),
    }


def likelihood_ratio(full: CoxPHFitter, reduced: CoxPHFitter, df: int = 1) -> float:
    statistic = max(0.0, 2.0 * (full.log_likelihood_ - reduced.log_likelihood_))
    return float(chi2.sf(statistic, df=df))


def fit_model_suite(
    frame: pd.DataFrame, pathway_id: str
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    specs = {
        "clinical": CLINICAL_FEATURES,
        "clinical_plus_deltasignal": [*CLINICAL_FEATURES, "ds_z"],
        "clinical_plus_expression": [*CLINICAL_FEATURES, "baseline_z"],
        "clinical_plus_expression_plus_deltasignal": [
            *CLINICAL_FEATURES,
            "baseline_z",
            "ds_z",
        ],
        "clinical_plus_mki67": [*CLINICAL_FEATURES, "mki67_z"],
        "clinical_plus_mki67_plus_deltasignal": [
            *CLINICAL_FEATURES,
            "mki67_z",
            "ds_z",
        ],
    }
    fits = {name: fit_cox(frame, features) for name, features in specs.items()}
    comparisons = {
        "clinical_plus_deltasignal": ("clinical", "ds_added_to_clinical"),
        "clinical_plus_expression": ("clinical", "expression_added_to_clinical"),
        "clinical_plus_expression_plus_deltasignal": (
            "clinical_plus_expression",
            "deltasignal_added_beyond_expression",
        ),
        "clinical_plus_mki67": ("clinical", "mki67_added_to_clinical"),
        "clinical_plus_mki67_plus_deltasignal": (
            "clinical_plus_mki67",
            "deltasignal_added_beyond_mki67",
        ),
    }
    model_rows: list[dict[str, object]] = []
    for name, fitter in fits.items():
        reduced_name, test_name = comparisons.get(name, (None, None))
        model_rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": PATHWAY_SHORT_NAMES[pathway_id],
                "model": name,
                "n": len(frame),
                "events": int(frame["survival_event"].sum()),
                "concordance_index_in_sample": float(fitter.concordance_index_),
                "partial_aic": float(fitter.AIC_partial_),
                "log_likelihood": float(fitter.log_likelihood_),
                "nested_test": test_name,
                "nested_lrt_p": (
                    likelihood_ratio(fitter, fits[reduced_name]) if reduced_name else np.nan
                ),
            }
        )

    coefficient_rows = [
        coefficient_row(pathway_id, "unadjusted_deltasignal", fit_cox(frame, ["ds_z"]), "ds_z"),
        coefficient_row(
            pathway_id,
            "clinical_plus_deltasignal",
            fits["clinical_plus_deltasignal"],
            "ds_z",
        ),
        coefficient_row(
            pathway_id,
            "clinical_plus_expression_plus_deltasignal",
            fits["clinical_plus_expression_plus_deltasignal"],
            "ds_z",
        ),
        coefficient_row(
            pathway_id,
            "clinical_plus_expression_plus_deltasignal",
            fits["clinical_plus_expression_plus_deltasignal"],
            "baseline_z",
        ),
        coefficient_row(
            pathway_id,
            "clinical_plus_mki67_plus_deltasignal",
            fits["clinical_plus_mki67_plus_deltasignal"],
            "ds_z",
        ),
        coefficient_row(
            pathway_id,
            "clinical_plus_mki67_plus_deltasignal",
            fits["clinical_plus_mki67_plus_deltasignal"],
            "mki67_z",
        ),
    ]
    return model_rows, coefficient_rows


def stratified_folds(events: pd.Series, folds: int, rng: np.random.Generator) -> np.ndarray:
    assignments = np.empty(len(events), dtype=int)
    event_array = events.to_numpy(dtype=int)
    for event_value in [0, 1]:
        indices = np.flatnonzero(event_array == event_value)
        rng.shuffle(indices)
        assignments[indices] = np.arange(len(indices)) % folds
    return assignments


def cross_validated_models(
    frame: pd.DataFrame,
    pathway_id: str,
    folds: int,
    repeats: int,
    seed: int,
) -> list[dict[str, object]]:
    specs = {
        "clinical": CLINICAL_FEATURES,
        "clinical_plus_deltasignal": [*CLINICAL_FEATURES, "ds_z"],
        "clinical_plus_expression": [*CLINICAL_FEATURES, "baseline_z"],
        "clinical_plus_expression_plus_deltasignal": [
            *CLINICAL_FEATURES,
            "baseline_z",
            "ds_z",
        ],
        "clinical_plus_mki67": [*CLINICAL_FEATURES, "mki67_z"],
        "clinical_plus_mki67_plus_deltasignal": [
            *CLINICAL_FEATURES,
            "mki67_z",
            "ds_z",
        ],
    }
    rows: list[dict[str, object]] = []
    for repeat in range(repeats):
        rng = np.random.default_rng(seed + repeat)
        assignments = stratified_folds(frame["survival_event"], folds, rng)
        for model_name, features in specs.items():
            fold_cindices: list[float] = []
            for fold in range(folds):
                train = frame.iloc[assignments != fold]
                test = frame.iloc[assignments == fold]
                fitter = fit_cox(train, features)
                risks = fitter.predict_partial_hazard(test[features]).to_numpy()
                fold_cindices.append(
                    float(
                        concordance_index(
                            test["survival_days"], -risks, test["survival_event"]
                        )
                    )
                )
            rows.append(
                {
                    "pathway_id": pathway_id,
                    "pathway_name": PATHWAY_SHORT_NAMES[pathway_id],
                    "repeat": repeat,
                    "folds": folds,
                    "model": model_name,
                    "concordance_index": float(np.mean(fold_cindices)),
                    "fold_cindex_min": min(fold_cindices),
                    "fold_cindex_max": max(fold_cindices),
                }
            )
    return rows


def sensitivity_models(frame: pd.DataFrame, pathway_id: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for metric in ACTIVITY_METRICS:
        feature = f"{metric}_z"
        frame = frame.copy()
        frame[feature] = standardize(frame[metric])
        fitter = fit_cox(frame, [*CLINICAL_FEATURES, feature])
        row = coefficient_row(pathway_id, "clinical_adjusted", fitter, feature)
        row["metric"] = metric
        rows.append(row)
    return rows


def score_correlations(frame: pd.DataFrame, pathway_id: str) -> dict[str, object]:
    stage_numeric = frame["stage_group"].map({"I": 1, "II": 2, "III": 3, "IV": 4})

    def rho(left: str, right: pd.Series | str) -> tuple[float, float]:
        other = frame[right] if isinstance(right, str) else right
        result = spearmanr(frame[left], other)
        return float(result.statistic), float(result.pvalue)

    ds_expression_rho, ds_expression_p = rho("terminal_mean_activity", "expression_baseline")
    ds_mki67_rho, ds_mki67_p = rho("terminal_mean_activity", "mki67_percentile")
    ds_stage_rho, ds_stage_p = rho("terminal_mean_activity", stage_numeric)
    expression_pearson = pearsonr(
        frame["terminal_mean_activity"], frame["expression_baseline"]
    )
    expression_vif = 1.0 / (1.0 - float(expression_pearson.statistic) ** 2)
    return {
        "pathway_id": pathway_id,
        "pathway_name": PATHWAY_SHORT_NAMES[pathway_id],
        "n": len(frame),
        "deltasignal_vs_expression_spearman": ds_expression_rho,
        "deltasignal_vs_expression_p": ds_expression_p,
        "deltasignal_vs_expression_pearson": float(expression_pearson.statistic),
        "two_predictor_vif": expression_vif,
        "deltasignal_vs_mki67_spearman": ds_mki67_rho,
        "deltasignal_vs_mki67_p": ds_mki67_p,
        "deltasignal_vs_stage_spearman": ds_stage_rho,
        "deltasignal_vs_stage_p": ds_stage_p,
    }


def plot_adjusted_forest(coefficients: pd.DataFrame, output_dir: Path) -> None:
    data = coefficients[
        (coefficients["model"] == "clinical_plus_expression_plus_deltasignal")
        & (coefficients["feature"] == "ds_z")
    ].copy()
    data["pathway_name"] = pd.Categorical(
        data["pathway_name"],
        [PATHWAY_SHORT_NAMES[pathway] for pathway in PATHWAY_ORDER],
        ordered=True,
    )
    data = data.sort_values("pathway_name")
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    y = np.arange(len(data))
    ax.errorbar(
        data["hazard_ratio"],
        y,
        xerr=np.vstack(
            [
                data["hazard_ratio"] - data["hazard_ratio_ci95_lower"],
                data["hazard_ratio_ci95_upper"] - data["hazard_ratio"],
            ]
        ),
        fmt="o",
        color="#176b87",
        ecolor="#5d7f8c",
        capsize=4,
    )
    ax.axvline(1.0, color="#444444", linewidth=1, linestyle="--")
    ax.set_yticks(y, data["pathway_name"])
    ax.set_xlabel("Hazard ratio per 1 SD DeltaSignal increase")
    ax.set_title("Adjusted for clinical covariates and expression baseline")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "adjusted_deltasignal_forest.png", dpi=180)
    fig.savefig(output_dir / "adjusted_deltasignal_forest.svg")
    plt.close(fig)


def plot_cv(cv: pd.DataFrame, output_dir: Path) -> None:
    display = {
        "clinical": "Clinical",
        "clinical_plus_deltasignal": "Clinical + DS",
        "clinical_plus_expression": "Clinical + expression",
        "clinical_plus_expression_plus_deltasignal": "Clinical + expression + DS",
    }
    data = cv[cv["model"].isin(display)].copy()
    data["model_label"] = data["model"].map(display)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.7), sharey=True)
    for ax, pathway_id in zip(axes, PATHWAY_ORDER):
        sub = data[data["pathway_id"] == pathway_id]
        sns.boxplot(
            data=sub,
            x="model_label",
            y="concordance_index",
            color="#8fb8c8",
            width=0.6,
            ax=ax,
        )
        ax.set_title(PATHWAY_SHORT_NAMES[pathway_id])
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Repeated 5-fold CV concordance")
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_dir / "cross_validated_concordance.png", dpi=180)
    fig.savefig(output_dir / "cross_validated_concordance.svg")
    plt.close(fig)


def plot_score_comparison(frames: dict[str, pd.DataFrame], output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
    for ax, pathway_id in zip(axes, PATHWAY_ORDER):
        frame = frames[pathway_id]
        sns.scatterplot(
            data=frame,
            x="expression_baseline",
            y="terminal_mean_activity",
            hue="survival_event",
            palette={0: "#4f7cac", 1: "#c75146"},
            alpha=0.65,
            s=30,
            linewidth=0,
            ax=ax,
            legend=pathway_id == PATHWAY_ORDER[-1],
        )
        ax.set_title(PATHWAY_SHORT_NAMES[pathway_id])
        ax.set_xlabel("Mean mapped-gene percentile")
        ax.set_ylabel("DeltaSignal terminal mean")
    if axes[-1].get_legend() is not None:
        axes[-1].legend(title="Death event", loc="best")
    fig.tight_layout()
    fig.savefig(output_dir / "deltasignal_vs_expression_baseline.png", dpi=180)
    fig.savefig(output_dir / "deltasignal_vs_expression_baseline.svg")
    plt.close(fig)


def summarize_cv(cv: pd.DataFrame) -> pd.DataFrame:
    summary = (
        cv.groupby(["pathway_id", "pathway_name", "model"])["concordance_index"]
        .agg(
            cv_cindex_median="median",
            cv_cindex_min="min",
            cv_cindex_max="max",
        )
        .reset_index()
    )
    clinical = summary[summary["model"] == "clinical"][
        ["pathway_id", "cv_cindex_median"]
    ].rename(columns={"cv_cindex_median": "clinical_cv_cindex_median"})
    summary = summary.merge(clinical, on="pathway_id", how="left")
    summary["median_delta_vs_clinical"] = (
        summary["cv_cindex_median"] - summary["clinical_cv_cindex_median"]
    )
    expression = summary[summary["model"] == "clinical_plus_expression"][
        ["pathway_id", "cv_cindex_median"]
    ].rename(columns={"cv_cindex_median": "expression_cv_cindex_median"})
    summary = summary.merge(expression, on="pathway_id", how="left")
    summary["median_delta_vs_expression"] = (
        summary["cv_cindex_median"] - summary["expression_cv_cindex_median"]
    )
    return summary


def write_report(
    output_dir: Path,
    coefficients: pd.DataFrame,
    models: pd.DataFrame,
    cv_summary: pd.DataFrame,
    correlations: pd.DataFrame,
    sensitivity: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    lines = [
        "# TCGA LUAD clinical-utility analysis",
        "",
        "## Scope",
        "",
        "This is a secondary external-association analysis of the completed historical "
        "502-tumor run. It is not perturbation ground truth and it does not represent "
        "the later current-stack LNG networks, whose pathway-level aggregation remains "
        "unresolved.",
        "",
        "DeltaSignal terminal-mean scores were tested continuously. Clinical models "
        "adjust for age, sex, pathologic stage, and ever-smoking history. The expression "
        "baseline is the mean cohort percentile across unique mapped genes used as "
        "DeltaSignal inputs, so position-aware UUID duplication does not reweight genes.",
        "",
        "## Results",
        "",
    ]
    for pathway_id in PATHWAY_ORDER:
        coef = coefficients[
            (coefficients["pathway_id"] == pathway_id)
            & (coefficients["model"] == "clinical_plus_expression_plus_deltasignal")
            & (coefficients["feature"] == "ds_z")
        ].iloc[0]
        model = models[
            (models["pathway_id"] == pathway_id)
            & (models["model"] == "clinical_plus_expression_plus_deltasignal")
        ].iloc[0]
        cv_clinical = cv_summary[
            (cv_summary["pathway_id"] == pathway_id)
            & (cv_summary["model"] == "clinical")
        ].iloc[0]
        cv_ds = cv_summary[
            (cv_summary["pathway_id"] == pathway_id)
            & (cv_summary["model"] == "clinical_plus_deltasignal")
        ].iloc[0]
        cv_expression = cv_summary[
            (cv_summary["pathway_id"] == pathway_id)
            & (cv_summary["model"] == "clinical_plus_expression")
        ].iloc[0]
        cv_full = cv_summary[
            (cv_summary["pathway_id"] == pathway_id)
            & (cv_summary["model"] == "clinical_plus_expression_plus_deltasignal")
        ].iloc[0]
        corr = correlations[correlations["pathway_id"] == pathway_id].iloc[0]
        lines.extend(
            [
                f"### {PATHWAY_SHORT_NAMES[pathway_id]}",
                "",
                f"- Complete-case cohort: {int(coef['n'])} tumors, {int(coef['events'])} deaths.",
                f"- DeltaSignal after clinical and expression adjustment: HR {coef['hazard_ratio']:.3f} "
                f"(95% CI {coef['hazard_ratio_ci95_lower']:.3f}-{coef['hazard_ratio_ci95_upper']:.3f}), "
                f"p={coef['p']:.3g}, three-pathway FDR={coef['fdr_bh_within_model_feature']:.3g} "
                "per 1 SD increase.",
                f"- Nested test for adding DeltaSignal beyond clinical plus expression: "
                f"p={model['nested_lrt_p']:.3g}.",
                f"- Repeated 5-fold CV median concordance: clinical "
                f"{cv_clinical['cv_cindex_median']:.3f}; plus DeltaSignal "
                f"{cv_ds['cv_cindex_median']:.3f}; plus expression "
                f"{cv_expression['cv_cindex_median']:.3f}; plus both "
                f"{cv_full['cv_cindex_median']:.3f}.",
                f"- Adding DeltaSignal to the expression model changed median held-out "
                f"concordance by {cv_full['median_delta_vs_expression']:+.3f}.",
                f"- DeltaSignal/expression-baseline Spearman rho: "
                f"{corr['deltasignal_vs_expression_spearman']:.3f}; two-predictor VIF: "
                f"{corr['two_predictor_vif']:.1f}.",
                "",
            ]
        )

    checkpoint = sensitivity[sensitivity["pathway_id"] == "R-HSA-69620"].copy()
    lines.extend(
        [
            "## Sensitivity and interpretation",
            "",
            "Cell Cycle Checkpoints is the prespecified TCGA signal of interest. Its "
            "clinical-adjusted estimates across terminal mean, median, lower quartile, "
            "and upper quartile are retained in `tables/aggregation_sensitivity.tsv`; "
            "disagreement between these summaries is evidence that the pathway-level "
            "readout is aggregation-sensitive.",
            "",
            "Cell Cycle Checkpoints retains a nominal coefficient after clinical and "
            "expression adjustment, but that coefficient does not survive the "
            "three-pathway FDR correction. The DeltaSignal and expression scores are "
            "strongly collinear, and adding both does not improve repeated-CV "
            "concordance over expression alone. The defensible conclusion is that the "
            "historical DeltaSignal score captures clinically relevant cell-cycle "
            "variation, but this analysis does not establish incremental predictive "
            "value from network propagation.",
            "",
            f"The analysis used {int(coverage['unique_mapped_genes_in_expression'].sum())} "
            "pathway-specific unique-gene entries across the three expression baselines. "
            "These counts overlap between pathways and must not be interpreted as that "
            "many distinct genes overall.",
            "",
            "A survival association can demonstrate clinical relevance, but not causal "
            "perturbation accuracy. Any in-sample likelihood-ratio result is explanatory; "
            "the repeated cross-validation is the more conservative discrimination check. "
            "No result here repairs the current-stack terminal aggregation problem.",
            "",
            "## Files",
            "",
            "- `tables/cox_coefficients.tsv`: continuous hazard-ratio estimates.",
            "- `tables/model_comparison.tsv`: nested tests and in-sample fit.",
            "- `tables/repeated_cv_concordance.tsv`: all repeated-CV results.",
            "- `tables/repeated_cv_summary.tsv`: compact discrimination summary.",
            "- `tables/score_correlations.tsv`: score/baseline/stage correlations.",
            "- `tables/aggregation_sensitivity.tsv`: alternate output summaries.",
            "- `analysis_manifest.json`: revisions, parameters, and input hashes.",
            "",
        ]
    )
    if checkpoint.empty:
        raise ValueError("Cell Cycle Checkpoints sensitivity results are missing")
    (output_dir / "TCGA_CLINICAL_UTILITY_REPORT.md").write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    if args.cv_folds < 2 or args.cv_repeats < 1:
        raise ValueError("Cross-validation requires at least 2 folds and 1 repeat")
    tables_dir = args.output_dir / "tables"
    figures_dir = args.output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    activity = pd.read_csv(args.activity_tsv, sep="\t")
    missing_columns = {
        "sample_id",
        "pathway_id",
        "survival_days",
        "survival_event",
        *ACTIVITY_METRICS,
    } - set(activity.columns)
    if missing_columns:
        raise ValueError(f"Activity table is missing columns: {sorted(missing_columns)}")
    activity = activity[activity["pathway_id"].isin(PATHWAY_ORDER)].copy()
    sample_ids = sorted(activity["sample_id"].unique())
    if len(sample_ids) != 502:
        raise ValueError(f"Expected the locked 502-sample cohort, found {len(sample_ids)}")

    clinical = load_clinical(args.clinical_tsv, set(sample_ids))
    baselines, coverage = load_expression_baselines(
        args.expression_tsv, args.root_audit_tsv, sample_ids
    )

    frames: dict[str, pd.DataFrame] = {}
    model_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    cv_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []
    correlation_rows: list[dict[str, object]] = []
    for pathway_id in PATHWAY_ORDER:
        frame = prepare_pathway_frame(activity, baselines, clinical, pathway_id)
        frames[pathway_id] = frame
        pathway_models, pathway_coefficients = fit_model_suite(frame, pathway_id)
        model_rows.extend(pathway_models)
        coefficient_rows.extend(pathway_coefficients)
        cv_rows.extend(
            cross_validated_models(
                frame, pathway_id, args.cv_folds, args.cv_repeats, args.seed
            )
        )
        sensitivity_rows.extend(sensitivity_models(frame, pathway_id))
        correlation_rows.append(score_correlations(frame, pathway_id))

    models = pd.DataFrame(model_rows)
    coefficients = pd.DataFrame(coefficient_rows)
    cv = pd.DataFrame(cv_rows)
    sensitivity = pd.DataFrame(sensitivity_rows)
    correlations = pd.DataFrame(correlation_rows)
    coefficients["fdr_bh_within_model_feature"] = coefficients.groupby(
        ["model", "feature"]
    )["p"].transform(benjamini_hochberg)
    sensitivity["fdr_bh_within_metric"] = sensitivity.groupby("metric")["p"].transform(
        benjamini_hochberg
    )
    cv_summary = summarize_cv(cv)

    models.to_csv(tables_dir / "model_comparison.tsv", sep="\t", index=False)
    coefficients.to_csv(tables_dir / "cox_coefficients.tsv", sep="\t", index=False)
    cv.to_csv(tables_dir / "repeated_cv_concordance.tsv", sep="\t", index=False)
    cv_summary.to_csv(tables_dir / "repeated_cv_summary.tsv", sep="\t", index=False)
    sensitivity.to_csv(tables_dir / "aggregation_sensitivity.tsv", sep="\t", index=False)
    correlations.to_csv(tables_dir / "score_correlations.tsv", sep="\t", index=False)
    coverage.to_csv(tables_dir / "expression_baseline_coverage.tsv", sep="\t", index=False)
    pd.concat(frames.values(), ignore_index=True).to_csv(
        tables_dir / "analysis_cohort.tsv", sep="\t", index=False
    )

    plot_adjusted_forest(coefficients, figures_dir)
    plot_cv(cv, figures_dir)
    plot_score_comparison(frames, figures_dir)
    write_report(
        args.output_dir,
        coefficients,
        models,
        cv_summary,
        correlations,
        sensitivity,
        coverage,
    )

    manifest = {
        "analysis": "historical_tcga_luad_clinical_utility",
        "scope": (
            "Secondary external association; not perturbation ground truth and not "
            "a current-stack pathway score."
        ),
        "sample_count_locked": 502,
        "pathways": PATHWAY_ORDER,
        "clinical_covariates": ["age", "sex", "pathologic_stage", "ever_smoker"],
        "expression_baseline": "mean cohort percentile across unique mapped genes",
        "cv_folds": args.cv_folds,
        "cv_repeats": args.cv_repeats,
        "seed": args.seed,
        "git": git_revision(Path(__file__).resolve().parent),
        "inputs": {
            "activity": {"name": args.activity_tsv.name, "sha256": sha256_file(args.activity_tsv)},
            "clinical": {"name": args.clinical_tsv.name, "sha256": sha256_file(args.clinical_tsv)},
            "expression": {
                "name": args.expression_tsv.name,
                "sha256": sha256_file(args.expression_tsv),
            },
            "root_audit": {
                "name": args.root_audit_tsv.name,
                "sha256": sha256_file(args.root_audit_tsv),
            },
        },
        "script_sha256": sha256_file(Path(__file__)),
    }
    (args.output_dir / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


if __name__ == "__main__":
    main()
