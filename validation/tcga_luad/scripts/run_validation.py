#!/usr/bin/env python3
"""Run the TCGA LUAD DeltaSignal validation pipeline.

This script is intentionally strict about inputs and failures:

- consumes fixed LNG outputs and a frozen observability audit;
- parses each LNG pathway into DeltaSignal JSON;
- scores a TCGA LUAD sample cohort with real DeltaSignal solves;
- compares terminal aggregation choices;
- computes median-split survival/log-rank sanity checks;
- computes sample-label shuffled survival controls;
- records solve/mapping failures rather than silently skipping them.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_PATHWAYS = [
    "R-HSA-1257604",
    "R-HSA-453279",
    "R-HSA-69620",
]

SCRIPT_DIR = Path(__file__).resolve().parent
TCGA_VALIDATION_DIR = SCRIPT_DIR.parent
DELTASIGNAL_REPO = TCGA_VALIDATION_DIR.parents[1]


@dataclass(frozen=True)
class PathwaySpec:
    pathway_id: str
    pathway_name: str
    pathway_dir: str
    suitability: str
    logic_network: Path
    uuid_map: Path
    parsed_network: Path


@dataclass
class SolveRecord:
    sample_id: str
    pathway_id: str
    pathway_name: str
    observed_root_uuid_count: int
    terminal_uuid_count: int
    terminal_mean_activity: float
    terminal_median_activity: float
    terminal_q25_activity: float
    terminal_q75_activity: float
    survival_days: float
    survival_event: int
    vital_status: str
    converged: bool | None
    iterations: int | None
    final_residual: float | None
    solve_time_reported: float | None
    solve_wall_seconds: float
    result_json: str
    observations_csv: str
    solve_log: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--lng-output-root",
        type=Path,
        required=True,
        help=(
            "Directory containing Logic Network Generator pathway output folders. "
            "Each pathway folder must contain logic_network.csv and stid_to_uuid_mapping.csv."
        ),
    )
    parser.add_argument(
        "--observability-summary",
        type=Path,
        required=True,
        help=(
            "TSV with at least pathway_id, pathway_name, and suitability columns. "
            "This is produced by the pathway observability audit."
        ),
    )
    parser.add_argument(
        "--root-audit",
        type=Path,
        required=True,
        help=(
            "TSV with root UUID to gene observability rows. "
            "Required columns include pathway_id, root_uuid, mapping_status, and gene_symbol."
        ),
    )
    parser.add_argument(
        "--expression-tsv",
        type=Path,
        required=True,
        help="TCGA LUAD expression matrix TSV with genes as rows and sample IDs as columns.",
    )
    parser.add_argument(
        "--clinical-tsv",
        type=Path,
        required=True,
        help=(
            "TCGA LUAD clinical matrix TSV containing sampleID, days_to_death, "
            "days_to_last_followup, and vital_status."
        ),
    )
    parser.add_argument(
        "--deltasignal-dir",
        type=Path,
        default=DELTASIGNAL_REPO,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=TCGA_VALIDATION_DIR / "outputs" / "tcga_luad_three_pathways",
    )
    parser.add_argument(
        "--pathway-id",
        action="append",
        dest="pathway_ids",
        help="Reactome pathway ID to include. May be repeated. Defaults to three good candidates.",
    )
    parser.add_argument(
        "--sample-count",
        type=int,
        default=80,
        help="Number of eligible TCGA tumor samples to score. Use 0 for all eligible samples.",
    )
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--shuffle-count", type=int, default=200)
    parser.add_argument("--max-iters", type=int, default=2000)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument(
        "--solve-mode",
        choices=["batch", "per-sample"],
        default="batch",
        help="Use one Julia process per pathway, or one Julia process per sample.",
    )
    parser.add_argument("--force-parse", action="store_true")
    parser.add_argument("--force-solve", action="store_true")
    return parser.parse_args()


def normalize_text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", "None", "nan", "NaN", "<NA>"}:
        return None
    return text


def prepare_output_dirs(output_dir: Path) -> dict[str, Path]:
    dirs = {
        "parsed": output_dir / "parsed_networks",
        "observations": output_dir / "observations",
        "results": output_dir / "solve_results",
        "logs": output_dir / "logs",
        "tables": output_dir / "tables",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def load_pathway_specs(args: argparse.Namespace, parsed_dir: Path) -> list[PathwaySpec]:
    requested = set(args.pathway_ids or DEFAULT_PATHWAYS)
    summary = pd.read_csv(args.observability_summary, sep="\t")
    specs: list[PathwaySpec] = []
    for row in summary.itertuples(index=False):
        pathway_id = str(row.pathway_id)
        if pathway_id not in requested:
            continue
        pathway_dir = find_pathway_dir(args.lng_output_root, pathway_id)
        logic_network = pathway_dir / "logic_network.csv"
        uuid_map = pathway_dir / "stid_to_uuid_mapping.csv"
        missing = [str(path) for path in [logic_network, uuid_map] if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing required LNG outputs for {pathway_id}: {missing}")
        parsed_network = parsed_dir / f"{pathway_id.replace('-', '_')}.parsed_network.json"
        specs.append(
            PathwaySpec(
                pathway_id=pathway_id,
                pathway_name=str(row.pathway_name),
                pathway_dir=pathway_dir.name,
                suitability=str(row.suitability),
                logic_network=logic_network,
                uuid_map=uuid_map,
                parsed_network=parsed_network,
            )
        )
    missing_ids = requested - {spec.pathway_id for spec in specs}
    if missing_ids:
        raise ValueError(f"Requested pathway IDs not found in suitability table: {sorted(missing_ids)}")
    return specs


def find_pathway_dir(output_root: Path, pathway_id: str) -> Path:
    suffix = "_" + pathway_id
    matches = [path for path in output_root.iterdir() if path.is_dir() and path.name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one output dir for {pathway_id}, found {len(matches)}")
    return matches[0]


def run_command(cmd: list[str], cwd: Path, log_path: Path) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False)
    log_path.write_text(result.stdout + "\n" + result.stderr)
    return result


def parse_pathway(spec: PathwaySpec, deltasignal_dir: Path, log_dir: Path, force: bool) -> None:
    if spec.parsed_network.exists() and not force:
        return
    cmd = [
        "julia",
        "--project=.",
        "cli/deltasignal.jl",
        "parse",
        "--logic",
        str(spec.logic_network),
        "--uuid-map",
        str(spec.uuid_map),
        "--output",
        str(spec.parsed_network),
        "--validate",
    ]
    log_path = log_dir / f"{spec.pathway_id.replace('-', '_')}.parse.log"
    result = run_command(cmd, deltasignal_dir, log_path)
    if result.returncode != 0:
        raise RuntimeError(f"DeltaSignal parse failed for {spec.pathway_id}. See {log_path}")


def load_network_sets(parsed_network: Path) -> tuple[list[str], list[str]]:
    network = json.loads(parsed_network.read_text())
    parents = {edge["parent_uuid"] for edge in network["edges"]}
    children = {edge["child_uuid"] for edge in network["edges"]}
    return sorted(parents - children), sorted(children - parents)


def load_root_mappings(root_audit_path: Path, pathway_id: str) -> pd.DataFrame:
    audit = pd.read_csv(root_audit_path, sep="\t")
    sub = audit[audit["pathway_id"] == pathway_id].copy()
    if sub.empty:
        raise ValueError(f"No root audit rows found for {pathway_id}")
    return sub


def load_expression(expression_path: Path, genes: set[str]) -> pd.DataFrame:
    expr = pd.read_csv(expression_path, sep="\t")
    first_col = expr.columns[0]
    expr = expr.rename(columns={first_col: "gene"})
    expr = expr[expr["gene"].isin(genes)].copy()
    expr = expr.set_index("gene")
    tumor_cols = [col for col in expr.columns if col.endswith("-01")]
    expr = expr[tumor_cols]
    if expr.empty:
        raise ValueError("No requested genes were found in the expression matrix.")
    return expr


def load_clinical_survival(clinical_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(clinical_path, sep="\t")
    required = {"sampleID", "days_to_death", "days_to_last_followup", "vital_status"}
    missing = required - set(clinical.columns)
    if missing:
        raise ValueError(f"Clinical file is missing required columns: {sorted(missing)}")
    clinical = clinical[list(required)].rename(columns={"sampleID": "sample_id"}).copy()
    clinical["days_to_death"] = pd.to_numeric(clinical["days_to_death"], errors="coerce")
    clinical["days_to_last_followup"] = pd.to_numeric(
        clinical["days_to_last_followup"], errors="coerce"
    )
    clinical["survival_days"] = clinical["days_to_death"].fillna(clinical["days_to_last_followup"])
    clinical["survival_event"] = clinical["vital_status"].map({"DECEASED": 1, "LIVING": 0})
    clinical = clinical.dropna(subset=["sample_id", "survival_days", "survival_event"])
    clinical = clinical[clinical["survival_days"] > 0].copy()
    clinical["survival_event"] = clinical["survival_event"].astype(int)
    return clinical.drop_duplicates(subset=["sample_id"]).set_index("sample_id")


def select_samples(expression: pd.DataFrame, clinical: pd.DataFrame, sample_count: int) -> list[str]:
    eligible = sorted(set(expression.columns).intersection(clinical.index))
    if sample_count > 0:
        eligible = eligible[:sample_count]
    if not eligible:
        raise ValueError("No eligible expression/clinical-overlap samples found.")
    return eligible


def build_observation_rows(
    root_mapping: pd.DataFrame,
    percentiles: pd.DataFrame,
    sample_id: str,
) -> list[dict[str, object]]:
    mapped = root_mapping[root_mapping["mapping_status"] == "mapped"]
    rows: list[dict[str, object]] = []
    for row in mapped.itertuples(index=False):
        gene_symbol = normalize_text(row.gene_symbol)
        if gene_symbol is None or gene_symbol not in percentiles.index:
            continue
        rows.append(
            {
                "node_uuid": row.root_uuid,
                "activity": round(float(percentiles.at[gene_symbol, sample_id]), 6),
                "confidence": 1.0,
            }
        )
    return rows


def write_observations(rows: list[dict[str, object]], output_path: Path) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["node_uuid", "activity", "confidence"])
        writer.writeheader()
        writer.writerows(rows)


def solve_sample(
    spec: PathwaySpec,
    sample_id: str,
    observations_path: Path,
    result_path: Path,
    log_path: Path,
    deltasignal_dir: Path,
    force: bool,
    max_iters: int,
    tolerance: float,
) -> tuple[dict[str, object], float]:
    if result_path.exists() and not force:
        return json.loads(result_path.read_text()), 0.0
    cmd = [
        "julia",
        "--project=.",
        "cli/deltasignal.jl",
        "solve",
        "--network",
        str(spec.parsed_network),
        "--observations",
        str(observations_path),
        "--output",
        str(result_path),
        "--max-iters",
        str(max_iters),
        "--tolerance",
        str(tolerance),
    ]
    start = time.time()
    result = run_command(cmd, deltasignal_dir, log_path)
    wall = time.time() - start
    if result.returncode != 0:
        raise RuntimeError(f"Solve failed for {spec.pathway_id} {sample_id}. See {log_path}")
    return json.loads(result_path.read_text()), wall


def run_batch_solver(
    spec: PathwaySpec,
    manifest_path: Path,
    summary_path: Path,
    log_path: Path,
    deltasignal_dir: Path,
    max_iters: int,
    tolerance: float,
) -> dict[str, object]:
    script_path = Path(__file__).resolve().with_name("batch_solve_deltasignal.jl")
    cmd = [
        "julia",
        str(script_path),
        "--network",
        str(spec.parsed_network),
        "--manifest",
        str(manifest_path),
        "--summary",
        str(summary_path),
        "--max-iters",
        str(max_iters),
        "--tolerance",
        str(tolerance),
    ]
    env = dict(**os.environ, DELTASIGNAL_DIR=str(deltasignal_dir))
    result = subprocess.run(
        cmd,
        cwd=deltasignal_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    log_path.write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"Batch solve failed for {spec.pathway_id}. See {log_path}")
    return json.loads(summary_path.read_text())


def terminal_summary(result_json: dict[str, object], terminals: Iterable[str]) -> dict[str, float]:
    activities = result_json["node_activities"]
    values = np.array([float(activities[uuid]) for uuid in terminals], dtype=float)
    return {
        "terminal_mean_activity": float(np.mean(values)),
        "terminal_median_activity": float(np.median(values)),
        "terminal_q25_activity": float(np.quantile(values, 0.25)),
        "terminal_q75_activity": float(np.quantile(values, 0.75)),
    }


def solver_info(result_json: dict[str, object]) -> dict[str, object]:
    info = result_json.get("solver_info", {})
    if not isinstance(info, dict):
        return {}
    return info


def score_pathway(
    spec: PathwaySpec,
    root_mapping: pd.DataFrame,
    percentiles: pd.DataFrame,
    clinical: pd.DataFrame,
    samples: list[str],
    dirs: dict[str, Path],
    deltasignal_dir: Path,
    force_solve: bool,
    max_iters: int,
    tolerance: float,
    solve_mode: str,
) -> tuple[list[SolveRecord], list[dict[str, object]]]:
    _, terminals = load_network_sets(spec.parsed_network)
    pathway_obs_dir = dirs["observations"] / spec.pathway_id
    pathway_res_dir = dirs["results"] / spec.pathway_id
    pathway_log_dir = dirs["logs"] / spec.pathway_id
    for path in [pathway_obs_dir, pathway_res_dir, pathway_log_dir]:
        path.mkdir(exist_ok=True)

    records: list[SolveRecord] = []
    failures: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []
    for sample_id in samples:
        obs_rows = build_observation_rows(root_mapping, percentiles, sample_id)
        obs_path = pathway_obs_dir / f"{sample_id}.observations.csv"
        result_path = pathway_res_dir / f"{sample_id}.result.json"
        log_path = pathway_log_dir / f"{sample_id}.solve.log"
        write_observations(obs_rows, obs_path)
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "observations_csv": str(obs_path),
                "result_json": str(result_path),
                "solve_log": str(log_path),
            }
        )
        if solve_mode == "batch":
            continue

        try:
            result_json, wall_seconds = solve_sample(
                spec,
                sample_id,
                obs_path,
                result_path,
                log_path,
                deltasignal_dir,
                force_solve,
                max_iters,
                tolerance,
            )
            summary = terminal_summary(result_json, terminals)
            info = solver_info(result_json)
            survival = clinical.loc[sample_id]
            records.append(
                SolveRecord(
                    sample_id=sample_id,
                    pathway_id=spec.pathway_id,
                    pathway_name=spec.pathway_name,
                    observed_root_uuid_count=len(obs_rows),
                    terminal_uuid_count=len(terminals),
                    survival_days=float(survival["survival_days"]),
                    survival_event=int(survival["survival_event"]),
                    vital_status=str(survival["vital_status"]),
                    converged=bool(info["converged"]) if "converged" in info else None,
                    iterations=int(info["iterations"]) if "iterations" in info else None,
                    final_residual=float(info["final_residual"]) if "final_residual" in info else None,
                    solve_time_reported=(
                        float(info["solve_time"]) if "solve_time" in info else None
                    ),
                    solve_wall_seconds=wall_seconds,
                    result_json=str(result_path),
                    observations_csv=str(obs_path),
                    solve_log=str(log_path),
                    **summary,
                )
            )
        except Exception as exc:
            failures.append(
                {
                    "sample_id": sample_id,
                    "pathway_id": spec.pathway_id,
                    "pathway_name": spec.pathway_name,
                    "observations_csv": str(obs_path),
                    "result_json": str(result_path),
                    "solve_log": str(log_path),
                    "error": str(exc),
                }
            )

    if solve_mode == "batch":
        manifest_path = pathway_log_dir / "batch_manifest.tsv"
        batch_summary_path = pathway_log_dir / "batch_summary.json"
        batch_log_path = pathway_log_dir / "batch_solver.log"
        pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
        batch_summary: dict[str, object] | None = None
        try:
            if force_solve or not all(Path(row["result_json"]).exists() for row in manifest_rows):
                batch_summary = run_batch_solver(
                    spec,
                    manifest_path,
                    batch_summary_path,
                    batch_log_path,
                    deltasignal_dir,
                    max_iters,
                    tolerance,
                )
            elif batch_summary_path.exists():
                batch_summary = json.loads(batch_summary_path.read_text())
        except Exception as exc:
            for row in manifest_rows:
                failures.append(
                    {
                        "sample_id": row["sample_id"],
                        "pathway_id": spec.pathway_id,
                        "pathway_name": spec.pathway_name,
                        "observations_csv": row["observations_csv"],
                        "result_json": row["result_json"],
                        "solve_log": row["solve_log"],
                        "error": str(exc),
                    }
                )
            return records, failures

        wall_by_result: dict[str, float] = {}
        if batch_summary is not None:
            for row in batch_summary.get("results", []):
                if isinstance(row, dict):
                    wall_by_result[str(row.get("result_json"))] = float(row.get("wall_seconds", 0.0))

        for row in manifest_rows:
            result_path = Path(row["result_json"])
            obs_path = Path(row["observations_csv"])
            log_path = Path(row["solve_log"])
            sample_id = str(row["sample_id"])
            try:
                if not result_path.exists():
                    raise RuntimeError(f"Batch result JSON was not created: {result_path}")
                result_json = json.loads(result_path.read_text())
                summary = terminal_summary(result_json, terminals)
                info = solver_info(result_json)
                survival = clinical.loc[sample_id]
                obs_count = max(sum(1 for _ in obs_path.open()) - 1, 0)
                records.append(
                    SolveRecord(
                        sample_id=sample_id,
                        pathway_id=spec.pathway_id,
                        pathway_name=spec.pathway_name,
                        observed_root_uuid_count=obs_count,
                        terminal_uuid_count=len(terminals),
                        survival_days=float(survival["survival_days"]),
                        survival_event=int(survival["survival_event"]),
                        vital_status=str(survival["vital_status"]),
                        converged=bool(info["converged"]) if "converged" in info else None,
                        iterations=int(info["iterations"]) if "iterations" in info else None,
                        final_residual=float(info["final_residual"]) if "final_residual" in info else None,
                        solve_time_reported=(
                            float(info["solve_time"]) if "solve_time" in info else None
                        ),
                        solve_wall_seconds=wall_by_result.get(str(result_path), 0.0),
                        result_json=str(result_path),
                        observations_csv=str(obs_path),
                        solve_log=str(log_path),
                        **summary,
                    )
                )
            except Exception as exc:
                failures.append(
                    {
                        "sample_id": sample_id,
                        "pathway_id": spec.pathway_id,
                        "pathway_name": spec.pathway_name,
                        "observations_csv": str(obs_path),
                        "result_json": str(result_path),
                        "solve_log": str(log_path),
                        "error": str(exc),
                    }
                )
    return records, failures


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


def km_median_survival(durations: Iterable[float], events: Iterable[int]) -> float | None:
    df = pd.DataFrame({"duration": list(durations), "event": list(events)}).sort_values("duration")
    survival = 1.0
    for event_time in sorted(df.loc[df["event"] == 1, "duration"].unique()):
        at_risk = int((df["duration"] >= event_time).sum())
        events_at_time = int(((df["duration"] == event_time) & (df["event"] == 1)).sum())
        if at_risk == 0:
            continue
        survival *= 1.0 - (events_at_time / at_risk)
        if survival <= 0.5:
            return float(event_time)
    return None


def survival_split_stats(activity: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for pathway_id, sub in activity.groupby("pathway_id"):
        values = sub[metric].astype(float)
        threshold = float(values.median())
        groups = (values > threshold).astype(int)
        if groups.nunique() < 2:
            continue
        stats = logrank_statistic(sub["survival_days"], sub["survival_event"], groups)
        low = sub[groups == 0]
        high = sub[groups == 1]
        rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": sub["pathway_name"].iloc[0],
                "metric": metric,
                "threshold": threshold,
                "low_n": int(len(low)),
                "high_n": int(len(high)),
                "low_events": int(low["survival_event"].sum()),
                "high_events": int(high["survival_event"].sum()),
                "low_km_median_survival": km_median_survival(low["survival_days"], low["survival_event"]),
                "high_km_median_survival": km_median_survival(high["survival_days"], high["survival_event"]),
                **stats,
            }
        )
    return pd.DataFrame(rows)


def shuffled_survival_controls(
    activity: pd.DataFrame,
    metric: str,
    shuffle_count: int,
    seed: int,
) -> pd.DataFrame:
    rng = random.Random(seed)
    rows: list[dict[str, object]] = []
    for pathway_id, sub in activity.groupby("pathway_id"):
        true_stats = survival_split_stats(sub, metric)
        true_p = float(true_stats["logrank_p"].iloc[0]) if not true_stats.empty else float("nan")
        p_values: list[float] = []
        values = sub[metric].astype(float).tolist()
        for _ in range(shuffle_count):
            shuffled = values[:]
            rng.shuffle(shuffled)
            threshold = float(median(shuffled))
            groups = [int(value > threshold) for value in shuffled]
            stats = logrank_statistic(sub["survival_days"], sub["survival_event"], groups)
            p_values.append(float(stats["logrank_p"]))
        rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": sub["pathway_name"].iloc[0],
                "metric": metric,
                "shuffle_count": shuffle_count,
                "observed_logrank_p": true_p,
                "shuffle_p_median": float(np.median(p_values)) if p_values else float("nan"),
                "shuffle_p05": float(np.quantile(p_values, 0.05)) if p_values else float("nan"),
                "shuffle_p95": float(np.quantile(p_values, 0.95)) if p_values else float("nan"),
                "shuffle_fraction_p_le_observed": (
                    float(np.mean(np.array(p_values) <= true_p)) if p_values else float("nan")
                ),
            }
        )
    return pd.DataFrame(rows)


def write_key_findings(
    activity: pd.DataFrame,
    mapping: pd.DataFrame,
    survival_stats: pd.DataFrame,
    shuffle_controls: pd.DataFrame,
    output_path: Path,
) -> None:
    rows: list[dict[str, object]] = []
    sample_count = int(activity["sample_id"].nunique())
    for pathway_id, sub in activity.groupby("pathway_id"):
        mapping_row = mapping[mapping["pathway_id"] == pathway_id].iloc[0]
        if {"pathway_id", "metric"}.issubset(survival_stats.columns):
            survival_row = survival_stats[
                (survival_stats["pathway_id"] == pathway_id)
                & (survival_stats["metric"] == "terminal_mean_activity")
            ]
        else:
            survival_row = pd.DataFrame()
        if "pathway_id" in shuffle_controls.columns:
            shuffle_row = shuffle_controls[shuffle_controls["pathway_id"] == pathway_id]
        else:
            shuffle_row = pd.DataFrame()
        rows.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": sub["pathway_name"].iloc[0],
                "samples_scored": int(len(sub)),
                "solve_failures": sample_count - int(len(sub)),
                "converged": int((sub["converged"] == True).sum()),
                "nonconverged": int((sub["converged"] != True).sum()),
                "mapped_roots": int(mapping_row["mapped"]),
                "stable_id_missing": int(mapping_row["stable_id_missing"]),
                "mean_reported_solve_time": float(sub["solve_time_reported"].dropna().mean()),
                "mean_wall_seconds": float(sub["solve_wall_seconds"].dropna().mean()),
                "terminal_mean_min": float(sub["terminal_mean_activity"].min()),
                "terminal_mean_median": float(sub["terminal_mean_activity"].median()),
                "terminal_mean_max": float(sub["terminal_mean_activity"].max()),
                "median_split_logrank_p": (
                    float(survival_row["logrank_p"].iloc[0])
                    if not survival_row.empty
                    else float("nan")
                ),
                "shuffle_fraction_p_le_observed": (
                    float(shuffle_row["shuffle_fraction_p_le_observed"].iloc[0])
                    if not shuffle_row.empty
                    else float("nan")
                ),
            }
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def write_mapping_readiness(root_audit_path: Path, specs: list[PathwaySpec], output_path: Path) -> None:
    audit = pd.read_csv(root_audit_path, sep="\t")
    rows: list[dict[str, object]] = []
    for spec in specs:
        sub = audit[audit["pathway_id"] == spec.pathway_id]
        status_counts = sub["mapping_status"].value_counts().to_dict()
        rows.append(
            {
                "pathway_id": spec.pathway_id,
                "pathway_name": spec.pathway_name,
                "root_uuid_count": int(len(sub)),
                "stable_id_missing": int(status_counts.get("stable_id_missing", 0)),
                "mapped": int(status_counts.get("mapped", 0)),
                "gene_not_in_expression": int(status_counts.get("gene_not_in_expression", 0)),
                "no_gene_mapping": int(status_counts.get("no_gene_mapping", 0)),
            }
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    random.seed(args.seed)
    dirs = prepare_output_dirs(args.output_dir)
    specs = load_pathway_specs(args, dirs["parsed"])

    for spec in specs:
        parse_pathway(spec, args.deltasignal_dir, dirs["logs"], args.force_parse)

    root_audits = {spec.pathway_id: load_root_mappings(args.root_audit, spec.pathway_id) for spec in specs}
    mapped_genes = {
        normalize_text(row.gene_symbol)
        for root_mapping in root_audits.values()
        for row in root_mapping[root_mapping["mapping_status"] == "mapped"].itertuples(index=False)
    }
    mapped_genes = {gene for gene in mapped_genes if gene is not None}
    expression = load_expression(args.expression_tsv, mapped_genes)
    clinical = load_clinical_survival(args.clinical_tsv)
    samples = select_samples(expression, clinical, args.sample_count)
    percentiles = expression.rank(axis=1, method="average", pct=True) * 100.0

    mapping_readiness_path = dirs["tables"] / "mapping_readiness.tsv"
    write_mapping_readiness(args.root_audit, specs, mapping_readiness_path)

    all_records: list[SolveRecord] = []
    all_failures: list[dict[str, object]] = []
    for spec in specs:
        records, failures = score_pathway(
            spec,
            root_audits[spec.pathway_id],
            percentiles,
            clinical,
            samples,
            dirs,
            args.deltasignal_dir,
            args.force_solve,
            args.max_iters,
            args.tolerance,
            args.solve_mode,
        )
        all_records.extend(records)
        all_failures.extend(failures)

    activity = pd.DataFrame([asdict(record) for record in all_records])
    if activity.empty:
        raise RuntimeError("No successful DeltaSignal solves were produced.")
    activity.to_csv(dirs["tables"] / "sample_pathway_activity.tsv", sep="\t", index=False)
    pd.DataFrame(all_failures).to_csv(dirs["tables"] / "solve_failures.tsv", sep="\t", index=False)

    matrix = activity.pivot(index="sample_id", columns="pathway_id", values="terminal_mean_activity")
    matrix.to_csv(dirs["tables"] / "pathway_activity_matrix_mean.tsv", sep="\t")

    metric_rows = []
    metrics = [
        "terminal_mean_activity",
        "terminal_median_activity",
        "terminal_q25_activity",
        "terminal_q75_activity",
    ]
    for metric in metrics:
        metric_stats = survival_split_stats(activity, metric)
        if not metric_stats.empty:
            metric_rows.append(metric_stats)
    survival_stats = pd.concat(metric_rows, ignore_index=True) if metric_rows else pd.DataFrame()
    survival_stats.to_csv(dirs["tables"] / "survival_median_split_stats.tsv", sep="\t", index=False)

    shuffled = shuffled_survival_controls(
        activity,
        metric="terminal_mean_activity",
        shuffle_count=args.shuffle_count,
        seed=args.seed,
    )
    shuffled.to_csv(dirs["tables"] / "sample_label_shuffle_controls.tsv", sep="\t", index=False)
    mapping_readiness = pd.read_csv(mapping_readiness_path, sep="\t")
    write_key_findings(
        activity,
        mapping_readiness,
        survival_stats,
        shuffled,
        dirs["tables"] / "readiness_key_findings.tsv",
    )

    summary = {
        "sample_count_requested": args.sample_count,
        "sample_count_scored": len(samples),
        "samples": samples,
        "pathways": [asdict(spec) | {
            "logic_network": str(spec.logic_network),
            "uuid_map": str(spec.uuid_map),
            "parsed_network": str(spec.parsed_network),
        } for spec in specs],
        "solve_record_count": len(all_records),
        "solve_failure_count": len(all_failures),
        "solver_max_iters": args.max_iters,
        "solver_tolerance": args.tolerance,
        "solve_mode": args.solve_mode,
        "all_solves_converged": bool(activity["converged"].fillna(False).all()),
        "mean_reported_solve_time": float(activity["solve_time_reported"].dropna().mean()),
        "max_reported_solve_time": float(activity["solve_time_reported"].dropna().max()),
        "outputs": {name: str(path) for name, path in dirs.items()},
    }
    (args.output_dir / "readiness_run_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
