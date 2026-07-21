#!/usr/bin/env python3
"""Benchmark DeltaSignal against fixed MP-BioPath perturbation cases.

The benchmark keeps the network catalog fixed while solver configurations or
DeltaSignal commits change. It joins the published supplementary case table to
LNG UUIDs through Reactome database identifiers and reports all mapping losses.
Missing mappings are never converted into synthetic NORMAL predictions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import re
import signal
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


DOWN, NORMAL, UP = 0, 1, 2
LABELS = (DOWN, NORMAL, UP)
STATUS_TO_CLASS = {
    "downregulation/inactivation": DOWN,
    "unchanged": NORMAL,
    "upregulation/activation": UP,
}
DEFAULT_PATHWAYS = ("R-HSA-1257604", "R-HSA-453279", "R-HSA-69620")


@dataclass(frozen=True)
class Case:
    pathway_id: str
    pathway_name: str
    gene: str
    direction: int
    key_output_dbid: str
    expected: int
    curator_prediction: int
    mpbiopath_prediction: int
    evidence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplementary-workbook", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--deltasignal-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pathway-id", action="append", dest="pathway_ids")
    parser.add_argument("--ground-truth", choices=("experimental", "curator"), default="experimental")
    parser.add_argument("--down-cutoff", type=float, default=0.85)
    parser.add_argument("--up-cutoff", type=float, default=1.15)
    parser.add_argument("--output-aggregation", choices=("max", "mean", "min", "extreme"), default="max")
    parser.add_argument("--perturbation-up", type=float, default=80.0)
    parser.add_argument("--port", type=int, default=18080)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--server-timeout", type=float, default=180.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state(path: Path) -> dict[str, object]:
    def run(*args: str) -> str:
        return subprocess.check_output(
            ["git", "-C", str(path), *args], text=True, stderr=subprocess.DEVNULL
        ).strip()

    try:
        status = run("status", "--porcelain").splitlines()
        return {
            "commit": run("rev-parse", "HEAD"),
            "branch": run("branch", "--show-current"),
            "dirty": bool(status),
            "status_porcelain": status,
        }
    except (OSError, subprocess.CalledProcessError):
        return {"available": False}


def normalize_cell(value: object) -> str:
    return "" if value is None else str(value).strip()


def load_cases(workbook: Path, ground_truth: str) -> list[Case]:
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise RuntimeError("openpyxl is required to read the supplementary workbook") from exc

    sheet = load_workbook(workbook, read_only=True, data_only=True)["Supplementary Table S1"]
    rows = sheet.iter_rows(values_only=True)
    header = next(rows)
    index = {str(name): idx for idx, name in enumerate(header)}
    required = {
        "Reactome_pathway_ID",
        "Reactome_pathway_name",
        "RootInput_Perturbation",
        "KeyOutput_ID",
        "CuratorPredictedKeyOutputStatus",
        "MP-BioPathPredictedKeyOutputStatus",
        "ExperimentalKeyOutputStatus",
        "PublishedEvidence",
    }
    missing = sorted(required - set(index))
    if missing:
        raise ValueError(f"Supplementary Table S1 is missing columns: {missing}")

    cases: list[Case] = []
    truth_column = (
        "ExperimentalKeyOutputStatus"
        if ground_truth == "experimental"
        else "CuratorPredictedKeyOutputStatus"
    )
    perturbation_pattern = re.compile(r"^(.+)_(downregulation|upregulation)$")
    for row in rows:
        expected_text = normalize_cell(row[index[truth_column]])
        if expected_text not in STATUS_TO_CLASS:
            continue
        perturbation = normalize_cell(row[index["RootInput_Perturbation"]])
        match = perturbation_pattern.fullmatch(perturbation)
        if match is None:
            raise ValueError(f"Unexpected perturbation label: {perturbation!r}")
        gene, direction_text = match.groups()
        curator_text = normalize_cell(row[index["CuratorPredictedKeyOutputStatus"]])
        mpbiopath_text = normalize_cell(row[index["MP-BioPathPredictedKeyOutputStatus"]])
        if curator_text not in STATUS_TO_CLASS or mpbiopath_text not in STATUS_TO_CLASS:
            raise ValueError(f"Unexpected prediction status in case {perturbation}")
        pathway_numeric = normalize_cell(row[index["Reactome_pathway_ID"]])
        cases.append(
            Case(
                pathway_id=f"R-HSA-{pathway_numeric}",
                pathway_name=normalize_cell(row[index["Reactome_pathway_name"]]),
                gene=gene,
                direction=DOWN if direction_text == "downregulation" else UP,
                key_output_dbid=normalize_cell(row[index["KeyOutput_ID"]]),
                expected=STATUS_TO_CLASS[expected_text],
                curator_prediction=STATUS_TO_CLASS[curator_text],
                mpbiopath_prediction=STATUS_TO_CLASS[mpbiopath_text],
                evidence=normalize_cell(row[index["PublishedEvidence"]]),
            )
        )
    return cases


def load_gene_dbids(id_map: Path, genes: set[str]) -> dict[str, set[str]]:
    mapped: dict[str, set[str]] = defaultdict(set)
    with id_map.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Database_Identifier", "Display_Name", "Reference_Entity_Name"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"ID map is missing columns: {missing}")
        for row in reader:
            names = {normalize_cell(row.get("Display_Name")), normalize_cell(row.get("Reference_Entity_Name"))}
            for gene in names & genes:
                mapped[gene].add(normalize_cell(row["Database_Identifier"]))
    return mapped


def find_pathway_dir(catalog: Path, pathway_id: str) -> Path | None:
    matches = sorted(path for path in catalog.iterdir() if path.is_dir() and path.name.endswith(pathway_id))
    if len(matches) > 1:
        raise ValueError(f"Multiple catalog directories match {pathway_id}: {matches}")
    return matches[0] if matches else None


def load_dbid_to_uuids(pathway_dir: Path) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    nodes_file = pathway_dir / "nodes.csv"
    if nodes_file.exists():
        with nodes_file.open(newline="") as handle:
            for row in csv.DictReader(handle):
                stable_ids = {normalize_cell(row.get("diagram_entity_id"))}
                stable_ids.update(normalize_cell(item) for item in normalize_cell(row.get("member_leaves")).split("|"))
                for stable_id in stable_ids - {""}:
                    dbid = stable_id.rsplit("-", 1)[-1]
                    pair = (dbid, row["uuid"])
                    if pair not in seen:
                        seen.add(pair)
                        mapping[dbid].append(row["uuid"])
    else:
        with (pathway_dir / "stid_to_uuid_mapping.csv").open(newline="") as handle:
            for row in csv.DictReader(handle):
                dbid = row["stable_id"].rsplit("-", 1)[-1]
                pair = (dbid, row["uuid"])
                if pair not in seen:
                    seen.add(pair)
                    mapping[dbid].append(row["uuid"])
    return mapping


def request_json(url: str, payload: dict | None = None, timeout: float = 600.0) -> dict:
    data = None if payload is None else json.dumps(payload).encode()
    request = Request(url, data=data, headers={"Content-Type": "application/json"})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read())


def start_server(args: argparse.Namespace) -> tuple[subprocess.Popen, Path]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.output_dir / "deltasignal_server.log"
    log_handle = log_path.open("w")
    env = os.environ.copy()
    env["DS_PATHWAY_CATALOG"] = str(args.catalog.resolve())
    process = subprocess.Popen(
        [
            "julia",
            f"--project={args.deltasignal_dir.resolve()}",
            str(args.deltasignal_dir.resolve() / "src/api/server.jl"),
            "--port",
            str(args.port),
        ],
        cwd=args.deltasignal_dir,
        env=env,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
    )
    process._deltasignal_log_handle = log_handle  # type: ignore[attr-defined]
    base = f"http://127.0.0.1:{args.port}"
    deadline = time.time() + args.server_timeout
    while time.time() < deadline:
        if process.poll() is not None:
            log_handle.flush()
            raise RuntimeError(f"DeltaSignal server exited early; see {log_path}")
        try:
            if request_json(f"{base}/api/health", timeout=2).get("status") == "ok":
                return process, log_path
        except (HTTPError, URLError, TimeoutError):
            time.sleep(0.5)
    process.terminate()
    raise TimeoutError(f"DeltaSignal server did not become ready; see {log_path}")


def stop_server(process: subprocess.Popen) -> None:
    if process.poll() is None:
        process.send_signal(signal.SIGINT)
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.terminate()
            process.wait(timeout=10)
    log_handle = getattr(process, "_deltasignal_log_handle", None)
    if log_handle is not None:
        log_handle.close()


def aggregate(values: list[float], mode: str) -> float:
    if mode == "mean":
        return sum(values) / len(values)
    if mode == "min":
        return min(values)
    if mode == "extreme":
        return max(values, key=lambda value: abs(value - 1.0))
    return max(values)


def classify(value: float, down_cutoff: float, up_cutoff: float) -> int:
    if value < down_cutoff:
        return DOWN
    if value >= up_cutoff:
        return UP
    return NORMAL


def metric_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    scored = [row for row in rows if row["prediction"] is not None]
    confusion = Counter((int(row["prediction"]), int(row["expected"])) for row in scored)
    per_class = {}
    for label in LABELS:
        tp = confusion[(label, label)]
        fp = sum(confusion[(label, expected)] for expected in LABELS if expected != label)
        fn = sum(confusion[(predicted, label)] for predicted in LABELS if predicted != label)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        per_class[str(label)] = {"precision": precision, "recall": recall, "f1": f1, "support": tp + fn}
    correct = sum(int(row["prediction"]) == int(row["expected"]) for row in scored)
    return {
        "eligible_cases": len(rows),
        "scored_cases": len(scored),
        "coverage": len(scored) / len(rows) if rows else 0.0,
        "correct": correct,
        "accuracy": correct / len(scored) if scored else 0.0,
        "coverage_adjusted_accuracy": correct / len(rows) if rows else 0.0,
        "macro_f1": sum(per_class[str(label)]["f1"] for label in LABELS) / len(LABELS),
        "balanced_accuracy": sum(per_class[str(label)]["recall"] for label in LABELS) / len(LABELS),
        "change_f1": (per_class[str(DOWN)]["f1"] + per_class[str(UP)]["f1"]) / 2,
        "per_class": per_class,
        "confusion": {f"pred_{predicted}_expected_{expected}": count for (predicted, expected), count in sorted(confusion.items())},
        "unscored_reasons": dict(Counter(str(row["mapping_status"]) for row in rows if row["prediction"] is None)),
        "nonconverged_solves": sum(not bool(row["converged"]) for row in scored),
    }


def baseline_accuracy(rows: list[dict[str, object]], field: str) -> dict[str, float | int]:
    scored = [row for row in rows if row["prediction"] is not None]
    correct = sum(int(row[field]) == int(row["expected"]) for row in scored)
    return {"correct": correct, "total": len(scored), "accuracy": correct / len(scored) if scored else 0.0}


def paired_bootstrap(rows: list[dict[str, object]], field: str, iterations: int, seed: int) -> dict[str, float]:
    scored = [row for row in rows if row["prediction"] is not None]
    if not scored or iterations <= 0:
        return {}
    rng = random.Random(seed)
    differences = []
    for _ in range(iterations):
        sample = [scored[rng.randrange(len(scored))] for _ in scored]
        ds = sum(int(row["prediction"]) == int(row["expected"]) for row in sample) / len(sample)
        baseline = sum(int(row[field]) == int(row["expected"]) for row in sample) / len(sample)
        differences.append(ds - baseline)
    differences.sort()
    return {
        "mean_accuracy_difference": sum(differences) / len(differences),
        "ci95_low": differences[math.floor(0.025 * (len(differences) - 1))],
        "ci95_high": differences[math.ceil(0.975 * (len(differences) - 1))],
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    requested = set(args.pathway_ids or DEFAULT_PATHWAYS)
    all_cases = load_cases(args.supplementary_workbook, args.ground_truth)
    cases = [case for case in all_cases if case.pathway_id in requested]
    found_ids = {case.pathway_id for case in cases}
    if requested - found_ids:
        raise ValueError(f"Requested pathways absent from benchmark cases: {sorted(requested - found_ids)}")

    gene_dbids = load_gene_dbids(args.id_map, {case.gene for case in cases})
    pathway_dirs = {pathway_id: find_pathway_dir(args.catalog, pathway_id) for pathway_id in requested}
    absent_networks = sorted(pathway_id for pathway_id, path in pathway_dirs.items() if path is None)
    if absent_networks:
        raise FileNotFoundError(f"No LNG network for pathways: {absent_networks}")
    dbid_maps = {pathway_id: load_dbid_to_uuids(path) for pathway_id, path in pathway_dirs.items() if path}

    process, log_path = start_server(args)
    base = f"http://127.0.0.1:{args.port}"
    parsed_ids = {}
    rows: list[dict[str, object]] = []
    solve_cache: dict[tuple[str, str, int], dict] = {}
    try:
        for pathway_id, pathway_dir in pathway_dirs.items():
            parsed = request_json(f"{base}/api/parse", {"pathway_id": pathway_dir.name})
            if parsed.get("status") != "success":
                raise RuntimeError(f"Parse failed for {pathway_id}: {parsed}")
            parsed_ids[pathway_id] = parsed["network_id"]

        for case in cases:
            dbid_map = dbid_maps[case.pathway_id]
            gene_uuids = sorted({uuid for dbid in gene_dbids.get(case.gene, ()) for uuid in dbid_map.get(dbid, ())})
            output_uuids = sorted(set(dbid_map.get(case.key_output_dbid, ())))
            row = asdict(case)
            row.update(
                {
                    "gene_dbid_count": len(gene_dbids.get(case.gene, ())),
                    "gene_uuid_count": len(gene_uuids),
                    "output_uuid_count": len(output_uuids),
                    "prediction": None,
                    "predicted_ui": None,
                    "converged": None,
                    "iterations": None,
                    "mapping_status": "mapped",
                }
            )
            if not gene_dbids.get(case.gene):
                row["mapping_status"] = "gene_absent_from_legacy_id_map"
            elif not gene_uuids:
                row["mapping_status"] = "gene_dbids_absent_from_network"
            elif not output_uuids:
                row["mapping_status"] = "key_output_dbid_absent_from_network"
            else:
                cache_key = (case.pathway_id, case.gene, case.direction)
                if cache_key not in solve_cache:
                    perturbation = 0.0 if case.direction == DOWN else args.perturbation_up
                    observations = {uuid: [perturbation, 1.0] for uuid in gene_uuids}
                    solve_cache[cache_key] = request_json(
                        f"{base}/api/solve",
                        {"network_id": parsed_ids[case.pathway_id], "observations": observations},
                    )
                result = solve_cache[cache_key]
                if result.get("status") != "success":
                    row["mapping_status"] = "solve_failed"
                else:
                    values = [float(result["node_activities"].get(uuid, 0.01)) * 100.0 for uuid in output_uuids]
                    predicted_ui = aggregate(values, args.output_aggregation)
                    row.update(
                        {
                            "prediction": classify(predicted_ui, args.down_cutoff, args.up_cutoff),
                            "predicted_ui": predicted_ui,
                            "converged": bool(result.get("converged")),
                            "iterations": result.get("iterations"),
                        }
                    )
            rows.append(row)
    finally:
        stop_server(process)

    summary = metric_summary(rows)
    summary["curator_on_scored_cases"] = baseline_accuracy(rows, "curator_prediction")
    summary["mpbiopath_on_scored_cases"] = baseline_accuracy(rows, "mpbiopath_prediction")
    summary["paired_bootstrap_vs_curator"] = paired_bootstrap(rows, "curator_prediction", args.bootstrap, args.seed)
    summary["paired_bootstrap_vs_mpbiopath"] = paired_bootstrap(rows, "mpbiopath_prediction", args.bootstrap, args.seed + 1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    case_path = args.output_dir / "benchmark_cases.tsv"
    with case_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (args.output_dir / "benchmark_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    network_files = {}
    for pathway_id, pathway_dir in pathway_dirs.items():
        network_files[pathway_id] = {
            name: {"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size}
            for name in ("logic_network.csv", "stid_to_uuid_mapping.csv")
            if (path := pathway_dir / name).exists()
        }
    solver_config = {
        key: value for key, value in sorted(os.environ.items()) if key.startswith("DS_")
    }
    solver_config["DS_PATHWAY_CATALOG"] = str(args.catalog.resolve())
    manifest = {
        "created_unix": time.time(),
        "ground_truth": args.ground_truth,
        "pathway_ids": sorted(requested),
        "thresholds": {"down": args.down_cutoff, "up": args.up_cutoff},
        "output_aggregation": args.output_aggregation,
        "perturbation_up": args.perturbation_up,
        "solver_environment_overrides": solver_config,
        "deltasignal_git": git_state(args.deltasignal_dir),
        "inputs": {
            "supplementary_workbook": {"path": str(args.supplementary_workbook), "sha256": sha256(args.supplementary_workbook)},
            "id_map": {"path": str(args.id_map), "sha256": sha256(args.id_map)},
            "networks": network_files,
        },
        "server_log": str(log_path),
        "outputs": {"cases": str(case_path), "summary": str(args.output_dir / "benchmark_summary.json")},
    }
    (args.output_dir / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return summary


def main() -> None:
    args = parse_args()
    summary = run(args)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
