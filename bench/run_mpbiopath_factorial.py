#!/usr/bin/env python3
"""Run the perturbation benchmark across frozen graph and solver conditions."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


SOLVER_CONFIGS = {
    "legacy_flat": {
        "DS_SCC_SOLVE": "0",
        "DS_INHIBITION_MODE": "spec",
        "DS_AND_MODE": "geomean",
        "DS_OR_MODE": "max",
        "DS_ASSEMBLY_LIMITING": "0",
        "DS_HILL_LOG_ZMAX": "4.605170185988092",
    },
    "legacy_scc": {
        "DS_SCC_SOLVE": "1",
        "DS_INHIBITION_MODE": "spec",
        "DS_AND_MODE": "geomean",
        "DS_OR_MODE": "max",
        "DS_ASSEMBLY_LIMITING": "0",
        "DS_HILL_LOG_ZMAX": "4.605170185988092",
    },
    "current_flat": {
        "DS_SCC_SOLVE": "0",
        "DS_INHIBITION_MODE": "divide",
        "DS_AND_MODE": "hill_log",
        "DS_OR_MODE": "mean",
        "DS_ASSEMBLY_LIMITING": "1",
        "DS_HILL_LOG_ZMAX": "10.0",
    },
    "current_scc": {
        "DS_SCC_SOLVE": "1",
        "DS_INHIBITION_MODE": "divide",
        "DS_AND_MODE": "hill_log",
        "DS_OR_MODE": "mean",
        "DS_ASSEMBLY_LIMITING": "1",
        "DS_HILL_LOG_ZMAX": "10.0",
    },
}


def named_path(value: str) -> tuple[str, Path]:
    try:
        name, path = value.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected NAME=/path/to/catalog") from exc
    return name, Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", action="append", type=named_path, required=True)
    parser.add_argument("--supplementary-workbook", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True)
    parser.add_argument("--reactome-id-audit", type=Path, required=True)
    parser.add_argument(
        "--deltasignal-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--base-port", type=int, default=18100)
    parser.add_argument("--allow-output-proxies", action="store_true")
    return parser.parse_args()


def flatten(condition: str, summary: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    scopes = {
        "all": summary,
        **summary["by_split"],
    }
    for scope, metrics in scopes.items():
        rows.append(
            {
                "condition": condition,
                "scope": scope,
                "eligible_cases": metrics["eligible_cases"],
                "scored_cases": metrics["scored_cases"],
                "correct": metrics["correct"],
                "accuracy": metrics["accuracy"],
                "coverage_adjusted_accuracy": metrics[
                    "coverage_adjusted_accuracy"
                ],
                "macro_f1": metrics["macro_f1"],
                "balanced_accuracy": metrics["balanced_accuracy"],
                "change_f1": metrics["change_f1"],
                "nonconverged_solves": metrics["nonconverged_solves"],
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    benchmark = Path(__file__).with_name("benchmark_mpbiopath_cases.py")
    consolidated: list[dict[str, object]] = []
    run_index = 0
    for catalog_name, catalog_path in args.catalog:
        for solver_name, overrides in SOLVER_CONFIGS.items():
            condition = f"{catalog_name}__{solver_name}"
            output = args.output_dir / condition
            command = [
                sys.executable,
                str(benchmark),
                "--supplementary-workbook",
                str(args.supplementary_workbook),
                "--id-map",
                str(args.id_map),
                "--reactome-id-audit",
                str(args.reactome_id_audit),
                "--catalog",
                str(catalog_path),
                "--deltasignal-dir",
                str(args.deltasignal_dir),
                "--output-dir",
                str(output),
                "--bootstrap",
                str(args.bootstrap),
                "--port",
                str(args.base_port + run_index),
            ]
            if args.allow_output_proxies:
                command.append("--allow-output-proxies")
            env = os.environ.copy()
            for key in tuple(env):
                if key.startswith("DS_"):
                    del env[key]
            env.update(overrides)
            subprocess.run(command, env=env, check=True)
            summary = json.loads((output / "benchmark_summary.json").read_text())
            consolidated.extend(flatten(condition, summary))
            run_index += 1

    output = args.output_dir / "factorial_summary.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(consolidated[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(consolidated)
    print(output)


if __name__ == "__main__":
    main()
