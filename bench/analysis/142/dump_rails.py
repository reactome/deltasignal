#!/usr/bin/env python3
"""specs/142 parity gate, bench side: the rail flag on every scored case fold.

Reads the frozen build's recorded benchmark results (curator and experimental
cases TSVs) and computes, for every valid case's recorded fold (pred_ui), the
bench's own rail test — railed() IMPORTED from bench/analysis/score_124.py,
the one function behind every recorded rails% figure, never a reimplementation
(specs/142 FR-003). Writes bench_rails.json: the folds and the bench's flags,
which rails_142.jl replays through the engine's served function and
compare_142.py diffs.

Usage:
  python3 bench/analysis/142/dump_rails.py RESULTS_DIR [--out bench_rails.json]
where RESULTS_DIR is builds/<id>/results/<solver-sha>/ holding
curator_cases.tsv and experimental_cases.tsv (the frozen record is
20261026-1349_8b6527e / 05dc867).
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from score_124 import railed  # noqa: E402  the bench's own function

CEILING = 100.0  # the frozen build ran the default DS_BASELINE=0.01


def cases(path):
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["valid"] != "1":
                continue
            yield row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--out", default=str(Path(__file__).parent / "bench_rails.json"))
    args = ap.parse_args()
    results_dir = Path(args.results_dir)

    out = {"results_dir": str(results_dir), "ceiling": CEILING, "cases": []}
    for axis in ("curator", "experimental"):
        tsv = results_dir / f"{axis}_cases.tsv"
        for row in cases(tsv):
            fold = float(row["pred_ui"])
            out["cases"].append({
                "axis": axis,
                "pathway": row["pathway"],
                "gene": row["gene"],
                "direction": row["direction"],
                "key_output": row["key_output"],
                "fold": fold,
                "railed": bool(railed(fold, CEILING)),
            })

    Path(args.out).write_text(json.dumps(out, indent=1))
    n = len(out["cases"])
    n_railed = sum(c["railed"] for c in out["cases"])
    print(f"{n} valid case folds dumped ({n_railed} railed) -> {args.out}")


if __name__ == "__main__":
    main()
