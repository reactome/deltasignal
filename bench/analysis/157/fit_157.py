#!/usr/bin/env python3
"""specs/157: build the served fragility table from the ensemble record.

Descriptive rates, nothing fitted: per class — the registered level walk
full (call|mag|conv) -> no_converged (call|mag) -> call_only (call), the
specs/146 magnitude bin on the canonical fold — the share of labelled ensemble
pairs that were FRAGILE (any contributing variant flipped the call), pooled
across the three record axes with per-axis counts kept. Provenance (catalog
build, solver commit, spec) travels with the table and the loader refuses a
mismatch (the specs/146/155 pattern).

Usage: fit_157.py RUNS_DIR OUT_JSON [SOLVER_COMMIT]  (reads labels.tsv + calls_canonical.tsv[.gz]; SOLVER_COMMIT names
the code the ensemble runs were solved with, default this checkout's HEAD)
"""
import gzip
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dataset_roles  # noqa: E402
from fold_rules import record_fold, untouched  # noqa: E402

MIN_CELL = 30
BUILD_JSON = Path.home() / "deltasignal-catalogs/builds/20261113-r98fix-plus-leads/BUILD.json"


def mag_bin(fold):
    fold = record_fold(fold)
    if fold == 1.0:
        return "baseline"
    a = abs(math.log2(fold)) if fold > 0 else float("inf")
    return "le1" if a <= 1 else ("le3" if a <= 3 else "gt3")


def opengz(path):
    p = Path(path)
    if p.exists():
        return open(p)
    return gzip.open(str(p) + ".gz", "rt")


def main():
    runs, out = Path(sys.argv[1]), Path(sys.argv[2])
    dataset_roles.require_role("korkut2015", "fit")

    import csv
    canon = {}
    with opengz(runs / "calls_canonical.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (r["axis"], r["pathway"], r["case_a"], r["case_b"], r["readout"])
            canon[key] = (float(r["fold"]) if r["fold"] else None, r["converged"] == "1")

    full = defaultdict(lambda: [0, 0])      # call|mag|conv -> [fragile, n]
    no_conv = defaultdict(lambda: [0, 0])   # call|mag
    call_only = defaultdict(lambda: [0, 0])
    per_axis = Counter()
    per_axis_fragile = Counter()
    skipped_untouched = 0
    with open(runs / "labels.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            key = (r["axis"], r["pathway"], r["case_a"], r["case_b"], r["readout"])
            fold, conv = canon.get(key, (None, None))
            if fold is None:
                continue
            if untouched(fold):
                skipped_untouched += 1   # served nothing (the untouched rule)
                continue
            frag = r["fragile"] == "1"
            call = r["canonical_call"]
            mag = mag_bin(fold)
            cv = "True" if conv else "False"
            for table, k in ((full, f"{call}|{mag}|{cv}"),
                             (no_conv, f"{call}|{mag}"),
                             (call_only, call)):
                table[k][0] += frag
                table[k][1] += 1
            per_axis[r["axis"]] += 1
            per_axis_fragile[r["axis"]] += frag

    def cells(table):
        return {k: {"fragile_rate": round(f / n, 4), "n": n}
                for k, (f, n) in sorted(table.items()) if n >= MIN_CELL}

    build = json.loads(BUILD_JSON.read_text())
    commit = sys.argv[3] if len(sys.argv) > 3 else subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
        cwd=Path(__file__).resolve().parents[2]).stdout.strip()
    doc = {
        # 2026-11-20: refitted with the record-precision untouched rule; 2026-11-21: korkut rows re-run with
        # phospho readouts valued as the record values them (`sites` mode, issue #278; runs_278/)
        "version": "2026-11-21",
        "spec": "specs/157-structural-fragility",
        "catalog_build": build["build_id"],
        "solver_commit": commit,
        "min_cell": MIN_CELL,
        "levels": ["full", "no_converged", "call_only"],
        "untouched_rule": "a fold that reads 1 at 6 significant figures serves no fragility",
        "claim": ("the share of recorded ensemble pairs in this class whose call "
                  "flipped under at least one recorded near-equivalent rule variant "
                  "(the registered 20-variant ensemble) — stability across contested "
                  "structure, never correctness"),
        "validation": {
            "gate": ("specs/157 research.md, amendment 2026-11-21 (korkut re-run, sites mode): curator +0.2818 "
                     "[+0.1298, +0.4573], korkut +0.6016 [+0.5024, +0.6796]"),
            "per_axis_pairs": dict(per_axis),
            "per_axis_fragile": dict(per_axis_fragile),
            "untouched_skipped": skipped_untouched,
        },
        "table": {"full": cells(full), "no_converged": cells(no_conv),
                  "call_only": cells(call_only)},
    }
    out.write_text(json.dumps(doc, indent=1, sort_keys=False) + "\n")
    print(f"{sum(n for _, n in call_only.values())} pairs -> {out} "
          f"({len(doc['table']['full'])} full cells, "
          f"{len(doc['table']['no_converged'])} no_converged, "
          f"{len(doc['table']['call_only'])} call_only; "
          f"{skipped_untouched} untouched skipped)")


if __name__ == "__main__":
    main()
