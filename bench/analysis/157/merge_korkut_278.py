#!/usr/bin/env python3
"""Issue #278: splice re-run korkut rows into the specs/157 ensemble record.

The 2026-11-15 runs valued korkut phospho readouts in `forms` mode; the re-run (ensemble_157.py
--axes korkut, phospho in `sites` mode) replaces every korkut row. Curator and experimental rows are
untouched by the readout mode, so they are carried over from the old files verbatim. Refuses unless
every variant has both an old file and a re-run file, and the re-run has a korkut row for every old
korkut key (the same population, revalued).

Usage: merge_korkut_278.py OLD_RUNS_DIR KORKUT_RERUN_DIR NEW_RUNS_DIR
(writes NEW_RUNS_DIR/calls_<variant>.tsv, uncompressed, for score_157.py)
"""
import csv
import gzip
import sys
from pathlib import Path

HEADER = ["axis", "pathway", "case_a", "case_b", "readout", "fold", "call", "converged"]


def rows(path):
    op = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    with op as f:
        r = csv.reader(f, delimiter="\t")
        head = next(r)
        if head != HEADER:
            raise SystemExit(f"{path}: unexpected header {head}")
        return list(r)


def key(row):
    return tuple(row[:5])


def main():
    old_dir, rerun_dir, new_dir = (Path(p) for p in sys.argv[1:4])
    new_dir.mkdir(parents=True, exist_ok=True)
    olds = {p.name[len("calls_"):].split(".")[0]: p for p in old_dir.glob("calls_*.tsv*")}
    reruns = {p.name[len("calls_"):].split(".")[0]: p for p in rerun_dir.glob("calls_*.tsv*")}
    if set(olds) != set(reruns):
        raise SystemExit(f"variant sets differ: only old {sorted(set(olds) - set(reruns))}, "
                         f"only re-run {sorted(set(reruns) - set(olds))}")
    for v in sorted(olds):
        old = rows(olds[v])
        new_k = [r for r in rows(reruns[v]) if r[0] == "korkut"]
        old_k = {key(r) for r in old if r[0] == "korkut"}
        got = {key(r) for r in new_k}
        if old_k != got:
            raise SystemExit(f"{v}: korkut keys differ (old {len(old_k)}, re-run {len(got)}, "
                             f"missing {len(old_k - got)}, extra {len(got - old_k)})")
        kept = [r for r in old if r[0] != "korkut"]
        with open(new_dir / f"calls_{v}.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(HEADER)
            w.writerows(kept + new_k)
        print(f"{v}: {len(kept)} curator/experimental rows kept, {len(new_k)} korkut rows replaced")


if __name__ == "__main__":
    main()
