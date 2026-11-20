#!/usr/bin/env python3
"""2026-11-20 review, finding 1: do the engine's labels, read off REAL solved folds, match the record's?

Joins the specs/157 canonical replay (folds straight from /api/readouts, unrounded) with the recorded
case files (folds printed at 6 figures) and counts disagreements in the untouched and railed labels,
under the old exact rule and under the record-precision rule (`fold_rules.py`, the engine's
`record_fold`).

Usage: untouched_parity.py CALLS_CANONICAL[.gz] CURATOR_CASES EXPERIMENTAL_CASES KORKUT_CASES
"""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fold_rules import record_fold, untouched  # noqa: E402


def railed_exact(f):
    return f >= 99.0 or f <= 0.01


def railed_record(f):
    return railed_exact(record_fold(f))


def recorded(cur, exp, kork):
    out = {}
    for axis, p in (("curator", cur), ("experimental", exp)):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r.get("valid") == "1" and r["pred_ui"]:
                out[(axis, r["pathway"], r["gene"], r["direction"], r["key_output"])] = float(r["pred_ui"])
    for r in csv.DictReader(open(kork), delimiter="\t"):
        if r["status"] == "scored" and r["pred_ui"]:
            out[("korkut", r["pathway_id"], r["drug_set"], "", r["readout"])] = float(r["pred_ui"])
    return out


def compare(replay, rec):
    c = Counter()
    for k, f in replay.items():
        if k not in rec:
            continue
        g = rec[k]
        c["joined"] += 1
        c["record_untouched"] += g == 1.0
        c["untouched_mismatch_exact"] += (f == 1.0) != (g == 1.0)
        c["untouched_mismatch_record"] += untouched(f) != (g == 1.0)
        c["railed_mismatch_exact"] += railed_exact(f) != railed_exact(g)
        c["railed_mismatch_record"] += railed_record(f) != railed_exact(g)
    return dict(c)


def main():
    path = sys.argv[1]
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with op as f:
        replay = {(r["axis"], r["pathway"], r["case_a"], r["case_b"], r["readout"]): float(r["fold"])
                  for r in csv.DictReader(f, delimiter="\t") if r["fold"]}
    for k, v in compare(replay, recorded(*sys.argv[2:5])).items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
