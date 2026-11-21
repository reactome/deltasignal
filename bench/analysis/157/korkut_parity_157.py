#!/usr/bin/env python3
"""Issue #278 parity gate: does a replay's korkut leg reproduce the korkut record fold for fold?

Joins a replay's calls_*.tsv (ensemble_157.py output: folds straight from /api/readouts) with the
recorded korkut_cases.tsv (korkut_bench.py output, pred_ui printed `:.6g`; at baseline 0.01 the fold
over rest IS pred_ui) and requires, for EVERY scored record case: the replay has it, its fold at the
record's precision (`fold_rules.record_fold`, 6 significant figures) equals the recorded pred_ui, and
its call equals the recorded predicted call. Counts are printed per readout kind; any difference is
listed and the exit is non-zero.

Usage: korkut_parity_157.py CALLS_TSV[.gz] KORKUT_CASES_TSV [READOUTS_TSV]
"""
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from fold_rules import record_fold  # noqa: E402

READOUTS = HERE.parents[1] / "062" / "readouts.tsv"


def read_kinds(path):
    return {r["readout"]: r["kind"]
            for r in csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t")}


def read_record(path):
    """{(pathway, drug set, readout): (fold, call)} for the record's scored korkut cases."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] == "scored":
            out[(r["pathway_id"], r["drug_set"], r["readout"])] = (float(r["pred_ui"]), r["predicted"])
    return out


def read_replay(path):
    op = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    out = {}
    with op as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["axis"] == "korkut":
                out[(r["pathway"], r["case_a"], r["readout"])] = (
                    float(r["fold"]) if r["fold"] else None, r["call"])
    return out


def compare(replay, record, kinds):
    """-> (Counter of totals / equal per kind, list of differences)."""
    c, diffs = Counter(), []
    for k, (rfold, rcall) in sorted(record.items()):
        kind = kinds.get(k[2], "?")
        c[f"{kind}_total"] += 1
        got = replay.get(k)
        if got is None:
            diffs.append((k, kind, "missing from replay", rfold, None))
            continue
        fold, call = got
        if fold is None:
            diffs.append((k, kind, "no fold", rfold, None))
        elif record_fold(fold) != rfold:
            diffs.append((k, kind, "fold", rfold, fold))
        elif call != rcall:
            diffs.append((k, kind, "call", rcall, call))
        else:
            c[f"{kind}_equal"] += 1
    return c, diffs


def main():
    replay = read_replay(sys.argv[1])
    record = read_record(sys.argv[2])
    kinds = read_kinds(sys.argv[3] if len(sys.argv) > 3 else READOUTS)
    c, diffs = compare(replay, record, kinds)
    equal = sum(v for k, v in c.items() if k.endswith("_equal"))
    for kind in sorted({k.rsplit("_", 1)[0] for k in c}):
        print(f"{kind}: {c[kind + '_equal']} of {c[kind + '_total']} equal")
    print(f"korkut: {equal} of {len(record)} recorded folds and calls reproduced at print precision")
    extra = len(set(replay) - set(record))
    if extra:
        print(f"note: {extra} replay rows are not scored record cases")
    for d in diffs[:40]:
        print("DIFF", *d, sep="\t")
    if diffs:
        print(f"{len(diffs)} differences: the parity gate FAILS")
        sys.exit(1)
    print("no differences: the parity gate passes")


if __name__ == "__main__":
    main()
