#!/usr/bin/env python3
"""specs/146: fit the per-call confidence table on the curator TUNING split and
validate its calibration on curator held-out and experimental cases.

confidence = the Laplace-smoothed share of same-stratum tuning calls that
matched the curators. Strata: call x railed x basis x magnitude-bin x
converged, with the registered >=50-case fallback hierarchy (drop basis, then
magnitude, then converged, then railed; root = per-call). Every reused
definition is imported from the record's own code (FR-001): TUNING_PATHWAYS
(holdout_report), railed (score_124), basis (the specs/141 dump).

Usage:
  fit_confidence.py RESULTS_DIR [--out src/api/confidence.json] [--seed 146]
RESULTS_DIR holds curator_cases.tsv and experimental_cases.tsv (the frozen
record: builds/20261026-1349_8b6527e/results/05dc867).
"""

import argparse
import csv
import json
import math
import random
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from holdout_report import TUNING_PATHWAYS  # noqa: E402
from score_124 import railed  # noqa: E402

CEILING = 100.0
MIN_CELL = 50
CALL = {"0": "down", "1": "not", "2": "up"}
LEVELS = ("full", "no_basis", "no_magnitude", "no_converged", "call_only")


def magnitude_bin(fold):
    if fold == 1.0:
        return "baseline"           # untouched: the no_path-shaped class
    a = abs(math.log2(fold)) if fold > 0 else float("inf")
    return "le1" if a <= 1 else ("le3" if a <= 3 else "gt3")


def stratum_key(level, call, is_railed, basis, mag, conv):
    if level == "full":
        return (call, is_railed, basis, mag, conv)
    if level == "no_basis":
        return (call, is_railed, mag, conv)
    if level == "no_magnitude":
        return (call, is_railed, conv)
    if level == "no_converged":
        return (call, is_railed)
    return (call,)


def load_basis(dump_path):
    d = json.loads(Path(dump_path).read_text())
    out = {}
    for pw in d["curator"]:
        for r in pw["readouts"]:
            out[(pw["pathway_name"], r["ko"])] = r["basis"]
    return out


def features(row, basis_by_key):
    fold = float(row["pred_ui"])
    return {
        "call": CALL[row["predicted"]],
        "railed": bool(railed(fold, CEILING)),
        "basis": basis_by_key.get((row["pathway"], row["key_output"])),
        "mag": magnitude_bin(fold),
        "conv": row["converged"] == "1",
        "correct": row["predicted"] == row["expected"],
        "pathway": row["pathway"],
    }


def cases(path, basis_by_key):
    return [features(r, basis_by_key) for r in csv.DictReader(open(path), delimiter="\t")
            if r["valid"] == "1"]


def fit(train):
    """Counts (n, correct) per stratum at every level of the hierarchy."""
    counts = {lv: defaultdict(lambda: [0, 0]) for lv in LEVELS}
    for c in train:
        for lv in LEVELS:
            if lv == "full" and c["basis"] is None:
                continue
            k = stratum_key(lv, c["call"], c["railed"], c["basis"], c["mag"], c["conv"])
            cell = counts[lv][k]
            cell[0] += 1
            cell[1] += c["correct"]
    return counts


def lookup(counts, c, converged_known=True):
    """(confidence, resolved level); the registered >=MIN_CELL walk."""
    for lv in LEVELS:
        if lv == "full" and c["basis"] is None:
            continue
        if not converged_known and lv in ("full", "no_basis", "no_magnitude"):
            continue  # those levels key on converged
        k = stratum_key(lv, c["call"], c["railed"], c["basis"], c["mag"], c["conv"])
        n, right = counts[lv].get(k, (0, 0))
        if n >= MIN_CELL:
            return (right + 1) / (n + 2), lv
    k = stratum_key("call_only", c["call"], None, None, None, None)
    n, right = counts["call_only"].get(k, (0, 0))
    return (right + 1) / (n + 2), "call_only"


def ece(scored):
    """Confidence-weighted calibration error over the resolved strata."""
    buckets = defaultdict(lambda: [0, 0])
    for conf, correct in scored:
        b = buckets[round(conf, 4)]
        b[0] += 1
        b[1] += correct
    total = sum(n for n, _ in buckets.values())
    return sum(n * abs(conf - right / n) for conf, (n, right) in buckets.items()) / total


def auc(scored):
    pos = sorted(c for c, ok in scored if ok)
    neg = sorted(c for c, ok in scored if not ok)
    if not pos or not neg:
        return None
    import bisect
    wins = 0.0
    for c in pos:
        lo = bisect.bisect_left(neg, c)
        hi = bisect.bisect_right(neg, c)
        wins += lo + 0.5 * (hi - lo)
    return wins / (len(pos) * len(neg))


def evaluate(name, split, counts, seed):
    scored = [(lookup(counts, c)[0], c["correct"]) for c in split]
    e, a = ece(scored), auc(scored)
    by_pw = defaultdict(list)
    for c, (conf, ok) in zip(split, scored):
        by_pw[c["pathway"]].append((conf, ok))
    pws = sorted(by_pw)
    rng = random.Random(seed)
    aucs = []
    for _ in range(2000):
        sample = [p for p in (rng.choice(pws) for _ in pws) for p in by_pw[p]]
        v = auc(sample)
        if v is not None:
            aucs.append(v)
    aucs.sort()
    lo, hi = aucs[int(0.025 * len(aucs))], aucs[int(0.975 * len(aucs)) - 1]
    acc = sum(ok for _, ok in scored) / len(scored)
    mean_conf = sum(c for c, _ in scored) / len(scored)
    print(f"{name}: n={len(scored)} accuracy {acc:.4f} mean confidence {mean_conf:.4f} "
          f"ECE {e:.4f} AUC {a:.4f} [{lo:.4f}, {hi:.4f}]")
    return {"n": len(scored), "accuracy": round(acc, 4), "mean_confidence": round(mean_conf, 4),
            "ece": round(e, 4), "auc": round(a, 4), "auc_95": [round(lo, 4), round(hi, 4)]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results_dir")
    ap.add_argument("--out", default=str(HERE.parents[3] / "src" / "api" / "confidence.json"))
    ap.add_argument("--seed", type=int, default=146)
    ap.add_argument("--basis", default=str(HERE.parents[1] / "141" / "bench_readouts.json"),
                    help="the specs/141-format dump whose basis matches the results' build")
    ap.add_argument("--catalog-build", default="20261026-1349_8b6527e",
                    help="provenance recorded in the table (must match the results dir's build)")
    ap.add_argument("--solver-commit", default="05dc867")
    args = ap.parse_args()
    rd = Path(args.results_dir)

    basis_by_key = load_basis(args.basis)
    curator = cases(rd / "curator_cases.tsv", basis_by_key)
    experimental = cases(rd / "experimental_cases.tsv", basis_by_key)
    train = [c for c in curator if c["pathway"] in TUNING_PATHWAYS]
    heldout = [c for c in curator if c["pathway"] not in TUNING_PATHWAYS]
    print(f"fit on tuning: {len(train)}; validate on held-out {len(heldout)} "
          f"+ experimental {len(experimental)}")
    no_basis = sum(c["basis"] is None for c in curator + experimental)
    print(f"cases without a joinable basis (start at no_basis level): {no_basis}")

    counts = fit(train)
    results = {
        "curator_held_out": evaluate("curator held-out", heldout, counts, args.seed),
        "experimental": evaluate("experimental", experimental, counts, args.seed + 1),
    }
    gate = (results["curator_held_out"]["ece"] <= 0.05
            and results["experimental"]["ece"] <= 0.10
            and results["curator_held_out"]["auc_95"][0] > 0.5
            and results["experimental"]["auc_95"][0] > 0.5)
    print("H1 " + ("PASSES" if gate else "FAILS"))

    exp_counts = fit(experimental)   # descriptive both-axes rates (amendment 1)
    table = {}
    for lv in LEVELS:
        table[lv] = {}
        for k, (n, r) in sorted(counts[lv].items()):
            en, er = exp_counts[lv].get(k, (0, 0))
            table[lv]["|".join(map(str, k))] = {
                "n": n, "correct": r, "confidence": round((r + 1) / (n + 2), 4),
                "experimental_n": en, "experimental_correct": er,
            }
    out = {
        "version": str(date.today()),
        "spec": "specs/146-call-confidence",
        "fit_split": "curator tuning (the 11 TUNING_PATHWAYS)",
        "catalog_build": args.catalog_build,
        "solver_commit": args.solver_commit,
        "min_cell": MIN_CELL,
        "levels": list(LEVELS),
        # amendment 1: an untouched readout (fold exactly 1.0) is served
        # confidence null, stratum "untouched" — research.md has the trace.
        "untouched_rule": "fold == 1.0 -> confidence null, stratum 'untouched'",
        "validation": results,
        "table": table,
    }
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"table -> {args.out}")


if __name__ == "__main__":
    main()
