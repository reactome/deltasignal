#!/usr/bin/env python3
"""specs/137: rescore a korkut case file with the model's call made at the dataset's own threshold, as registered.

From the stored `pred_ui` (UI scale, 1 = baseline): M calls DOWN if log2(pred) <= -log2 t, UP if >= log2 t (the truth's
own `call()`, korkut_bench.py), with t = 1.5; C is the file as scored (0.85 / 1.15). Gate: score_134's comparison
(balanced accuracy M - C, 95% interval by resampling drug sets) on all cases and with RAF set aside. Reported beside it:
S (C's calls against korkut's `measured_secondary`, truth at 1.15), threshold-free AUROC of the predicted log fold, and
a leave-one-drug-set-out choice of t from GRID.

Usage: score_137.py CASES [T] > out.json
"""
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_134 as s134  # noqa: E402

RAF = "R-HSA-5673001"
GRID = (1.15, 1.25, 1.5, 2.0, 3.0)
FLOOR = 1e-6   # a predicted 0 (a knockout read through) gets log2 of this, below every cutoff


def log_fold(ui):
    return math.log2(max(float(ui), FLOOR))


def call_at(ui, t):
    lf, thr = log_fold(ui), math.log2(t)
    return "down" if lf <= -thr else "up" if lf >= thr else "not"


def rescored(cases, t):
    """{key: row with `predicted` re-derived at t}; unscored rows unchanged."""
    out = {}
    for k, r in cases.items():
        r = dict(r)
        if r["status"] == "scored":
            r["predicted"] = call_at(r["pred_ui"], t)
        out[k] = r
    return out


def auroc(pos, neg):
    """P(score of a positive > score of a negative), ties counting half (Mann-Whitney)."""
    if not pos or not neg:
        return None
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum, i = 0.0, 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        mid = (i + j + 1) / 2   # mean 1-based rank of the tie block
        rank_sum += mid * sum(1 for _, p in allv[i:j] if p)
        i = j
    return (rank_sum - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def threshold_free(keys, cases):
    lf = {k: log_fold(cases[k]["pred_ui"]) for k in keys}
    m = {k: cases[k]["measured"] for k in keys}
    return {"down_vs_rest": auroc([-lf[k] for k in keys if m[k] == "down"], [-lf[k] for k in keys if m[k] != "down"]),
            "up_vs_rest": auroc([lf[k] for k in keys if m[k] == "up"], [lf[k] for k in keys if m[k] != "up"])}


def balanced_at(keys, cases, t):
    return s134.metrics([(call_at(cases[k]["pred_ui"], t), cases[k]["measured"]) for k in keys])["balanced"]


def cv_choice(keys, cases, grid=GRID):
    """Leave one drug set out: choose t on the rest by balanced accuracy (ties to the smaller t), call the held-out
    drug set with it; pooled balanced accuracy of those calls against the same cases at grid[0]."""
    by = defaultdict(list)
    for k in keys:
        by[k[1]].append(k)
    chosen, pairs = Counter(), []
    for g, held in sorted(by.items()):
        rest = [k for k in keys if k[1] != g]
        t = max(grid, key=lambda t: (balanced_at(rest, cases, t), -t))
        chosen[t] += len(held)
        pairs += [(call_at(cases[k]["pred_ui"], t), cases[k]["measured"]) for k in held]
    cv = s134.metrics(pairs)["balanced"]
    base = balanced_at(keys, cases, grid[0])
    return {"chosen_t_by_cases": {str(t): n for t, n in sorted(chosen.items())}, "cv_balanced": cv,
            "balanced_at_" + str(grid[0]): base, "cv_gain": cv - base,
            "in_sample": {str(t): balanced_at(keys, cases, t) for t in grid}}


def secondary(keys, cases):
    return s134.metrics([(cases[k]["predicted"], cases[k]["measured_secondary"]) for k in keys])


def report(cases, t=1.5):
    keys = sorted(k for k, r in cases.items() if r["status"] == "scored" and r["measured"] in s134.CLASSES)
    arm = rescored(cases, t)
    out = {"t": t}
    for name, ks in (("all", keys), ("raf_set_aside", [k for k in keys if k[0] != RAF])):
        out[name] = {"gate": s134.compare(ks, cases, arm),
                     "S_secondary_truth": secondary([k for k in ks if cases[k].get("measured_secondary") in s134.CLASSES],
                                                    cases),
                     "threshold_free": threshold_free(ks, cases),
                     "cv": cv_choice(ks, cases)}
    return out


def main():
    cases = s134.load(sys.argv[1])
    json.dump(report(cases, float(sys.argv[2]) if len(sys.argv) > 2 else 1.5), sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
