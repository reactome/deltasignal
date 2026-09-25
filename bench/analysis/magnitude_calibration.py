#!/usr/bin/env python3
"""Does the SIZE of a predicted change carry case-level signal about whether
its DIRECTION is right? Pre-registered in specs/021.

The observed U-shape -- direction accuracy high for large predicted folds and
low for small ones -- could be composition rather than signal:

  * predicted class: extreme folds may be mostly one direction, and that
    direction may simply be easier;
  * pathway: some pathways saturate more AND score higher. That is the
    between-pathway confound that killed nine earlier levers in this project.

So the decisive statistic is a STRATIFIED AUC: within each (pathway, predicted
class) stratum, how well does strength rank correct cases above incorrect ones,
combined across strata weighted by their number of correct/incorrect pairs.
That is a stratified Mann-Whitney U. The null permutes strength WITHIN strata,
which preserves every stratum's composition and destroys only the case-level
association.

The naive falsification test -- replace each magnitude with its class mean and
watch calibration collapse -- is deliberately NOT used: it removes within-class
variation by construction, so it cannot fail.

Usage:
  python bench/analysis/magnitude_calibration.py --cases curator_cases.tsv
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from holdout_report import TUNING_PATHWAYS  # noqa: E402

BASELINE_UI = 1.0
ZERO_FLOOR = 1e-6
CHANGED = {"0", "2"}          # predicted DOWN or UP; NORMAL ("1") is excluded


def strength(pred_ui: float) -> float:
    """|log10 fold| from baseline, with a predicted zero floored so it ranks as
    the strongest possible knockdown rather than as undefined."""
    return abs(math.log10(max(pred_ui, ZERO_FLOOR) / BASELINE_UI))


def load(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in rows:
        if r["predicted"] not in CHANGED:
            continue
        try:
            ui = float(r["pred_ui"])
        except (KeyError, ValueError):
            continue
        out.append({"pathway": r["pathway"], "cls": r["predicted"],
                    "correct": r["predicted"] == r["expected"], "s": strength(ui)})
    return out


def mann_whitney_parts(s: np.ndarray, correct: np.ndarray) -> tuple[float, int]:
    """(U, pairs) for strength ranking correct above incorrect; ties count 1/2."""
    n_c = int(correct.sum())
    n_w = len(correct) - n_c
    if n_c == 0 or n_w == 0:
        return 0.0, 0
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    sorted_s = s[order]
    # average ranks over ties
    i = 0
    while i < len(sorted_s):
        j = i
        while j + 1 < len(sorted_s) and sorted_s[j + 1] == sorted_s[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    u = ranks[correct].sum() - n_c * (n_c + 1) / 2.0
    return float(u), n_c * n_w


def build_strata(rows: list[dict], key) -> list[tuple[np.ndarray, np.ndarray]]:
    groups: dict = defaultdict(list)
    for r in rows:
        groups[key(r)].append(r)
    strata = []
    for g in groups.values():
        s = np.array([r["s"] for r in g], dtype=float)
        c = np.array([r["correct"] for r in g], dtype=bool)
        if 0 < c.sum() < len(c):          # informative only with both outcomes
            strata.append((s, c))
    return strata


def stratified_auc(strata) -> tuple[float, int]:
    total_u = total_pairs = 0
    for s, c in strata:
        u, pairs = mann_whitney_parts(s, c)
        total_u += u
        total_pairs += pairs
    return (total_u / total_pairs if total_pairs else float("nan")), total_pairs


def permutation_p(strata, observed: float, n_perm: int, seed: int) -> float:
    """One-sided: how often a within-stratum shuffle does at least as well."""
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_perm):
        shuffled = [(rng.permutation(s), c) for s, c in strata]
        if stratified_auc(shuffled)[0] >= observed:
            hits += 1
    return (hits + 1) / (n_perm + 1)


def report(label: str, rows: list[dict], n_perm: int, seed: int) -> dict:
    out = {"label": label, "cases": len(rows)}
    for name, key in (("unstratified", lambda r: 0),
                      ("within predicted class", lambda r: r["cls"]),
                      ("within (pathway, class)", lambda r: (r["pathway"], r["cls"]))):
        strata = build_strata(rows, key)
        auc, pairs = stratified_auc(strata)
        entry = {"auc": auc, "pairs": pairs, "strata": len(strata)}
        if name == "within (pathway, class)":
            entry["p"] = permutation_p(strata, auc, n_perm, seed)
            entry["cases_in_informative_strata"] = int(sum(len(c) for _, c in strata))
        out[name] = entry
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260925)
    a = ap.parse_args()

    rows = load(a.cases)
    held = [r for r in rows if r["pathway"] not in TUNING_PATHWAYS]
    tune = [r for r in rows if r["pathway"] in TUNING_PATHWAYS]

    for label, subset in (("HELD-OUT (decides)", held), ("TUNING (reported only)", tune)):
        r = report(label, subset, a.permutations, a.seed)
        print(f"\n{label}: {r['cases']:,} changed predictions")
        for name in ("unstratified", "within predicted class", "within (pathway, class)"):
            e = r[name]
            extra = ""
            if "p" in e:
                extra = (f"  permutation p = {e['p']:.4f}  "
                         f"({e['cases_in_informative_strata']:,} cases in informative strata)")
            print(f"  {name:<26} AUC {e['auc']:.4f}  over {e['strata']:>3} strata, "
                  f"{e['pairs']:>10,} pairs{extra}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
