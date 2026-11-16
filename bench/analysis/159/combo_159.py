#!/usr/bin/env python3
"""specs/159: the combination axis — the recorded pair solves vs registered
combination-naive nulls, entirely offline on the recorded korkut TSV.

Nulls (registered in spec.md US1, parameter-free):
- product: the pair's predicted fold is the product of the two singles'
  recorded predicted folds, clamped to the model's own [0, 100] UI bounds,
  called at the standard cutoffs.
- extreme: the single prediction with the larger |log2 fold| stands for the
  pair (a zero fold counts as infinitely extreme).

Gate (FR-003): recorded-pair accuracy vs the PRODUCT null's, 95% bootstrap
over pairs (10k resamples, seed 159), interval excluding zero. Balanced
accuracy, the extreme null, and the combination-specific class (the pair's
measured call differs from BOTH singles') reported whatever they say.

Usage: combo_159.py KORKUT_CASES_TSV
"""
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import dataset_roles  # noqa: E402

DOWN_CUT, UP_CUT = 0.85, 1.15
N_BOOT = 10_000
SEED = 159


def call(fold):
    return "down" if fold < DOWN_CUT else ("up" if fold >= UP_CUT else "not")


def extreme(fa, fb):
    ea = abs(math.log2(fa)) if fa > 0 else float("inf")
    eb = abs(math.log2(fb)) if fb > 0 else float("inf")
    return fa if ea >= eb else fb


def load(path):
    rows = [r for r in csv.DictReader(open(path), delimiter="\t")
            if r["status"] == "scored"]
    by_key = {(r["drug_set"], r["pathway_id"], r["readout"]): r for r in rows}
    pairs = sorted({r["drug_set"] for r in rows if "+" in r["drug_set"]})
    cases = []
    unjoinable = 0
    for p in pairs:
        a, b = p.split("+")
        for (ds, pw, ro), r in by_key.items():
            if ds != p:
                continue
            ra, rb = by_key.get((a, pw, ro)), by_key.get((b, pw, ro))
            if ra is None or rb is None:
                unjoinable += 1
                continue
            fp = float(r["pred_ui"])
            fa, fb = float(ra["pred_ui"]), float(rb["pred_ui"])
            cases.append({
                "pair": p, "pathway": pw, "readout": ro,
                "measured": r["measured"],
                "pair_call": r["predicted"],
                "product_call": call(min(100.0, max(0.0, fa * fb))),
                "extreme_call": call(extreme(fa, fb)),
                "combo_specific": r["measured"] != ra["measured"] and
                                  r["measured"] != rb["measured"],
                "pair_fold": fp,
            })
    return cases, unjoinable


def acc(cases, col):
    return sum(c[col] == c["measured"] for c in cases) / len(cases)


def ba(cases, col):
    per = []
    for m in ("down", "not", "up"):
        sub = [c for c in cases if c["measured"] == m]
        if sub:
            per.append(sum(c[col] == m for c in sub) / len(sub))
    return sum(per) / len(per)


def main():
    dataset_roles.require_role("korkut2015", "fit")
    cases, unjoinable = load(sys.argv[1])
    pairs = sorted({c["pair"] for c in cases})
    by_pair = defaultdict(list)
    for c in cases:
        by_pair[c["pair"]].append(c)
    print(f"{len(cases)} joinable pair cases over {len(pairs)} pairs "
          f"({unjoinable} unjoinable)")

    for col, name in (("pair_call", "recorded pair solve"),
                      ("product_call", "product null"),
                      ("extreme_call", "extreme null")):
        print(f"  {name:>19}: acc {acc(cases, col):.4f}  balanced {ba(cases, col):.4f}")

    rng = random.Random(SEED)
    deltas = []
    for _ in range(N_BOOT):
        picked = [by_pair[rng.choice(pairs)] for _ in pairs]
        flat = [c for g in picked for c in g]
        deltas.append(acc(flat, "pair_call") - acc(flat, "product_call"))
    deltas.sort()
    lo, hi = deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1]
    d = acc(cases, "pair_call") - acc(cases, "product_call")
    print(f"\nGATE (pair vs product null, accuracy): delta {d:+.4f} "
          f"[{lo:+.4f}, {hi:+.4f}] over {len(pairs)} pairs -> "
          f"{'PASS' if lo > 0 else 'FAIL'}")

    combo = [c for c in cases if c["combo_specific"]]
    if combo:
        print(f"\ncombination-specific cases (measured differs from BOTH "
              f"singles): {len(combo)}")
        for col, name in (("pair_call", "pair solve"),
                          ("product_call", "product null"),
                          ("extreme_call", "extreme null")):
            hit = sum(c[col] == c["measured"] for c in combo)
            print(f"  {name:>13}: {hit}/{len(combo)} right")
    agree = sum(c["pair_call"] == c["product_call"] for c in cases) / len(cases)
    print(f"\npair solve agrees with the product null on {agree:.1%} of cases "
          f"(where they differ, the solve's modelled interaction is the claim)")
    diff = [c for c in cases if c["pair_call"] != c["product_call"]]
    if diff:
        pr = sum(c["pair_call"] == c["measured"] for c in diff)
        nr = sum(c["product_call"] == c["measured"] for c in diff)
        print(f"  on those {len(diff)} disagreements: pair right {pr}, "
              f"product right {nr}")


if __name__ == "__main__":
    main()
