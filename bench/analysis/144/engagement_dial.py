#!/usr/bin/env python3
"""specs/144: the engagement dial, measured offline on the frozen korkut record.

For each drug set, one alpha fitted by least squares in log-fold space on the
set's own direct markers (bench/analysis/138/engagement.py MARKERS — imported,
never re-listed); conditioned fold = fold^alpha, classified at the bench's
cutoffs. Markers are excluded from scoring in BOTH arms. The gate (H1):
conditioned balanced accuracy over non-marker scored readouts beats the
control, 95% by-drug-set bootstrap on the difference excluding 0, at the
registered primary threshold.

Usage:
  engagement_dial.py CASES_TSV [--seed 144]
"""

import argparse
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "138"))
import korkut_bench as k  # noqa: E402
from engagement import MARKERS, table as engagement_table  # noqa: E402  the specs/138 code itself

DOWN, NORMAL, UP = "down", "not", "up"
EPS = 1e-8


def classify(fold):
    return DOWN if fold < k.bm.DOWN_CUTOFF else (UP if fold >= k.bm.UP_CUTOFF else NORMAL)


def measured_marker_folds():
    """Per drug code, per marker: mean log2 over that drug's own single-drug
    condition rows (both doses), via the bench's own readers."""
    names = [l.strip() for l in open(k.DATA / "node_Index.csv")]
    per_drug = defaultdict(list)
    for (i, ds), row in zip(k.conditions(), csv.reader(open(k.DATA / "expr.csv"))):
        per_drug[ds].append(dict(zip(names, map(float, row))))
    return per_drug


def set_rows(per_cond, drug_set):
    return per_cond.get(tuple(drug_set), [])


def markers_of(ds):
    """The drug set's markers: the union over its member drugs."""
    return sorted({m for d in ds for m, _ in MARKERS.get(d, [])})


def fit_alpha(pairs):
    """Least squares in natural-log space, clipped to [0, 1]; None when the
    predicted marker signal is below EPS (the registered undefined case)."""
    num = sum(m * f for m, f in pairs)
    den = sum(f * f for _, f in pairs)
    if den < EPS:
        return None
    return min(1.0, max(0.0, num / den))


def marker_predictions(rows):
    """(ds, pathway) -> marker -> predicted fold, from a korkut cases log."""
    pred = defaultdict(dict)
    for r in rows:
        ds = tuple(r["drug_set"].split("+"))
        if r["readout"] in markers_of(ds):
            pred[(ds, r["pathway_id"])][r["readout"]] = float(r["pred_ui"])
    return pred


def pooled_alphas(rows, per_cond):
    """The registered pooled estimator: per drug set, alpha over every
    pathway's scored markers (measured log2 and predicted fold both taken to
    natural logs)."""
    def measured_log2(ds, marker):
        vals = [v[marker] for v in set_rows(per_cond, ds) if marker in v]
        return sum(vals) / len(vals) if vals else None
    pred = marker_predictions(rows)
    pooled_pairs = defaultdict(list)
    for (ds, pw), marks in pred.items():
        for m, f in marks.items():
            m2 = measured_log2(ds, m)
            if m2 is not None and f > 0:
                pooled_pairs[ds].append((m2 * math.log(2), math.log(f)))
    return {ds: fit_alpha(pairs) for ds, pairs in pooled_pairs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--seed", type=int, default=144)
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.cases), delimiter="\t")
            if r["status"] == "scored" and r["measured"] != ""]
    per_cond = measured_marker_folds()

    def measured_log2(ds, marker):
        vals = [v[marker] for v in set_rows(per_cond, ds) if marker in v]
        return sum(vals) / len(vals) if vals else None

    pred = marker_predictions(rows)
    fit = fit_alpha

    # alpha per (ds, pathway), with the registered fallback chain
    fallbacks = {"pathway": 0, "pooled": 0, "one": 0}
    pooled_alpha = pooled_alphas(rows, per_cond)

    def alpha_for(ds, pw):
        marks = pred.get((ds, pw), {})
        pairs = []
        for m, f in marks.items():
            m2 = measured_log2(ds, m)
            if m2 is not None and f > 0:
                pairs.append((m2 * math.log(2), math.log(f)))
        a = fit(pairs) if pairs else None
        if a is not None:
            fallbacks["pathway"] += 1
            return a
        a = pooled_alpha.get(ds)
        if a is not None:
            fallbacks["pooled"] += 1
            return a
        fallbacks["one"] += 1
        return 1.0

    # score non-marker readouts in both arms, identical denominators
    cases = []
    for r in rows:
        ds = tuple(r["drug_set"].split("+"))
        if r["readout"] in markers_of(ds):
            continue
        fold = float(r["pred_ui"])
        a = alpha_for(ds, r["pathway_id"])
        cases.append({
            "ds": ds, "pathway": r["pathway_id"], "readout": r["readout"],
            "measured": r["measured"], "measured_secondary": r["measured_secondary"],
            "control": r["predicted"], "alpha": a,
            "arm": classify(fold ** a) if fold > 0 else classify(0.0),
            "fold": fold,
        })

    def balanced(cs, col, truth):
        per = defaultdict(lambda: [0, 0])
        for c in cs:
            per[c[truth]][0] += c[col] == c[truth]
            per[c[truth]][1] += 1
        recalls = [a / b for a, b in per.values() if b]
        return sum(recalls) / len(recalls)

    def raw(cs, col, truth):
        return sum(c[col] == c[truth] for c in cs) / len(cs)

    def always_not(cs, truth):
        return (balanced([dict(c, const=NORMAL) for c in cs], "const", truth),
                raw([dict(c, const=NORMAL) for c in cs], "const", truth))

    def report(label, cs):
        for truth, tname in (("measured", "primary |log2|>=log2 1.5"),
                             ("measured_secondary", "secondary |log2|>=log2 1.15")):
            bc, ba = balanced(cs, "control", truth), balanced(cs, "arm", truth)
            rc, ra = raw(cs, "control", truth), raw(cs, "arm", truth)
            nb, nr = always_not(cs, truth)
            print(f"{label} [{tname}] n={len(cs)}: balanced control {bc:.4f} arm {ba:.4f} "
                  f"(diff {ba-bc:+.4f}); raw control {rc:.4f} arm {ra:.4f}; always-NOT {nb:.4f}/{nr:.4f}")

    sets = sorted({c["ds"] for c in cases})
    by_set = {ds: [c for c in cases if c["ds"] == ds] for ds in sets}

    report("ALL", cases)
    noraf = [c for c in cases if c["pathway"] != "R-HSA-5673001"]
    report("RAF SET ASIDE", noraf)

    rng = random.Random(args.seed)
    diffs = []
    for _ in range(10000):
        sample = [by_set[rng.choice(sets)] for _ in sets]
        cs = [c for chunk in sample for c in chunk]
        diffs.append(balanced(cs, "arm", "measured") - balanced(cs, "control", "measured"))
    diffs.sort()
    lo, hi = diffs[249], diffs[9749]
    print(f"bootstrap 95% (drug sets, primary, balanced): [{lo:+.4f}, {hi:+.4f}]"
          f" -> gate {'PASSES' if lo > 0 else 'FAILS'}")

    # McNemar-style fixed/broken and the non-engaged drugs' trace
    fixed = [c for c in cases if c["arm"] == c["measured"] != c["control"]]
    broke = [c for c in cases if c["control"] == c["measured"] != c["arm"]]
    print(f"fixed {len(fixed)} / broke {len(broke)}; fallbacks {dict(fallbacks)}")
    nonengaged = {r["drug"] for r in engagement_table() if r["engaged"] != "yes"}
    print("non/unknown-engagement drugs (specs/138 table):", sorted(nonengaged))
    for group, name in ((fixed, "fixed"), (broke, "broke")):
        n = sum(1 for c in group if set(c["ds"]) & nonengaged)
        print(f"  {name} touching a non/unknown-engagement drug: {n} of {len(group)}")
    per_alpha = sorted({(c['ds'], round(c['alpha'], 3)) for c in cases})
    singles = [(ds, a) for ds, a in per_alpha if len(ds) == 1]
    print("single-drug alphas:", {"+".join(ds): a for ds, a in singles})


if __name__ == "__main__":
    main()
