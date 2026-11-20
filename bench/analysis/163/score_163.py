#!/usr/bin/env python3
"""specs/163: rank, don't call — precision@k of the model's per-condition
ranking of readouts vs reachability and random, on RECORDED predictions.

A condition is a list of readouts, each (pred_l2, reach, hit_model, hit_reach):
pred_l2 the mean clipped log2 fold over the pathways that score it; reach True
if a route reaches it; hit_model whether it measurably changed in the MODEL's
predicted direction; hit_reach the same for reachability's direction (down).
Ranking by |pred_l2| (model) or reach (reachability); a tied block at the
cut-off contributes its mean outcome; k = min(K, n).

Usage:
  score_163.py korkut KORKUT_CASES_TSV
  score_163.py cppa CPPA_UNPACKED SPLIT_TSV MARKERS_TSV DRUG_TARGETS_TSV OUT_JSON CASES_TSV [...]
"""
import csv
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))
sys.path.insert(0, str(HERE.parents[1] / "162"))

CLIP = math.log2(100)
T = math.log2(1.5)
N_BOOT = 10_000
SEED = 163


def l2(ui):
    return max(-CLIP, min(CLIP, math.log2(ui))) if ui > 0 else -CLIP


def precision_at_k(scores, hits, k):
    """Mean hit over the top k by score; the tied block at the cut-off counts at its mean."""
    n = len(scores)
    if n == 0:
        return None
    k = min(k, n)
    order = sorted(range(n), key=lambda i: -scores[i])
    taken, total, i = 0, 0.0, 0
    while i < n and taken < k:
        j = i
        while j < n and scores[order[j]] == scores[order[i]]:
            j += 1
        block = [hits[order[x]] for x in range(i, j)]
        room = k - taken
        if len(block) <= room:
            total += sum(block); taken += len(block)
        else:
            total += room * (sum(block) / len(block)); taken = k
        i = j
    return total / k


def evaluate(conditions, k_list=(10, 20)):
    """conditions: {drug: [[(pred_l2, reach, hit_model, hit_reach), ...], ...]} ->
    per-drug-averaged precision@k for model / reach / random, plus bootstrap gate on k=10."""
    def per_drug(fn):
        out = {}
        for d, conds in conditions.items():
            vals = [v for v in (fn(c) for c in conds) if v is not None]
            if vals:
                out[d] = sum(vals) / len(vals)
        return out
    res = {}
    for k in k_list:
        m = per_drug(lambda c: precision_at_k([abs(r[0]) for r in c], [r[2] for r in c], k))
        rch = per_drug(lambda c: precision_at_k([1.0 if r[1] else 0.0 for r in c], [r[3] for r in c], k))
        rnd = per_drug(lambda c: (sum(r[2] for r in c) / len(c)) if c else None)
        drugs = sorted(set(m) & set(rch))
        res[k] = {"model": sum(m[d] for d in drugs) / len(drugs),
                  "reach": sum(rch[d] for d in drugs) / len(drugs),
                  "random": sum(rnd[d] for d in drugs) / len(drugs),
                  "drugs": len(drugs), "_m": m, "_r": rch}
    m, rch = res[10]["_m"], res[10]["_r"]
    drugs = sorted(set(m) & set(rch))
    rng = random.Random(SEED)
    deltas = []
    for _ in range(N_BOOT):
        pick = [rng.choice(drugs) for _ in drugs]
        deltas.append(sum(m[d] - rch[d] for d in pick) / len(pick))
    deltas.sort()
    lo, hi = deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1]
    for k in res:
        res[k].pop("_m"); res[k].pop("_r")
    res["gate"] = {"delta": res[10]["model"] - res[10]["reach"], "lo": lo, "hi": hi,
                   "pass": lo > 0}
    return res


def korkut(cases_tsv):
    import dataset_roles
    from engagement_161 import marker_names
    dataset_roles.require_role("korkut2015", "fit")
    marks = marker_names()
    acc = defaultdict(list); reach = defaultdict(bool); meas = {}
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored":
            continue
        ds = r["drug_set"]
        if any(r["readout"] in marks.get(d, set()) for d in ds.split("+")):
            continue
        key = (ds, r["readout"])
        acc[key].append(l2(float(r["pred_ui"])))
        reach[key] |= r["reach_predicted"] == "down"
        meas[key] = r["measured"]
    conds = defaultdict(list)
    for (ds, ro), vals in acc.items():
        p = sum(vals) / len(vals)
        pdir = "down" if p < 0 else ("up" if p > 0 else None)
        conds[ds].append((p, reach[(ds, ro)], meas[(ds, ro)] == pdir, meas[(ds, ro)] == "down"))
    # korkut: a condition IS a drug set; "per drug" averaging is per drug set.
    return {ds: [rows] for ds, rows in conds.items()}


def cppa(unpacked, split_tsv, markers_tsv, drug_targets, case_files):
    import cppa_data_121 as cd
    import dataset_roles
    from engagement_162 import markers as load_markers
    for ds in ("cppa_train", "cppa_test"):       # as the korkut path: every portion is read
        dataset_roles.require_role(ds, "fit")
    marks = load_markers(markers_tsv, drug_targets)
    acc = defaultdict(list); reach = defaultdict(bool)
    for path in case_files:
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["status"] != "scored":
                continue
            key = (r["drug"], r["readout"])
            acc[key].append(l2(float(r["pred_ui"])))
            reach[key] |= r["reach_predicted"] == "down"
    pred = {k: sum(v) / len(v) for k, v in acc.items()}
    by_drug = defaultdict(set)
    for d, ro in pred:
        by_drug[d].add(ro)
    conds = defaultdict(list)
    for portion in cd.PORTIONS:
        for s in cd.load(unpacked, split_tsv, portion, final=True):
            d = s["drug"]
            rows = []
            for ro in by_drug.get(d, ()):
                if ro in marks.get(d, set()) or ro not in s["values"]:
                    continue
                p, dp = pred[(d, ro)], s["values"][ro]
                m = "up" if dp >= T else ("down" if dp <= -T else "not")
                pdir = "down" if p < 0 else ("up" if p > 0 else None)
                rows.append((p, reach[(d, ro)], m == pdir, m == "down"))
            if rows:
                conds[d].append(rows)
    return conds


def main():
    if sys.argv[1] == "korkut":
        res = evaluate(korkut(sys.argv[2]))
        print(json.dumps(res, indent=1))
    else:
        unpacked, split_tsv, markers_tsv, drug_targets, out_json, *files = sys.argv[2:]
        res = evaluate(cppa(unpacked, split_tsv, markers_tsv, drug_targets, files))
        Path(out_json).write_text(json.dumps(res, indent=1))
        print(json.dumps(res, indent=1))
    print(f"\nGATE: {'PASS' if res['gate']['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
