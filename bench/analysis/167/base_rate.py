#!/usr/bin/env python3
"""specs/167: is a change call better than the base rate among the readouts a route reaches?

Read from recorded case files, no solve. Per axis:
  reached      the predicted fold is not 1 (a route carries the perturbation to the readout)
  base rate    how often the measured (or curated) label is a change, among reached and among unreached
  call rate    change calls / reached
  precision    measured change among change calls (`change`: any change; `exact`: the called direction)
  lift         precision(change) - base rate among reached, 95% bootstrap over perturbation units
  declined     reached readouts called no change, and how often those truly changed

Usage: base_rate.py FROZEN_RUN KORKUT_RUN DECRYPTM_CASES OUT_JSON
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))
sys.path.insert(0, str(HERE.parents[1] / "166"))

SEED = 167
RESAMPLES = 10_000


def load_decryptm(path):
    """{(pathway, drug, uniprot, position): (call, truth, fold)} for scored cases."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "scored":
            continue
        out[(r["pathway_id"], r["drug"], r["uniprot"], r["position"])] = (
            r["predicted"], r["measured"], float(r["pred_ui"]))
    return out


def unit_counts(cases, unit):
    """{unit: [reached, reached changes, calls, calls on a change, calls exact, unreached, unreached changes]}."""
    u = defaultdict(lambda: [0] * 7)
    for k, (call, truth, fold) in cases.items():
        c = u[unit(k)]
        chg = truth != "not"
        if fold == 1.0:
            c[5] += 1
            c[6] += chg
            continue
        c[0] += 1
        c[1] += chg
        if call != "not":
            c[2] += 1
            c[3] += chg
            c[4] += call == truth
    return u


def ratio(a, b):
    return a / b if b else None


def lift(t):
    """precision(change) - base rate among reached, from summed counts; None when either is undefined."""
    if not t[2] or not t[0]:
        return None
    return t[3] / t[2] - t[1] / t[0]


def stats(cases, unit, resamples=RESAMPLES, seed=SEED):
    u = unit_counts(cases, unit)
    m = np.array([u[k] for k in sorted(u)], dtype=float).reshape(-1, 7)
    t = m.sum(axis=0)
    rng = np.random.default_rng(seed)
    boot = []
    for _ in range(resamples if len(m) else 0):
        b = lift(m[rng.integers(0, len(m), len(m))].sum(axis=0))
        if b is not None:
            boot.append(b)
    lo, hi = (float(x) for x in np.percentile(boot, [2.5, 97.5])) if boot else (None, None)
    declined = t[0] - t[2]
    return {
        "cases": int(t[0] + t[5]), "units": len(m),
        "reached": int(t[0]), "unreached": int(t[5]),
        "base_rate_reached": ratio(t[1], t[0]), "base_rate_unreached": ratio(t[6], t[5]),
        "calls": int(t[2]), "call_rate": ratio(t[2], t[0]),
        "precision_change": ratio(t[3], t[2]), "precision_exact": ratio(t[4], t[2]),
        "sign_right": ratio(t[4], t[3]),
        "lift": lift(t), "lift_ci": [lo, hi],
        "declined": int(declined), "declined_change_rate": ratio(t[1] - t[3], declined),
    }


def main():
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    import score_166 as s
    frozen, korkut, decryptm, out_json = (Path(a) for a in sys.argv[1:5])
    axes = {}
    for lab in ("curator", "experimental"):
        axes.update(s.load_mpb(frozen / f"{lab}_cases.tsv", set(TUNING_PATHWAYS), lab))
    mpb_unit = lambda k: k[:3]  # noqa: E731  (pathway, gene, direction)
    res = {ax: stats(axes[ax], mpb_unit) for ax in ("curator_heldout", "curator_tuning", "experimental")}
    res["korkut"] = stats(s.load_korkut(korkut / "korkut_cases.tsv", marker_names()), lambda k: k[1])
    res["korkut_with_own_markers"] = stats(s.load_korkut(korkut / "korkut_cases.tsv", {}), lambda k: k[1])
    res["decryptm"] = stats(load_decryptm(decryptm), lambda k: k[1])
    out_json.write_text(json.dumps(res, indent=1))
    f = lambda x: "  n/a" if x is None else f"{x:.3f}"  # noqa: E731
    print("axis\tcases\treached\tbase(reached)\tbase(unreached)\tcall rate\tprecision\texact\tlift [95%]\t"
          "declined\tdeclined changed")
    for ax, r in res.items():
        print(f"{ax}\t{r['cases']}\t{r['reached']}\t{f(r['base_rate_reached'])}\t{f(r['base_rate_unreached'])}\t"
              f"{f(r['call_rate'])}\t{f(r['precision_change'])}\t{f(r['precision_exact'])}\t"
              f"{f(r['lift'])} [{f(r['lift_ci'][0])}, {f(r['lift_ci'][1])}]\t{r['declined']}\t"
              f"{f(r['declined_change_rate'])}")


if __name__ == "__main__":
    main()
