#!/usr/bin/env python3
"""specs/184: score the HPN-DREAM ligand-stimulation predictions (bench/hpn_stim_bench.py) against the training data.

Usage: stim_184.py PREDICTIONS HPN_CSV_DIR OUT_JSON
  HPN_CSV_DIR holds <line>_main.csv (x_experimental/experimental/CSV of the challenge download).

The measured call, the scored pairs and the gate are registered in specs/184-ligand-stimulation/spec.md.
"""
import csv
import json
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
import lun_score as L  # noqa: E402

SEED, DRAWS = 184, 10_000
CELL_LINES = ("BT20", "BT549", "MCF7", "UACC812")
LIGANDS = ("EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1")
EARLY, LATE = ("5min", "15min", "30min", "60min"), ("2hr", "4hr")
FLOOR = math.log2(1.15)
RECEPTOR = {"EGF": ("EGFR_pY1068", "EGFR_pY1173", "EGFR_pY992"), "HGF": ("c-Met_pY1235",),
            "NRG1": ("HER3_pY1298", "HER2_pY1248")}
ARMS = ("model", "reach", "never")


def load(path):
    """{(stimulus, timepoint): [{readout: log2 value}]} for the DMSO rows of one line; the time-0 rows have stimulus ''."""
    rows = list(csv.reader(open(path)))
    names = next(r for r in rows if "Antibody Name" in r)
    first = names.index("Antibody Name") + 1
    head = next(i for i, r in enumerate(rows) if r[:4] == ["Cell Line", "Inhibitor", "Stimulus", "Timepoint"])
    out = defaultdict(list)
    for r in rows[head + 1:]:
        if len(r) < first or r[1] != "DMSO":
            continue
        out[(r[2], r[3])].append({n: math.log2(float(v)) for n, v in zip(names[first:], r[first:]) if v and float(v) > 0})
    return out


def mean_of(samples, readout):
    vals = [s[readout] for s in samples if readout in s]
    return sum(vals) / len(vals) if vals else None


def measured(data, ligand, window, mult=1.0):
    """{readout: (E, delta)}: E the mean over the window's time points of log2(ligand / PBS); delta = max(2 sd, floor)."""
    out = {}
    for r in sorted({n for s in data.get(("", "0min"), []) for n in s}):
        zero = [s[r] for s in data[("", "0min")] if r in s]
        ds = []
        for t in window:
            a, b = mean_of(data.get((ligand, t), []), r), mean_of(data.get(("PBS", t), []), r)
            if a is not None and b is not None:
                ds.append(a - b)
        if len(ds) >= 2 and len(zero) >= 2:
            out[r] = (sum(ds) / len(ds), mult * max(2 * statistics.stdev(zero), FLOOR))
    return out


def call_of(e, delta):
    return "up" if e > delta else "down" if e < -delta else "not"


def predictions(pred_file):
    """{(condition, readout): {"model", "reach", "withdraw"}} combined over networks as on Lun."""
    cols = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if r["status"] != "scored":
            continue
        k = (r["condition"], r["readout"])
        cols[k]["model"].append(r["predicted"])
        cols[k]["reach"].append(r["reach_predicted"])
        cols[k]["withdraw"].append(r["predicted_withdraw"])
    return {k: {a: L.as_scored(L.combine_networks(v)) for a, v in c.items()} for k, c in cols.items()}


def build_cases(all_data, pred, ligands, window, mult=1.0):
    cases, unpredicted = [], 0
    for cell, data in all_data.items():
        for lig in ligands:
            for r, (e, d) in measured(data, lig, window, mult).items():
                p = pred.get((lig, r))
                if p is None:
                    unpredicted += 1
                    continue
                cases.append({"cell": cell, "ligand": lig, "readout": r, "E": e, "delta": d, "measured": call_of(e, d),
                              "never": "not", **p})
    return cases, unpredicted


def changed(c):
    return c in ("up", "down")


def arm_stats(cases, arm):
    real = [c for c in cases if changed(c["measured"])]
    calls = [c for c in cases if changed(c[arm])]
    both = [c for c in real if changed(c[arm])]
    return {"balanced_accuracy": L.balanced_accuracy([(c["measured"], c[arm]) for c in cases]),
            "calls": len(calls), "right": sum(c[arm] == c["measured"] for c in calls),
            "precision": sum(c[arm] == c["measured"] for c in calls) / len(calls) if calls else None,
            "recall": sum(c[arm] == c["measured"] for c in real) / len(real) if real else None,
            "direction": sum(c[arm] == c["measured"] for c in both) / len(both) if both else None,
            "direction_n": len(both)}


def enrichment(cases):
    """(share of reachable pairs that change) / (share of unreachable pairs that change); None if a group is empty."""
    r = [changed(c["measured"]) for c in cases if c["reach"] == "up"]
    n = [changed(c["measured"]) for c in cases if c["reach"] != "up"]
    if not r or not n or not sum(n):
        return None
    return (sum(r) / len(r)) / (sum(n) / len(n))


def ba_diff(a, b):
    return lambda cs: (L.balanced_accuracy([(c["measured"], c[a]) for c in cs])
                       - L.balanced_accuracy([(c["measured"], c[b]) for c in cs]))


def bootstrap(cases, stat, key, n=DRAWS, seed=SEED):
    """Point and 95% percentile interval of `stat`, resampling the groups `key` defines with replacement."""
    by = defaultdict(list)
    for c in cases:
        by[key(c)].append(c)
    groups = sorted(by)
    rng = random.Random(seed)
    draws = []
    for _ in range(n):
        v = stat([c for g in (rng.choice(groups) for _ in groups) for c in by[g]])
        if v is not None and not math.isnan(v):
            draws.append(v)
    lo, hi = np.percentile(draws, [2.5, 97.5]) if draws else (float("nan"), float("nan"))
    return {"point": stat(cases), "ci": [float(lo), float(hi)], "groups": len(groups)}


def engagement(all_data, window):
    """{(cell, ligand): True/False/None}: the ligand's own receptor readout is measured UP (None: no such readout)."""
    out = {}
    for cell, data in all_data.items():
        for lig in LIGANDS:
            m = measured(data, lig, window)
            got = [call_of(*m[r]) == "up" for r in RECEPTOR.get(lig, ()) if r in m]
            out[(cell, lig)] = any(got) if got else None
    return out


def summary(cases):
    s = {"pairs": len(cases), "measured": dict(Counter(c["measured"] for c in cases)),
         "reachable": sum(c["reach"] == "up" for c in cases),
         "reached_base_rate": (sum(changed(c["measured"]) for c in cases if c["reach"] == "up")
                               / max(1, sum(c["reach"] == "up" for c in cases))),
         "unreached_rate": (sum(changed(c["measured"]) for c in cases if c["reach"] != "up")
                            / max(1, sum(c["reach"] != "up" for c in cases))),
         "enrichment": enrichment(cases), "arms": {a: arm_stats(cases, a) for a in ARMS}}
    return s


def analyse(all_data, pred):
    cond = lambda c: (c["cell"], c["ligand"])  # noqa: E731
    cases, unpred = build_cases(all_data, pred, LIGANDS, EARLY)
    res = {"primary": summary(cases), "pairs_without_prediction": unpred,
           "conditions": len({cond(c) for c in cases}),
           "G1_model-reach": {"by_condition": bootstrap(cases, ba_diff("model", "reach"), cond),
                              "by_ligand": bootstrap(cases, ba_diff("model", "reach"), lambda c: c["ligand"])},
           "model-never": bootstrap(cases, ba_diff("model", "never"), cond),
           "enrichment": {"by_condition": bootstrap(cases, enrichment, cond),
                          "by_ligand": bootstrap(cases, enrichment, lambda c: c["ligand"])}}
    res["G1_passes"] = res["G1_model-reach"]["by_condition"]["ci"][0] > 0
    res["by_ligand"] = {lig: summary([c for c in cases if c["ligand"] == lig]) for lig in LIGANDS}
    res["by_line"] = {cell: summary([c for c in cases if c["cell"] == cell]) for cell in all_data}
    eng = engagement(all_data, EARLY)
    res["engagement"] = {f"{k[0]} {k[1]}": v for k, v in sorted(eng.items())}
    for label, want in (("engaged", True), ("not_engaged", False), ("unknown", None)):
        sub = [c for c in cases if eng[cond(c)] is want and c["readout"] not in RECEPTOR.get(c["ligand"], ())]
        res[f"split_{label}"] = summary(sub) if sub else None
    late, _ = build_cases(all_data, pred, LIGANDS, LATE)
    strict, _ = build_cases(all_data, pred, LIGANDS, EARLY, mult=2.0)
    serum, _ = build_cases(all_data, pred, ("Serum",), EARLY)
    wd = [dict(c, model=c["withdraw"]) for c in cases]
    for label, cs in (("late_window", late), ("delta_doubled", strict), ("serum", serum), ("withdrawal", wd)):
        res[label] = {**summary(cs), "model-reach": bootstrap(cs, ba_diff("model", "reach"), cond)} if cs else None
    return res, cases


def main():
    import dataset_roles
    pred_file, csv_dir, out = sys.argv[1:4]
    dataset_roles.require_role("hpn_dream", "fit")
    all_data = {cell: load(Path(csv_dir) / f"{cell}_main.csv") for cell in CELL_LINES}
    res, cases = analyse(all_data, predictions(pred_file))
    Path(out).write_text(json.dumps(res, indent=1, default=str))
    with open(Path(out).with_suffix(".cases.tsv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(cases[0]), delimiter="\t")
        w.writeheader()
        w.writerows(cases)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
