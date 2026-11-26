#!/usr/bin/env python3
"""specs/173 amendment 1: the Lun 2019 measured call from the authors' per-sample GFP bins, and the A/B split.

The deposited per-cell file is truncated (3 of 14 plate pairs), so the call is made from `bindat_perbin.csv`: per
sample, ten equal-width GFP bins (set within each sample), each readout's median and cell count per bin.

  split    stream the binned file once into an A file (controls + half-A proteins) and a B file (controls + half-B
           proteins); print row counts only
  score    the measured call on half A (refuses any half-B protein), then R1 and R2 against the frozen predictions

Usage: lun_a_173.py split BINNED OUT_A OUT_B
       lun_a_173.py score A_FILE PREDICTIONS OUT_JSON
"""
import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent))
from common_173 import split_of  # noqa: E402

MAPS = HERE.parents[1] / "127"
CONTROLS = ("FLAG-GFP", "untransfected")
BASAL = "0"
MIN_BIN_CELLS = 20      # specs/127 step 4
LOW_BINS, HIGH_BINS = (0, 1, 2), (7, 8, 9)
MIN_GROUP_BINS = 2
NULL_PCT = 95.0         # specs/127 step 5
MIN_CONTROLS = 10
BOOTSTRAP = 10_000
SEED = 173
# file target -> registered readout (bench/analysis/127/readouts.tsv)
TARGETS = {
    "b-catenin": "beta-catenin", "p-4EBP1": "p-4EBP1", "p-AKT": "p-AKT", "p-AMPKa": "p-AMPKa", "p-BTK": "p-BTK/ITK",
    "p-CREB": "p-CREB/ATF1", "p-ERK1-2": "p-ERK1/2", "p-FAK": "p-FAK", "p-GSK3B": "p-GSK3b", "p-JNK": "p-JNK",
    "p-MAPKAPK2": "p-MAPKAPK2", "p-MARCKS": "p-MARCKS", "p-MEK1-2": "p-MEK1/2", "p-MKK3-6": "p-MKK3/6",
    "p-MKK4": "p-MKK4", "p-NFkB": "p-NFkB", "p-PDK1": "p-PDK1", "p-PLCg2": "p-PLCg2", "p-RAF": "p-c-RAF",
    "p-S6": "p-S6", "p-SMAD2-3": "p-SMAD2/3", "p-SRC": "p-SRC", "p-STAT1": "p-STAT1", "p-STAT3": "p-STAT3",
    "p-STAT5": "p-STAT5", "p-mTOR": "p-mTOR", "p-p38": "p-p38", "p-p53": "p-p53", "p-p70S6K": "p-p70S6K",
    "p-p90RSK": "p-p90RSK",
}


def read_tsv(path):
    return list(csv.DictReader((ln for ln in open(path) if not ln.startswith("#")), delimiter="\t"))


def marker_genes(path=MAPS / "pois.tsv"):
    """Table S1 construct marker -> gene."""
    return {r["marker"]: r["gene"] for r in read_tsv(path)}


def is_control(marker):
    return marker.startswith(CONTROLS)


def split(binned, out_a, out_b, genes):
    """Stream once; controls go to both files, proteins to their half, anything unmapped to B. Returns counts."""
    counts = Counter()
    with open(binned) as src, open(out_a, "w") as fa, open(out_b, "w") as fb:
        header = src.readline()
        fa.write(header)
        fb.write(header)
        for line in src:
            marker = line.split(",", 2)[1]
            if is_control(marker):
                fa.write(line)
                fb.write(line)
                counts["control"] += 1
            elif marker in genes and split_of(genes[marker]) == "A":
                fa.write(line)
                counts["A"] += 1
            else:
                fb.write(line)
                counts["B" if marker in genes else "B_unmapped"] += 1
    return counts


def load_a(path, genes, timepoint=BASAL):
    """{(experiment, marker, row_col): {target: {bin: (counts, median)}}} for one timepoint; refuses half B."""
    out = defaultdict(lambda: defaultdict(dict))
    for r in csv.DictReader(open(path)):
        m = r["marker"]
        if not is_control(m) and (m not in genes or split_of(genes[m]) != "A"):
            raise SystemExit(f"a half-B or unmapped sample ({m}) reached the half-A reader: refusing (half B is sealed)")
        if r["timepoint"] != timepoint or r["target"] not in TARGETS:
            continue
        out[(r["experiment"], m, r["row_col"])][TARGETS[r["target"]]][int(r["bin"])] = (
            float(r["counts"]), float(r["median"]))
    return out


def dose_slope(bins):
    """Median readout in the high-GFP bins minus the low-GFP bins of one sample; None if too few cells."""
    def group(idx):
        v = [bins[i][1] for i in idx if i in bins and bins[i][0] >= MIN_BIN_CELLS]
        return float(np.median(v)) if len(v) >= MIN_GROUP_BINS else None
    hi, lo = group(HIGH_BINS), group(LOW_BINS)
    return None if hi is None or lo is None else hi - lo


def slopes(samples):
    """{sample: {readout: slope}} (None slopes dropped)."""
    out = {}
    for s, per in samples.items():
        out[s] = {r: v for r, v in ((r, dose_slope(b)) for r, b in per.items()) if v is not None}
    return out


def control_reference(sl, experiment, readout, exclude=None):
    """Median slope of FLAG-GFP samples on the same plate pair (all plates if none there), leaving one out."""
    def vals(same):
        return [v[readout] for s, v in sl.items() if s[1].startswith("FLAG-GFP") and s != exclude
                and readout in v and (not same or s[0] == experiment)]
    v = vals(True) or vals(False)
    return float(np.median(v)) if v else None


def thresholds(sl):
    """delta_r: 95th percentile (or the maximum, under 10 controls) of |E| over FLAG-GFP samples, leave one out."""
    effects = defaultdict(list)
    for s, per in sl.items():
        if not s[1].startswith("FLAG-GFP"):
            continue
        for r, v in per.items():
            ref = control_reference(sl, s[0], r, exclude=s)
            if ref is not None:
                effects[r].append(abs(v - ref))
    return {r: (float(np.percentile(e, NULL_PCT)) if len(e) >= MIN_CONTROLS else max(e)) for r, e in effects.items() if e}


def measured_calls(sl, genes, delta):
    """{(gene, readout): (call, E)}; a gene's E is the median over its samples (constructs, replicates)."""
    per = defaultdict(list)
    for s, rs in sl.items():
        if is_control(s[1]):
            continue
        for r, v in rs.items():
            ref = control_reference(sl, s[0], r)
            if ref is not None:
                per[(genes[s[1]], r)].append(v - ref)
    out = {}
    for k, es in per.items():
        e, d = float(np.median(es)), delta.get(k[1])
        if d is None:
            continue
        out[k] = ("up" if e > d else "down" if e < -d else "not", e)
    return out


def self_pairs(calls, readouts):
    rg = {r["readout"]: {x.split(":")[0] for x in r["sites"].split(",") if x} for r in readouts}
    return {k: v for k, v in calls.items() if k[0] in rg.get(k[1], set())}


def protein_bootstrap(cases, stat, n=BOOTSTRAP, seed=SEED):
    by = defaultdict(list)
    for c in cases:
        by[c["gene"]].append(c)
    genes = sorted(by)
    rng = random.Random(seed)
    draws = []
    for _ in range(n):
        v = stat([c for g in (rng.choice(genes) for _ in genes) for c in by[g]])
        if v is not None and not np.isnan(v):
            draws.append(v)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def enrichment(cases):
    reached = [c for c in cases if c["reach"] == "up"]
    unreached = [c for c in cases if c["reach"] != "up"]
    if not reached or not unreached:
        return None
    br = sum(c["measured"] != "not" for c in reached) / len(reached)
    bu = sum(c["measured"] != "not" for c in unreached) / len(unreached)
    return br / bu if bu else None


def lift(cases):
    reached = [c for c in cases if c["reach"] == "up"]
    calls = [c for c in reached if c["frozen"] in ("up", "down")]
    if not reached or not calls:
        return None
    return sum(c["measured"] != "not" for c in calls) / len(calls) - sum(c["measured"] != "not" for c in reached) / len(reached)


def main():
    cmd = sys.argv[1]
    genes = marker_genes()
    if cmd == "split":
        print(dict(split(sys.argv[2], sys.argv[3], sys.argv[4], genes)))
        return
    import dataset_roles
    import lun_score as L
    dataset_roles.require_role("lun2019_a", "fit")
    a_file, pred_file, out_json = sys.argv[2:5]
    sl = slopes(load_a(a_file, genes))
    delta = thresholds(sl)
    calls = measured_calls(sl, genes, delta)
    readouts = read_tsv(MAPS / "readouts.tsv")
    selfc = self_pairs(calls, readouts)
    passes, ups, nmeas = L.positive_control([c for c, _ in selfc.values()])
    res = {"positive_control": {"passes": passes, "up": ups, "measured": nmeas,
                                "pairs": {f"{g}|{r}": c for (g, r), (c, _) in sorted(selfc.items())}},
           "thresholds": delta, "measured_pairs": len(calls),
           "measured_classes": dict(Counter(c for c, _ in calls.values()))}
    if not passes:
        Path(out_json).write_text(json.dumps(res, indent=1))
        raise SystemExit("positive control FAILS: scoring stops (specs/127 step 7)")
    frozen, reach = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(pred_file), delimiter="\t"):
        if r["self"] == "1" or r["status"] != "scored" or split_of(r["gene"]) != "A":
            continue
        frozen[(r["gene"], r["readout"])].append(r["predicted"])
        reach[(r["gene"], r["readout"])].append(r["reach_predicted"])
    cases = [{"gene": k[0], "readout": k[1], "measured": calls[k][0], "E": calls[k][1],
              "frozen": L.combine_networks(frozen[k]), "reach": L.combine_networks(reach[k]), "never": "not"}
             for k in sorted(frozen) if k in calls]
    res["P_A_pairs"], res["P_A_measured"] = len(frozen), len(cases)
    res["P_A_measured_classes"] = dict(Counter(c["measured"] for c in cases))
    ba = {m: L.balanced_accuracy([(c["measured"], c[m]) for c in cases]) for m in ("frozen", "reach", "never")}
    res["R1"] = {"balanced_accuracy": ba}
    for a, b in (("frozen", "reach"), ("frozen", "never")):
        pt, lo, hi = L.bootstrap_difference(cases, a, b, seed=SEED)
        res["R1"][f"{a}-{b}"] = {"point": pt, "ci": [lo, hi], "verdict": L.verdict(lo, hi)}
    reached = [c for c in cases if c["reach"] == "up"]
    unreached = [c for c in cases if c["reach"] != "up"]
    fcalls = [c for c in reached if c["frozen"] in ("up", "down")]
    res["R2"] = {
        "reached": len(reached), "unreached": len(unreached),
        "base_rate_reached": sum(c["measured"] != "not" for c in reached) / len(reached),
        "base_rate_unreached": sum(c["measured"] != "not" for c in unreached) / len(unreached),
        "enrichment": enrichment(cases), "enrichment_ci": protein_bootstrap(cases, enrichment),
        "frozen_change_calls": len(fcalls),
        "precision_change": (sum(c["measured"] != "not" for c in fcalls) / len(fcalls)) if fcalls else None,
        "precision_exact": (sum(c["measured"] == c["frozen"] for c in fcalls) / len(fcalls)) if fcalls else None,
        "lift": lift(cases), "lift_ci": protein_bootstrap(cases, lift),
        "reached_measured_direction": dict(Counter(c["measured"] for c in reached)),
    }
    lib = {r["gene"]: r["library"] for r in read_tsv(MAPS / "pois.tsv")}
    res["by_library"] = {
        k: {"pairs": len(sub), "ba_frozen": L.balanced_accuracy([(c["measured"], c["frozen"]) for c in sub]),
            "ba_reach": L.balanced_accuracy([(c["measured"], c["reach"]) for c in sub])}
        for k in ("k", "p") for sub in [[c for c in cases if lib.get(c["gene"]) == k]] if sub}
    Path(out_json).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("positive_control", "measured_pairs", "P_A_pairs", "P_A_measured",
                                          "P_A_measured_classes", "R1", "R2", "by_library")}, indent=1, default=str))


if __name__ == "__main__":
    main()
