#!/usr/bin/env python3
"""specs/121: score predictions against CPPA delta-p, the registered primary endpoint.

Per sample: Pearson correlation, over the antibodies that have both a prediction and a measurement, between the
prediction and delta-p (0 when either is constant, as specs/118). Averaged per drug, then over drugs; 95% bootstrap
over drugs (2,000 resamples, seed 121). Arms from a cppa_bench case log, each the mean over the pathways that score
the (drug, readout):
  fold   the predicted log2 fold, clipped at log2(100) (specs/118's ds_fold; the registered primary prediction)
  call   the predicted class as -1 / 0 / +1 (specs/118's ds, for continuity)
  reach  -1 where a route links a target to the readout, else 0 (reachability's direction)
Development reads the training portion only (cppa_data_121.load refuses a test portion without final=True).

Usage: score_121.py CPPA_UNPACKED SPLIT_TSV READOUTS_TSV CASES_TSV [CASES_TSV ...]
"""
import collections
import csv
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_data_121  # noqa: E402

CLIP = math.log2(100)
CALL = {"down": -1.0, "not": 0.0, "up": 1.0}


def predictions(case_logs):
    """{arm: {(drug, readout): value}}, each the mean over the pathways that score it."""
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for path in case_logs:
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["status"] != "scored":
                continue
            k = (r["drug"], r["readout"])
            ui = float(r["pred_ui"])
            acc["fold"][k].append(max(-CLIP, min(CLIP, math.log2(ui))) if ui > 0 else -CLIP)
            acc["call"][k].append(CALL[r["predicted"]])
            acc["reach"][k].append(-1.0 if r["reach_predicted"] == "down" else 0.0)
    return {a: {k: sum(v) / len(v) for k, v in d.items()} for a, d in acc.items()}


def pearson(x, y):
    n = len(x)
    if n < 3:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    if sx == 0 or sy == 0:
        return 0.0
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def per_drug(samples, preds, readouts):
    """{arm: {drug: mean per-sample correlation}}; samples with fewer than 3 shared antibodies are left out."""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for s in samples:
        for arm, p in preds.items():
            pairs = [(p[(s["drug"], a)], v) for a, v in s["values"].items()
                     if a in readouts and (s["drug"], a) in p]
            c = pearson([x for x, _ in pairs], [y for _, y in pairs])
            if c is not None:
                out[arm][s["drug"]].append(c)
    return {arm: {d: sum(v) / len(v) for d, v in by.items()} for arm, by in out.items()}


def summarize(drug_corr, seed=121, n_boot=2000):
    arms = sorted(drug_corr)
    drugs = sorted(set.intersection(*(set(drug_corr[a]) for a in arms)))
    out = {"drugs": len(drugs)}
    for a in arms:
        out[a] = sum(drug_corr[a][d] for d in drugs) / len(drugs)
    rng = random.Random(seed)
    boots = collections.defaultdict(list)
    for _ in range(n_boot):
        sample = [rng.choice(drugs) for _ in drugs]
        m = {a: sum(drug_corr[a][d] for d in sample) / len(sample) for a in arms}
        for a in arms:
            boots[a].append(m[a])
        if "fold" in m and "reach" in m:
            boots["fold-reach"].append(m["fold"] - m["reach"])
    for k, v in boots.items():
        v.sort()
        out[f"{k}_ci95"] = [v[int(0.025 * n_boot)], v[int(0.975 * n_boot) - 1]]
    out["positive_drugs"] = {a: sum(drug_corr[a][d] > 0 for d in drugs) for a in arms}
    return out


def main():
    unpacked, split_tsv, readouts_tsv, *logs = sys.argv[1:]
    readouts = {r["readout"] for r in csv.DictReader((l for l in open(readouts_tsv) if not l.startswith("#")),
                                                     delimiter="\t") if r["kind"] != "exclude"}
    samples = cppa_data_121.load(unpacked, split_tsv, "train")
    preds = predictions(logs)
    dc = per_drug(samples, preds, readouts)
    print(json.dumps({"portion": "train", "samples": len(samples), **summarize(dc),
                      "per_drug": {a: {d: round(v, 4) for d, v in sorted(by.items())} for a, by in dc.items()}},
                     indent=1))


if __name__ == "__main__":
    main()
