#!/usr/bin/env python3
"""specs/121 amendment 6: the one scoring of the CPPA test portions, for every registered arm together.

The only code that reads a test portion (cppa_data_121.load with final=True). Arms: the frozen model's fold, the
expression-gated fold (a sample without a gated prediction keeps the frozen one) and reachability, from the prediction
case logs. Per portion A (unseen drug), B (unseen cell line), C (both), the registered endpoints: per-sample Pearson
over mapped antibodies (at least 3), averaged per drug, then over drugs, with a 95% bootstrap (2,000 resamples, seed
121) over drugs (A, C) or cell lines (B); three-class balanced accuracy at |delta-p| >= log2 1.5; and the primary on
samples treated 12 hours or less (amendment 3). Every number stays outside the repository (the CPPA terms).

Usage: score_121_final.py CPPA_UNPACKED SPLIT_TSV READOUTS_TSV GATED_DIR FROZEN_CASES_TSV [...] > out.json
       (GATED_DIR holds cppa_gated_{A,B,C}*.tsv from run_cppa_gated.sh)
"""
import collections
import csv
import glob
import json
import math
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_data_121  # noqa: E402
import score_121  # noqa: E402

SEED = 121
N_BOOT = 2000
THRESHOLD = math.log2(1.5)        # amendment 3
MAX_MINUTES = 12 * 60             # amendment 3: the early-time secondary
CLASSES = (-1, 0, 1)
ARMS = ("fold", "gated", "reach")


def parse_minutes(time_str):
    """A treatment time as minutes: '5min', '30 min', '4hr', '0.5hr', '24 hrs', '3d', '2 days'; None if unparsed."""
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*(min(?:ute)?s?|h(?:ou)?rs?|h|d(?:ay)?s?)\s*\.?\s*",
                     str(time_str or ""), re.IGNORECASE)
    if not m:
        return None
    v, unit = float(m.group(1)), m.group(2).lower()
    return v if unit.startswith("m") else v * 60 if unit.startswith("h") else v * 24 * 60


def arm_predictions(frozen_logs, gated_logs):
    """{arm: prediction lookup}: `fold` and `reach` keyed (drug, readout); `gated` a pair (by-cell lookup keyed
    (drug, cell line, readout), fallback fold)."""
    frozen = score_121.predictions(frozen_logs)
    gated = score_121.gated_predictions(gated_logs)["fold"] if gated_logs else {}
    return {"fold": frozen["fold"], "reach": frozen["reach"], "gated": (gated, frozen["fold"])}


def sample_pairs(sample, preds, arm, readouts):
    """[(prediction, delta-p)] for one sample under one arm."""
    out = []
    for a, v in sample["values"].items():
        if a not in readouts:
            continue
        if arm == "gated":
            by_cell, fallback = preds["gated"]
            p = by_cell.get((sample["drug"], sample["cell_line"], a), fallback.get((sample["drug"], a)))
        else:
            p = preds[arm].get((sample["drug"], a))
        if p is not None:
            out.append((p, v))
    return out


def cls(x):
    return 1 if x >= THRESHOLD else -1 if x <= -THRESHOLD else 0


def collect(samples, preds, readouts):
    """Per-sample records: drug, cell line, minutes, and per arm the correlation and the (measured, predicted)
    classes. A sample scores only where every arm has at least 3 shared antibodies, so the arms share denominators."""
    out = []
    for s in samples:
        by_arm = {arm: sample_pairs(s, preds, arm, readouts) for arm in ARMS}
        corr = {arm: score_121.pearson([x for x, _ in p], [y for _, y in p]) for arm, p in by_arm.items()}
        if any(c is None for c in corr.values()):
            continue
        out.append({"drug": s["drug"], "cell_line": s["cell_line"], "minutes": parse_minutes(s.get("time")),
                    "corr": corr, "calls": {arm: [(cls(y), cls(x)) for x, y in p] for arm, p in by_arm.items()}})
    return out


def per_drug(records, arm, keep=lambda r: True):
    by = collections.defaultdict(list)
    for r in records:
        if keep(r):
            by[r["drug"]].append(r["corr"][arm])
    return {d: sum(v) / len(v) for d, v in by.items()}


def primary(records, arm, unit):
    """(mean per-drug correlation, 95% bootstrap interval) resampling drugs or cell lines (the portion's unit)."""
    point_by = per_drug(records, arm)
    drugs = sorted(point_by)
    point = sum(point_by.values()) / len(drugs)
    rng = random.Random(SEED)
    units = sorted({r[unit] for r in records})
    by_unit = collections.defaultdict(list)
    for r in records:
        by_unit[r[unit]].append(r)
    draws = []
    for _ in range(N_BOOT):
        rs = [r for _ in units for r in by_unit[rng.choice(units)]]
        pd = per_drug(rs, arm)
        draws.append(sum(pd.values()) / len(pd))
    draws.sort()
    return point, [draws[int(0.025 * N_BOOT)], draws[int(0.975 * N_BOOT) - 1]], len(drugs), len(units)


def balanced(records, arm):
    tot, hit = collections.Counter(), collections.Counter()
    for r in records:
        for m, p in r["calls"][arm]:
            tot[m] += 1
            hit[m] += m == p
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def portion_report(records, unit):
    out = {"samples_scored": len(records), "unit": unit}
    for arm in ARMS:
        point, ci, n_drugs, n_units = primary(records, arm, unit)
        early = per_drug(records, arm, keep=lambda r: r["minutes"] is not None and r["minutes"] <= MAX_MINUTES)
        out[arm] = {"per_drug_corr": point, "ci95": ci, "drugs": n_drugs, "units": n_units,
                    "balanced": balanced(records, arm),
                    "early_per_drug_corr": (sum(early.values()) / len(early)) if early else None,
                    "early_drugs": len(early)}
    # the paired difference the registration asks about (fold - reach), on the same samples
    diff = per_drug(records, "fold")
    reach = per_drug(records, "reach")
    out["fold_minus_reach"] = sum(diff[d] - reach[d] for d in diff) / len(diff)
    return out


def main():
    unpacked, split_tsv, readouts_tsv, gated_dir, *frozen_logs = sys.argv[1:]
    readouts = {r["readout"] for r in csv.DictReader((ln for ln in open(readouts_tsv) if not ln.startswith("#")),
                                                     delimiter="\t") if r["kind"] != "exclude"}
    report = {}
    for portion, unit in (("A", "drug"), ("B", "cell_line"), ("C", "drug")):
        gated_logs = sorted(glob.glob(f"{gated_dir}/cppa_gated_{portion}*.tsv"))
        preds = arm_predictions(frozen_logs, gated_logs)
        samples = cppa_data_121.load(unpacked, split_tsv, portion, final=True)
        records = collect(samples, preds, readouts)
        report[portion] = dict(portion_report(records, unit), samples=len(samples), gated_logs=len(gated_logs))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
