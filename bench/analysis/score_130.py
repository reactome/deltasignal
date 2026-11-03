#!/usr/bin/env python3
"""specs/130: the readout's resting level in the cell, chosen by cross-validation over CPPA training drugs.

A sample's predicted measured change for antibody a is m = log2(w 2^F + 1 - w), F the frozen predicted log2 fold and
w = 1 / (1 + 2^(q_theta(a) - p0)) the share of a's resting level above background in that cell: p0 the mean of a over
the cell line's control samples in the same set, q_theta(a) the theta-quantile of p0 over training cell lines. theta
is chosen per fold of drugs on the other folds (specs/121's per-drug correlation); the gain is the CV correlation
minus the frozen model's, bootstrapped over drugs. Training portion only; CPPA values never leave this machine.

Usage: score_130.py CPPA_UNPACKED SPLIT_TSV READOUTS_TSV CASES_TSV [CASES_TSV ...] > out.json
"""
import collections
import csv
import glob
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cppa_data_121  # noqa: E402
import score_121  # noqa: E402
import split_121  # noqa: E402

CONTROLS = {"CONTROL", "DMSO", "PBS", "MOCK", "VEHICLE"}
THETAS = (None, 0.05, 0.25, 0.5, 0.75)            # None: the identity (w = 1)
SEED = 130
FOLDS = 5
N_BOOT = 2000
THRESHOLD = math.log2(1.5)
CLASSES = (-1, 0, 1)


def control_baselines(unpacked, split_tsv):
    """{(set, cell line): {antibody: mean control value}} for training cell lines, from the p0_p1 tables. Rows of
    other cell lines, and every non-control row, are skipped before their values are parsed."""
    train_lines = {n for n, r in cppa_data_121.roles(split_tsv)["cell_line"].items() if r == "train"}
    controls = {}
    for f in sorted(glob.glob(f"{unpacked}/p0_p1_TreatmentSummary/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        for r in csv.DictReader(open(f, encoding="utf-8", errors="replace"), delimiter="\t"):
            line = split_121.norm(r["cell_line_name"])
            if (r["compound_name_1"] or "").strip().upper() in CONTROLS and line in train_lines:
                controls[(r["set"].strip(), r["UID"].strip())] = line
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for f in sorted(glob.glob(f"{unpacked}/p0_p1_RPPA/**/*.tsv", recursive=True)):
        if "__MACOSX" in f:
            continue
        set_id = "set" + Path(f).stem.rsplit("_set", 1)[1]
        with open(f, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                cells = line.rstrip("\n").split("\t")
                cl = controls.get((set_id, cells[0].strip()))
                if cl is None:
                    continue
                for a, v in zip(header[1:], cells[1:]):
                    try:
                        x = float(v)
                    except ValueError:
                        continue
                    if math.isfinite(x):
                        acc[(set_id, cl)][a].append(x)
    return {k: {a: sum(v) / len(v) for a, v in by.items()} for k, by in acc.items()}


def quantile(values, theta):
    """Linear-interpolated quantile of a non-empty list."""
    v = sorted(values)
    pos = theta * (len(v) - 1)
    lo = math.floor(pos)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def backgrounds(baselines, theta):
    """{antibody: q_theta} over every (set, cell line) baseline."""
    by = collections.defaultdict(list)
    for vals in baselines.values():
        for a, v in vals.items():
            by[a].append(v)
    return {a: quantile(v, theta) for a, v in by.items()}


def measured_change(fold, p0, q):
    """log2(w 2^F + 1 - w), w = 1 / (1 + 2^(q - p0)); F may be -inf (a knockout)."""
    w = 1.0 / (1.0 + 2.0 ** (q - p0))
    lin = w * (2.0 ** fold if fold > -math.inf else 0.0) + 1.0 - w
    return math.log2(lin) if lin > 0 else -score_121.CLIP


def sample_predictions(sample, fold_pred, baselines, bg):
    """{antibody: prediction} for one sample; bg None is the identity."""
    base = baselines.get((sample["set"], sample["cell_line"]), {})
    out = {}
    for a in sample["values"]:
        F = fold_pred.get((sample["drug"], a))
        if F is None:
            continue
        if bg is None or a not in base or a not in bg:
            out[a] = F
        else:
            out[a] = measured_change(F, base[a], bg[a])
    return out


def sample_scores(samples, fold_pred, readouts, baselines, bg):
    """[(drug, correlation, [(measured class, predicted class)])] per sample with at least 3 shared antibodies."""
    out = []
    for s in samples:
        p = sample_predictions(s, fold_pred, baselines, bg)
        pairs = [(p[a], v) for a, v in s["values"].items() if a in readouts and a in p]
        c = score_121.pearson([x for x, _ in pairs], [y for _, y in pairs])
        if c is not None:
            calls = [(cls(y), cls(x)) for x, y in pairs]
            out.append((s["drug"], c, calls))
    return out


def cls(x):
    return 1 if x >= THRESHOLD else -1 if x <= -THRESHOLD else 0


def per_drug(scores):
    by = collections.defaultdict(list)
    for d, c, _ in scores:
        by[d].append(c)
    return {d: sum(v) / len(v) for d, v in by.items()}


def balanced_accuracy(scores, drugs):
    hit, tot = collections.Counter(), collections.Counter()
    for d, _, calls in scores:
        if d in drugs:
            for m, p in calls:
                tot[m] += 1
                hit[m] += m == p
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def assign_folds(drugs, k=FOLDS, seed=SEED):
    order = sorted(drugs)
    random.Random(seed).shuffle(order)
    return {d: i % k for i, d in enumerate(order)}


def choose(per_drug_by_theta, drugs):
    """The theta with the highest mean per-drug correlation over `drugs`; ties to the identity, then smaller theta."""
    def rank(t):
        vals = [per_drug_by_theta[t][d] for d in drugs if d in per_drug_by_theta[t]]
        return (round(sum(vals) / len(vals), 12), t is None, -(t or 0))
    return max(per_drug_by_theta, key=rank)


def cross_validate(per_drug_by_theta, drugs, k=FOLDS, seed=SEED):
    """({drug: CV correlation}, [theta chosen per fold], full-data theta)."""
    folds = assign_folds(drugs, k, seed)
    cv, chosen = {}, []
    for f in range(k):
        t = choose(per_drug_by_theta, [d for d in drugs if folds[d] != f])
        chosen.append(t)
        for d in drugs:
            if folds[d] == f:
                cv[d] = per_drug_by_theta[t][d]
    return cv, chosen, choose(per_drug_by_theta, drugs)


def bootstrap_gain(a, b, drugs, n=N_BOOT, seed=SEED):
    """Mean of a − b over drugs and a 95% percentile interval from resampling drugs."""
    drugs = sorted(drugs)
    diff = {d: a[d] - b[d] for d in drugs}
    rng = random.Random(seed)
    draws = sorted(sum(diff[rng.choice(drugs)] for _ in drugs) / len(drugs) for _ in range(n))
    return sum(diff.values()) / len(drugs), draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


def main():
    unpacked, split_tsv, readouts_tsv, *logs = sys.argv[1:]
    readouts = {r["readout"] for r in csv.DictReader((ln for ln in open(readouts_tsv) if not ln.startswith("#")),
                                                     delimiter="\t") if r["kind"] != "exclude"}
    samples = cppa_data_121.load(unpacked, split_tsv, "train")
    fold_pred = score_121.predictions(logs)["fold"]
    baselines = control_baselines(unpacked, split_tsv)
    with_ctrl = sum((s["set"], s["cell_line"]) in baselines for s in samples)

    scores = {t: sample_scores(samples, fold_pred, readouts, baselines,
                               None if t is None else backgrounds(baselines, t)) for t in THETAS}
    pd = {t: per_drug(v) for t, v in scores.items()}
    drugs = sorted(set.intersection(*(set(v) for v in pd.values())))
    cv, chosen, full = cross_validate(pd, drugs)
    gain = bootstrap_gain(cv, pd[None], drugs)

    folds = assign_folds(drugs)
    cv_scores = [x for t, f in zip(chosen, range(FOLDS)) for x in scores[t] if folds.get(x[0]) == f]
    key = lambda t: "identity" if t is None else f"theta_{t}"  # noqa: E731
    print(json.dumps({
        "portion": "train", "samples": len(samples), "samples_with_controls": with_ctrl,
        "baselines": len(baselines), "drugs": len(drugs),
        "chosen_per_fold": [key(t) for t in chosen], "full_data_choice": key(full),
        "cv_gain": dict(zip(("point", "low", "high"), gain)), "earns_place": gain[1] > 0,
        "grid": {key(t): {"per_drug_correlation": sum(pd[t][d] for d in drugs) / len(drugs),
                          "balanced_accuracy": balanced_accuracy(scores[t], set(drugs))} for t in THETAS},
        "cv_balanced_accuracy": balanced_accuracy(cv_scores, set(drugs)),
    }, indent=1))


if __name__ == "__main__":
    main()
