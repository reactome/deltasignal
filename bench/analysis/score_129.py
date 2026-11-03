#!/usr/bin/env python3
"""specs/129: choose one activating curve from the grid by cross-validation on measured FIT data.

Each grid point is an arm (`cf_off` the identity, `cf_K<K>_h<h>` the others) run on the experimental and korkut axes.
The objective is the mean over the two datasets of three-class balanced accuracy. Perturbations (experimental: gene;
korkut: drug set) are split into 5 folds; each fold is predicted by the point chosen on the other four (ties: larger
elasticity, then larger K, so the identity wins a tie). The gain is the CV objective minus the identity's on the same
cases, with a 95% interval from resampling perturbations within each dataset. Curator balanced accuracy is reported for
every point, never used to choose.

Usage: score_129.py RESULTS_DIR TABLE_DIR > out.json
"""
import csv
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from holdout_report import TUNING_PATHWAYS  # noqa: E402

CLASSES = ("up", "not", "down")
SEED = 129
FOLDS = 5
BOOTSTRAP = 10000
RAF = "RAF_MAP_kinase_cascade"
CODE = {"0": "down", "1": "not", "2": "up"}
IDENTITY = "cf_off"


def elasticity(K, h):
    """Slope at rest, in log terms, of g(f) = f^h (1 + K^h) / (f^h + K^h)."""
    if math.isinf(K):
        return h
    return h * K ** h / (1 + K ** h)


def read_table(path):
    """(K, h) of a grid table; every activating class carries the same curve."""
    rows = [r for r in csv.DictReader((ln for ln in open(path) if not ln.startswith("#")), delimiter="\t")]
    curves = {(float(r["K"]), float(r["h"])) for r in rows}
    if len(curves) != 1:
        raise ValueError(f"{path}: grid tables hold one shared curve, found {sorted(curves)}")
    return curves.pop()


def load_experimental(path, set_aside_raf=True):
    """{case: (gene, measured, predicted)} over valid cases."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["valid"] != "1" or (set_aside_raf and r["pathway"] == RAF):
                continue
            out[(r["pathway"], r["gene"], r["direction"], r["key_output"])] = (
                r["gene"], CODE[r["expected"]], CODE[r["predicted"]])
    return out


def load_korkut(path):
    """{case: (drug set, measured, predicted)} over scored cases with a measured call."""
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] != "scored" or r["measured"] not in CLASSES:
                continue
            out[(r["pathway_id"], r["drug_set"], r["readout"])] = (r["drug_set"], r["measured"], r["predicted"] or "not")
    return out


def load_curator(path):
    """{"heldout"|"tuning": [(measured, predicted)]} over valid cases."""
    out = defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["valid"] == "1":
                split = "tuning" if r["pathway"] in TUNING_PATHWAYS else "heldout"
                out[split].append((CODE[r["expected"]], CODE[r["predicted"]]))
    return out


def balanced_accuracy(pairs):
    hit, tot = defaultdict(int), defaultdict(int)
    for m, p in pairs:
        tot[m] += 1
        hit[m] += int(m == p)
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def assign_folds(groups, k=FOLDS, seed=SEED):
    """Each group to a fold: sorted, shuffled with the seed, dealt round-robin."""
    order = sorted(set(groups))
    random.Random(seed).shuffle(order)
    return {g: i % k for i, g in enumerate(order)}


def objective(preds, truth, keep):
    """Mean over datasets of balanced accuracy on the cases `keep(dataset, group)` admits."""
    scores = []
    for ds, cases in truth.items():
        pairs = [(m, preds[ds][c]) for c, (g, m) in cases.items() if keep(ds, g)]
        scores.append(balanced_accuracy(pairs))
    return sum(scores) / len(scores)


def choose(points, preds, truth, keep):
    """Best grid point on the admitted cases; ties go to the larger elasticity, then the larger K."""
    def rank(name):
        K, h = points[name]
        return (round(objective(preds[name], truth, keep), 12), elasticity(K, h), K)
    return max(sorted(points), key=rank)


def cross_validate(points, preds, truth, k=FOLDS, seed=SEED):
    """CV predictions {dataset: {case: call}}, the point chosen for each fold, and the full-data choice."""
    folds = {ds: assign_folds([g for g, _ in cases.values()], k, seed) for ds, cases in truth.items()}
    cv = {ds: {} for ds in truth}
    chosen = []
    for f in range(k):
        name = choose(points, preds, truth, lambda ds, g: folds[ds][g] != f)
        chosen.append(name)
        for ds, cases in truth.items():
            for c, (g, _) in cases.items():
                if folds[ds][g] == f:
                    cv[ds][c] = preds[name][ds][c]
    full = choose(points, preds, truth, lambda ds, g: True)
    return cv, chosen, full


def _group_counts(cases, calls):
    """Per group: hits and totals by measured class, as arrays [groups x classes]."""
    groups = sorted({g for g, _ in cases.values()})
    idx = {g: i for i, g in enumerate(groups)}
    hit = np.zeros((len(groups), len(CLASSES)))
    tot = np.zeros((len(groups), len(CLASSES)))
    for c, (g, m) in cases.items():
        j = CLASSES.index(m)
        tot[idx[g], j] += 1
        hit[idx[g], j] += calls[c] == m
    return hit, tot


def _balacc_weighted(hit, tot, w):
    """Balanced accuracy with groups weighted by `w` [draws x groups]; classes absent in a draw are skipped."""
    h, t = w @ hit, w @ tot
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(t > 0, h / t, np.nan)
    return np.nanmean(rate, axis=1)


def bootstrap_gain(truth, model_a, model_b, n=BOOTSTRAP, seed=SEED):
    """Mean over datasets of (balanced accuracy of A − B), with a 95% interval from resampling groups within each
    dataset. `model_*`: {dataset: {case: call}}. Returns (point, low, high)."""
    rng = np.random.default_rng(seed)
    point, draws = 0.0, np.zeros(n)
    for ds, cases in truth.items():
        ha, tot = _group_counts(cases, model_a[ds])
        hb, _ = _group_counts(cases, model_b[ds])
        ones = np.ones((1, tot.shape[0]))
        point += (_balacc_weighted(ha, tot, ones) - _balacc_weighted(hb, tot, ones))[0]
        w = rng.multinomial(tot.shape[0], np.full(tot.shape[0], 1 / tot.shape[0]), size=n).astype(float)
        draws += _balacc_weighted(ha, tot, w) - _balacc_weighted(hb, tot, w)
    k = len(truth)
    low, high = np.percentile(draws / k, [2.5, 97.5])
    return point / k, float(low), float(high)


def main():
    results, tables = Path(sys.argv[1]), Path(sys.argv[2])
    points = {IDENTITY: (math.inf, 1.0)}
    points.update({p.stem: read_table(p) for p in sorted(tables.glob("cf_*.tsv"))})
    missing = [n for n in points if not (results / n / "korkut_cases.tsv").exists()]
    if missing:
        sys.exit(f"score_129: no results for {missing}")

    raw = {n: {"experimental": load_experimental(results / n / "experimental_cases.tsv"),
               "korkut": load_korkut(results / n / "korkut_cases.tsv")} for n in points}
    common = {ds: set.intersection(*(set(raw[n][ds]) for n in points)) for ds in ("experimental", "korkut")}
    truth = {ds: {c: raw[IDENTITY][ds][c][:2] for c in sorted(keys)} for ds, keys in common.items()}
    for n in points:
        for ds, keys in common.items():
            if any(raw[n][ds][c][:2] != truth[ds][c] for c in keys):
                sys.exit(f"score_129: {n} {ds} disagrees with {IDENTITY} on a case's perturbation or measured call")
    preds = {n: {ds: {c: raw[n][ds][c][2] for c in keys} for ds, keys in common.items()} for n in points}

    cv, chosen, full = cross_validate(points, preds, truth)
    gain = bootstrap_gain(truth, cv, preds[IDENTITY])
    full_gain = bootstrap_gain(truth, preds[full], preds[IDENTITY])

    grid = {}
    for n, (K, h) in sorted(points.items()):
        cur = load_curator(results / n / "curator_cases.tsv")
        raf = load_experimental(results / n / "experimental_cases.tsv", set_aside_raf=False)
        grid[n] = {
            "K": "inf" if math.isinf(K) else K, "h": h, "elasticity": round(elasticity(K, h), 4),
            "objective": objective(preds[n], truth, lambda ds, g: True),
            "experimental": balanced_accuracy([(m, preds[n]["experimental"][c])
                                               for c, (_, m) in truth["experimental"].items()]),
            "experimental_with_raf": balanced_accuracy([(m, p) for _, m, p in raf.values()]),
            "korkut": balanced_accuracy([(m, preds[n]["korkut"][c]) for c, (_, m) in truth["korkut"].items()]),
            "curator_heldout": balanced_accuracy(cur["heldout"]),
            "curator_tuning": balanced_accuracy(cur["tuning"]),
        }
    print(json.dumps({
        "cases": {ds: len(c) for ds, c in truth.items()},
        "perturbations": {ds: len({g for g, _ in c.values()}) for ds, c in truth.items()},
        "chosen_per_fold": chosen, "full_data_choice": full,
        "cv_gain": dict(zip(("point", "low", "high"), gain)),
        "full_choice_in_sample_gain": dict(zip(("point", "low", "high"), full_gain)),
        "earns_contraction": gain[1] > 0,
        "grid": grid,
    }, indent=1))


if __name__ == "__main__":
    main()
