#!/usr/bin/env python3
"""specs/118: the best published korkut predictor (linear regression, Long et al. 2024) with and without DeltaSignal.

LR: X ~ R.D over the 87 responses, fitted by least squares, the held-out drug's column 0. LR + F adds c.F, F a
per-condition, per-response feature from the case log (DeltaSignal's call, its clipped log2 fold, or reachability's
call, each averaged over the pathways that score the case), c one scalar (or one per response) fitted with R by exact
least squares. Leave one drug out; Long et al.'s mean per-drug correlation, and balanced accuracy at log2(1.5).

Usage: lr_118.py KORKUT_DIR KORKUT_CASES [--dose]
"""
import collections
import csv
import json
import math
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dataset_roles  # noqa: E402
from train_116 import balanced_accuracy  # noqa: E402

N_RESP, N_PROT = 87, 82
ACTIVITY = ("901", "AK", "HN", "NT", "P6", "PLX", "RO", "ST", "Tm", "ZS", "RY", "SR")   # pert.csv's last 12 columns
CALL = {"down": -1.0, "not": 0.0, "up": 1.0}
CLIP = math.log2(100)
THRESHOLD = math.log2(1.5)


def read_matrix(path):
    return np.array([[float(v) for v in line.split(",")] for line in open(path).read().splitlines() if line.strip()])


def load(kdir, dose=False):
    kdir = Path(kdir)
    names = [n for n in open(kdir / "node_Index.csv").read().splitlines() if n.strip()]
    X = read_matrix(kdir / "expr.csv")[:, :N_RESP]
    sets, doses = [], []
    for line in open(kdir / "expr_index.txt"):
        parts = [p.split("|") for p in line.split()[0].split(",")]
        sets.append(tuple(sorted(c for c, _ in parts)))
        doses.append({c: float(v) for c, v in parts})
    if dose:
        D = np.array([[d.get(c, 0.0) for c in ACTIVITY] for d in doses])
    else:
        D = read_matrix(kdir / "pert.csv")[:, -len(ACTIVITY):]
    for i, s in enumerate(sets):          # the column order must match the drugs given
        assert {ACTIVITY[k] for k in np.flatnonzero(D[i])} == set(s), (i, s)
    assert X.shape[0] == D.shape[0] == len(sets)
    return names[:N_RESP], X, D, sets


def features(cases_tsv, names, sets):
    """{name: n x 87 matrix}: DS, DS-fold and Reach, each the mean over the pathways scoring (drug set, readout)."""
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored":
            continue
        key = (r["drug_set"], r["readout"])
        acc["ds"][key].append(CALL[r["predicted"]])
        acc["reach"][key].append(CALL[r["reach_predicted"]])
        p = float(r["pred_ui"])
        acc["ds_fold"][key].append(max(-CLIP, min(CLIP, math.log2(p))) if p > 0 else -CLIP)
    out = {}
    for f, d in acc.items():
        M = np.zeros((len(sets), N_RESP))
        for i, s in enumerate(sets):
            for j, n in enumerate(names):
                v = d.get(("+".join(s), n))
                if v:
                    M[i, j] = sum(v) / len(v)
        out[f] = M
    return out


def fit(D, X, F=None, per_response=False):
    """(R, c) minimising |X - D R - c F|^2; c a scalar, or a vector over responses."""
    pinv = np.linalg.pinv(D)
    if F is None:
        return pinv @ X, 0.0
    P = np.eye(len(D)) - D @ pinv
    MF, MX = P @ F, P @ X
    if per_response:
        den = (MF * MF).sum(0)
        c = np.where(den > 0, (MF * MX).sum(0) / np.where(den > 0, den, 1), 0.0)
    else:
        den = (MF * MF).sum()
        c = (MF * MX).sum() / den if den > 0 else 0.0
    return pinv @ (X - c * F), c


def corr(a, b):
    a, b = a.ravel(), b.ravel()
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def calls(M):
    return np.where(M >= THRESHOLD, "up", np.where(M <= -THRESHOLD, "down", "not"))


def lodo(X, D, sets, feats):
    """{arm: {drug: (per-drug correlation, single-drug correlation, [(pred, measured)])}} and the fitted c's."""
    arms = {"lr": (None, False)}
    for f in feats:
        arms[f"lr+{f}"] = (f, False)
        arms[f"lr+{f}/resp"] = (f, True)
    res, cs = collections.defaultdict(dict), collections.defaultdict(dict)
    for j, drug in enumerate(ACTIVITY):
        test = np.array([drug in s for s in sets])
        single = np.array([s == (drug,) for s in sets])
        keep = [k for k in range(len(ACTIVITY)) if k != j]
        Dtr, Dte = D[~test][:, keep], D[test][:, keep]
        for arm, (f, per) in arms.items():
            R, c = fit(Dtr, X[~test], None if f is None else feats[f][~test], per)
            pred = Dte @ R + (0.0 if f is None else c * feats[f][test])
            sp = pred[single[test]]
            res[arm][drug] = (corr(X[test], pred), corr(X[single], sp) if single.any() else None,
                              list(zip(calls(pred[:, :N_PROT]).ravel(), calls(X[test][:, :N_PROT]).ravel())))
            cs[arm][drug] = c if np.isscalar(c) else float(np.median(c))
    return res, cs


def boot(res, a, b, idx=0, seed=118, n=2000):
    drugs = [d for d in ACTIVITY if res[a][d][idx] is not None and (b is None or res[b][d][idx] is not None)]
    diff = {d: res[a][d][idx] - (0 if b is None else res[b][d][idx]) for d in drugs}
    rng = random.Random(seed)
    v = sorted(sum(diff[rng.choice(drugs)] for _ in drugs) / len(drugs) for _ in range(n))
    return [v[int(0.025 * n)], v[int(0.975 * n) - 1]]


def summarize(res, cs):
    out = {}
    for arm, by in res.items():
        r = [by[d][0] for d in ACTIVITY]
        s = [by[d][1] for d in ACTIVITY if by[d][1] is not None]
        out[arm] = {"mean_corr": sum(r) / len(r), "single_mean_corr": sum(s) / len(s) if s else None,
                    "balanced_accuracy": balanced_accuracy([p for d in ACTIVITY for p in by[d][2]]),
                    "per_drug_corr": {d: round(by[d][0], 4) for d in ACTIVITY}, "c_by_fold": cs[arm]}
    for f in ("ds", "ds_fold", "reach"):
        out[f"lr+{f}_minus_lr"] = boot(res, f"lr+{f}", "lr")
    out["lr+ds_minus_lr+reach"] = boot(res, "lr+ds", "lr+reach")
    out["lr+ds_fold_minus_lr+reach"] = boot(res, "lr+ds_fold", "lr+reach")
    out["lr+ds_single_corr"] = boot(res, "lr+ds", None, idx=1)
    return out


def main():
    kdir, cases = sys.argv[1:3]
    dataset_roles.require_role("korkut2015", "fit")
    names, X, D, sets = load(kdir, dose="--dose" in sys.argv)
    res, cs = lodo(X, D, sets, features(cases, names, sets))
    print(json.dumps(summarize(res, cs), indent=1))


if __name__ == "__main__":
    main()
