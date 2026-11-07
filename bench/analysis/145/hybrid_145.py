#!/usr/bin/env python3
"""specs/145: DeltaSignal as the unseen drug's response row, in the specs/118 LODO.

The harness's leave-one-drug-out zeroes the held-out drug's design column, so
the fitted LR predicts zero direct effect for it. Arms (everything imported
from lr_118 / the specs/144 module, never re-implemented):

  lr           the control, exactly specs/118
  lr+ds_fold   specs/118's single-feature arm, reproduced for continuity
  lr+dsrow     + D[:, j] x (c * dsrow_j): DS's frozen single-drug clipped
               log2 response row for the held-out drug, c fitted on the
               TRAINING drugs' fitted rows only
  lr+reachrow  the same construction from reachability's calls (H1' control)
  lr+dialrow   dsrow_j scaled by the specs/144 pooled alpha for drug j
               (markers excluded from scoring in the H2 comparison, both arms)

Usage: hybrid_145.py KORKUT_DIR KORKUT_CASES
"""

import collections
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))            # bench/
sys.path.insert(0, str(HERE.parents[1]))            # bench/analysis/
sys.path.insert(0, str(HERE.parents[1] / "144"))    # the specs/144 module
import dataset_roles  # noqa: E402
import lr_118 as h  # noqa: E402  the specs/118 harness itself
from engagement_dial import measured_marker_folds, pooled_alphas, markers_of  # noqa: E402


def single_rows(sets):
    """drug -> the index of one of its single-drug conditions."""
    out = {}
    for i, s in enumerate(sets):
        if len(s) == 1 and s[0] not in out:
            out[s[0]] = i
    return out


def row_scale(R, keep, rows_by_drug):
    """c = argmin over the TRAINING drugs of |R_k - c * row_k|^2."""
    num = den = 0.0
    for pos, k in enumerate(keep):
        drug = h.ACTIVITY[k]
        row = rows_by_drug.get(drug)
        if row is None:
            continue
        num += float(R[pos] @ row)
        den += float(row @ row)
    return num / den if den > 0 else 0.0


def main():
    kdir, cases = sys.argv[1:3]
    dataset_roles.require_role("korkut2015", "fit")
    names, X, D, sets = h.load(kdir)
    feats = h.features(cases, names, sets)
    singles = single_rows(sets)
    dsrow = {d: feats["ds_fold"][i] for d, i in singles.items()}
    # reach calls are -1/0; as a log2-fold row a reach DOWN is a nominal
    # just-called fold, the natural magnitude-free translation of a call
    reachrow = {d: feats["reach"][i] * h.THRESHOLD * 1.01 for d, i in singles.items()}

    # specs/144 pooled alphas on the same frozen cases
    rows144 = [r for r in csv.DictReader(open(cases), delimiter="\t")
               if r["status"] == "scored" and r["measured"] != ""]
    alphas = pooled_alphas(rows144, measured_marker_folds())
    alpha1 = {d: (alphas.get((d,)) if alphas.get((d,)) is not None else 1.0) for d in h.ACTIVITY}

    res = collections.defaultdict(dict)
    cs = collections.defaultdict(dict)
    marker_cols = {d: [j for j, n in enumerate(names) if n in markers_of((d,))] for d in h.ACTIVITY}

    for j, drug in enumerate(h.ACTIVITY):
        test = np.array([drug in s for s in sets])
        single = np.array([s == (drug,) for s in sets])
        keep = [k for k in range(len(h.ACTIVITY)) if k != j]
        Dtr, Dte = D[~test][:, keep], D[test][:, keep]
        dj = D[test][:, j:j + 1]

        R, _ = h.fit(Dtr, X[~test])
        base = Dte @ R
        Rf, cf = h.fit(Dtr, X[~test], feats["ds_fold"][~test], False)
        arms = {
            "lr": base,
            "lr+ds_fold": Dte @ Rf + cf * feats["ds_fold"][test],
        }
        for name, rows_by_drug in (("lr+dsrow", dsrow), ("lr+reachrow", reachrow)):
            c = row_scale(R, keep, rows_by_drug)
            row = rows_by_drug.get(drug)
            arms[name] = base + (dj @ (c * row)[None, :] if row is not None else 0.0)
            cs[name][drug] = c
        row = dsrow.get(drug)
        c = cs["lr+dsrow"][drug]
        arms["lr+dialrow"] = base + (dj @ (c * alpha1[drug] * row)[None, :]
                                     if row is not None else 0.0)

        mask = np.ones(h.N_RESP, bool)
        for col in marker_cols[drug]:
            mask[col] = False
        for name, pred in arms.items():
            sp = pred[single[test]]
            res[name][drug] = (
                h.corr(X[test], pred),
                h.corr(X[single], sp) if single.any() else None,
                list(zip(h.calls(pred[:, :h.N_PROT]).ravel(),
                         h.calls(X[test][:, :h.N_PROT]).ravel())),
                h.corr(X[test][:, mask], pred[:, mask]),   # marker-excluded, the H2 view
            )

    out = {"alphas": {d: round(a, 3) for d, a in alpha1.items()},
           "c_dsrow_by_fold": {d: round(c, 4) for d, c in cs["lr+dsrow"].items()}}
    for arm, by in res.items():
        r = [by[d][0] for d in h.ACTIVITY]
        s = [by[d][1] for d in h.ACTIVITY if by[d][1] is not None]
        from train_116 import balanced_accuracy
        out[arm] = {"mean_corr": round(sum(r) / len(r), 4),
                    "single_mean_corr": round(sum(s) / len(s), 4) if s else None,
                    "balanced_accuracy": round(balanced_accuracy(
                        [p for d in h.ACTIVITY for p in by[d][2]]), 4),
                    "per_drug_corr": {d: round(by[d][0], 4) for d in h.ACTIVITY}}
    out["H1 lr+dsrow_minus_lr"] = h.boot(res, "lr+dsrow", "lr")
    out["H1' lr+dsrow_minus_lr+reachrow"] = h.boot(res, "lr+dsrow", "lr+reachrow")
    out["H2 lr+dialrow_minus_lr+dsrow (markers excluded)"] = h.boot(res, "lr+dialrow", "lr+dsrow", idx=3)
    out["continuity lr+ds_fold_minus_lr"] = h.boot(res, "lr+ds_fold", "lr")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
