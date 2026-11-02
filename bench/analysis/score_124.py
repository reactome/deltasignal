#!/usr/bin/env python3
"""specs/124: how much the results depend on the 100x ceiling.

Each run directory holds `curator_cases.tsv` and `experimental_cases.tsv` (benchmark_vs_mpbiopath.py). `pred_ui` is the
readout's fold over rest in every arm (the bench divides by DS_BASELINE), so arms with different ceilings compare
directly.

  d0   CONTROL FROZEN                 control must reproduce the frozen run case for case (status, call, value)
  arm  CONTROL ARM CEILING MILD       an arm against the control: class changes, rails, magnitudes, convergence;
                                      MILD is the control at a mild dose (specs/108), for rail-dependence

Usage: score_124.py d0 CONTROL_DIR FROZEN_DIR
       score_124.py arm CONTROL_DIR ARM_DIR CEILING MILD_DIR
"""
import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import mcnemar_exact  # noqa: E402
from holdout_report import TUNING_PATHWAYS  # noqa: E402

AXES = ("curator", "experimental")
RAF = "RAF_MAP_kinase_cascade"
MILD_LO, MILD_HI = 0.01, 99.0   # specs/108: a mild-dose readout at a rail
MOVE = 0.10                     # P5: a fold that moves by more than 10%


def load(path):
    """All rows with a scorable expectation, keyed by case."""
    with open(path, newline="") as fh:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(fh, delimiter="\t") if r["expected"] in ("0", "1", "2")}


def d0(control, frozen):
    """Differences between two runs that should be identical: keys, valid, predicted, pred_ui (as printed), converged."""
    out = {}
    for axis in AXES:
        a, b = load(Path(control) / f"{axis}_cases.tsv"), load(Path(frozen) / f"{axis}_cases.tsv")
        diff = Counter()
        diff["only_control"], diff["only_frozen"] = len(a.keys() - b.keys()), len(b.keys() - a.keys())
        for k in a.keys() & b.keys():
            for col in ("valid", "predicted", "pred_ui", "converged"):
                if a[k].get(col) != b[k].get(col):
                    diff[col] += 1
        out[axis] = {"cases": len(a.keys() & b.keys()), **diff}
    return out


def split_of(axis, pathway):
    if axis == "experimental":
        return "experimental"
    return "tuning" if pathway in TUNING_PATHWAYS else "held-out"


def railed(fold, ceiling):
    return fold >= 0.99 * ceiling or fold <= 0.01


def moved(a, b):
    """Whether fold b differs from fold a by more than MOVE (relative); two zeros do not move."""
    if a == 0.0:
        return b != 0.0
    return abs(b / a - 1.0) > MOVE


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for t in range(i, j + 1):
            r[order[t]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(xs, ys):
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx, vy = sum((a - mx) ** 2 for a in rx), sum((b - my) ** 2 for b in ry)
    return cov / math.sqrt(vx * vy) if vx and vy else float("nan")


def compare(axis, ctrl, arm, ceiling, mild, set_aside=()):
    """Counts over cases valid in both runs, per split and direction ('0' knockout, '2' over-expression)."""
    keys = [k for k in ctrl.keys() & arm.keys()
            if ctrl[k]["valid"] == "1" and arm[k]["valid"] == "1" and k[0] not in set_aside]
    out = {}
    for k in keys:
        s = split_of(axis, k[0])
        for key in (s, f"{s}|dir{k[2]}"):
            g = out.setdefault(key, Counter())
            c, a = ctrl[k], arm[k]
            g["cases"] += 1
            g["control_right"] += c["predicted"] == c["expected"]
            g["arm_right"] += a["predicted"] == a["expected"]
            if c["predicted"] != a["predicted"]:
                g["changed"] += 1
                g["fixed"] += a["predicted"] == a["expected"]
                g["broke"] += c["predicted"] == c["expected"]
                if k in mild:
                    g["changed_mild_known"] += 1
                    g["changed_rail_dependent"] += railed_mild(mild[k])
            fc, fa = float(c["pred_ui"]), float(a["pred_ui"])
            g["control_railed"] += railed(fc, 100.0)
            g["arm_railed"] += railed(fa, ceiling)
            g["moved"] += moved(fc, fa)
            g["control_converged"] += c["converged"] == "1"
            g["arm_converged"] += a["converged"] == "1"
    for g in out.values():
        g["mcnemar_p"] = mcnemar_exact(g["fixed"], g["broke"])
    return out, keys


def railed_mild(row):
    v = float(row["pred_ui"])
    return v <= MILD_LO or v >= MILD_HI


def per_pathway(ctrl, arm, keys):
    net = Counter()
    for k in keys:
        c, a = ctrl[k], arm[k]
        net[k[0]] += (a["predicted"] == a["expected"]) - (c["predicted"] == c["expected"])
    return {p: n for p, n in sorted(net.items(), key=lambda x: x[1]) if n}


def magnitudes(ctrl, arm, keys):
    xs = [math.log(max(float(ctrl[k]["pred_ui"]), 1e-6)) for k in keys]
    ys = [math.log(max(float(arm[k]["pred_ui"]), 1e-6)) for k in keys]
    return {"spearman_log_fold": spearman(xs, ys) if len(keys) > 2 else float("nan"),
            "moved_gt_10pct": sum(moved(float(ctrl[k]["pred_ui"]), float(arm[k]["pred_ui"])) for k in keys),
            "readouts": len(keys)}


def arm_report(control, armdir, ceiling, milddir):
    out = {}
    for axis in AXES:
        ctrl, arm = load(Path(control) / f"{axis}_cases.tsv"), load(Path(armdir) / f"{axis}_cases.tsv")
        mild = {k: r for k, r in load(Path(milddir) / f"{axis}_cases.tsv").items() if r["valid"] == "1"}
        res, keys = compare(axis, ctrl, arm, ceiling, mild)
        res_raf, _ = compare(axis, ctrl, arm, ceiling, mild, set_aside=(RAF,))
        out[axis] = {"all": {k: dict(v) for k, v in res.items()},
                     "raf_set_aside": {k: dict(v) for k, v in res_raf.items()},
                     "per_pathway_net": per_pathway(ctrl, arm, keys),
                     "magnitudes": magnitudes(ctrl, arm, keys)}
    return out


def main():
    if len(sys.argv) >= 4 and sys.argv[1] == "d0":
        print(json.dumps(d0(sys.argv[2], sys.argv[3]), indent=1))
    elif len(sys.argv) == 6 and sys.argv[1] == "arm":
        print(json.dumps(arm_report(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]), indent=1))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
