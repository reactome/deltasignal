#!/usr/bin/env python3
"""specs/128: score an edge-curve arm against the control on the curator, experimental and korkut axes.

Curator and experimental use specs/124's comparison (class changes, McNemar, rails, convergence, magnitudes). korkut
compares the arm and the control on the cases scored in both with a measured call (specs/062): right, change calls,
macro-F1, and McNemar on the cases whose correctness differs.

Usage: score_128.py CONTROL_DIR ARM_DIR MILD_DIR
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_124  # noqa: E402
from _common import mcnemar_exact  # noqa: E402

CALLS = ("up", "down", "not")


def load_korkut(path):
    with open(path, newline="") as fh:
        return {(r["pathway_id"], r["drug_set"], r["readout"]): r for r in csv.DictReader(fh, delimiter="\t")
                if r["status"] == "scored" and r["measured"] in CALLS}


def macro_f1(pairs):
    f1 = []
    for c in CALLS:
        tp = sum(p == c and m == c for m, p in pairs)
        fp = sum(p == c and m != c for m, p in pairs)
        fn = sum(p != c and m == c for m, p in pairs)
        f1.append(2 * tp / (2 * tp + fp + fn) if tp else 0.0)
    return sum(f1) / len(CALLS)


def korkut(ctrl, arm):
    keys = sorted(ctrl.keys() & arm.keys())
    out = {"cases": len(keys), "always_not_right": sum(ctrl[k]["measured"] == "not" for k in keys)}
    for name, d in (("control", ctrl), ("arm", arm)):
        pairs = [(d[k]["measured"], d[k]["predicted"]) for k in keys]
        out[f"{name}_right"] = sum(m == p for m, p in pairs)
        out[f"{name}_change_calls"] = sum(p != "not" for _, p in pairs)
        out[f"{name}_macro_f1"] = macro_f1(pairs)
    fixed = sum(ctrl[k]["predicted"] != ctrl[k]["measured"] and arm[k]["predicted"] == arm[k]["measured"] for k in keys)
    broke = sum(ctrl[k]["predicted"] == ctrl[k]["measured"] and arm[k]["predicted"] != arm[k]["measured"] for k in keys)
    out.update(fixed=fixed, broke=broke, mcnemar_p=mcnemar_exact(fixed, broke))
    return out


def main():
    control, arm, mild = sys.argv[1:4]
    report = score_124.arm_report(control, arm, 100.0, mild)
    report["korkut"] = korkut(load_korkut(Path(control) / "korkut_cases.tsv"),
                              load_korkut(Path(arm) / "korkut_cases.tsv"))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
