#!/usr/bin/env python3
"""specs/134: score a total-readout arm against the control on korkut, as registered.

On cases scored in both arms (pathway, drug set, readout): three-class balanced accuracy, arm minus control, with a 95%
interval from resampling drug sets (10,000, seed 134); accuracy and macro-F1 beside it; each also on total and phospho
readouts separately; change calls, class moves, and cases each arm leaves unscored (by status).

Usage: score_134.py CONTROL_CASES ARM_CASES [READOUTS_TSV] > out.json
"""
import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

CLASSES = ("up", "not", "down")
READOUTS = Path(__file__).resolve().parents[1] / "bench" / "analysis" / "062" / "readouts.tsv"


def load(path):
    return {(r["pathway_id"], r["drug_set"], r["readout"]): r for r in csv.DictReader(open(path), delimiter="\t")}


def kinds(path):
    return {r["readout"]: r["kind"] for r in csv.DictReader((l for l in open(path) if not l.startswith("#")),
                                                             delimiter="\t")}


def metrics(pairs):
    """pairs: [(predicted, measured)] -> balanced accuracy (mean recall over measured classes present), accuracy,
    macro-F1 over the three classes."""
    if not pairs:
        return {"n": 0}
    n_meas, n_pred, hit = Counter(m for _, m in pairs), Counter(p for p, _ in pairs), Counter(m for p, m in pairs if p == m)
    recalls = [hit[c] / n_meas[c] for c in CLASSES if n_meas[c]]
    f1 = []
    for c in CLASSES:
        prec = hit[c] / n_pred[c] if n_pred[c] else 0.0
        rec = hit[c] / n_meas[c] if n_meas[c] else 0.0
        f1.append(2 * prec * rec / (prec + rec) if prec + rec else 0.0)
    return {"n": len(pairs), "balanced": sum(recalls) / len(recalls), "accuracy": sum(hit.values()) / len(pairs),
            "macro_f1": sum(f1) / 3, "change_calls": sum(n_pred[c] for c in ("up", "down")),
            "change_calls_wrong": sum(1 for p, m in pairs if p != "not" and p != m)}


def compare(keys, ctl, arm, n_boot=10000, seed=134):
    pc = [(ctl[k]["predicted"], ctl[k]["measured"]) for k in keys]
    pa = [(arm[k]["predicted"], arm[k]["measured"]) for k in keys]
    mc, ma = metrics(pc), metrics(pa)
    out = {"control": mc, "arm": ma}
    if not keys:
        return out
    out["diff"] = {m: ma[m] - mc[m] for m in ("balanced", "accuracy", "macro_f1")}
    by = defaultdict(list)
    for k in keys:
        by[k[1]].append(k)
    groups = sorted(by)
    rng = random.Random(seed)
    draws = []
    for _ in range(n_boot):
        ks = [k for g in (rng.choice(groups) for _ in groups) for k in by[g]]
        draws.append(metrics([(arm[k]["predicted"], arm[k]["measured"]) for k in ks])["balanced"]
                     - metrics([(ctl[k]["predicted"], ctl[k]["measured"]) for k in ks])["balanced"])
    draws.sort()
    out["balanced_ci95"] = [draws[int(0.025 * n_boot)], draws[int(0.975 * n_boot) - 1]]
    moves = Counter((ctl[k]["predicted"], arm[k]["predicted"], ctl[k]["measured"]) for k in keys
                    if ctl[k]["predicted"] != arm[k]["predicted"])
    out["moved"] = sum(moves.values())
    out["fixed"] = sum(v for (c, a, m), v in moves.items() if a == m)
    out["broken"] = sum(v for (c, a, m), v in moves.items() if c == m)
    out["moves"] = {f"{c}->{a} (measured {m})": v for (c, a, m), v in moves.most_common()}
    return out


def main():
    ctl, arm = load(sys.argv[1]), load(sys.argv[2])
    kind = kinds(sys.argv[3] if len(sys.argv) > 3 else READOUTS)
    both = sorted(k for k in ctl if k in arm and ctl[k]["status"] == arm[k]["status"] == "scored"
                  and ctl[k]["measured"] in CLASSES)
    report = {"status_control": Counter(r["status"] for r in ctl.values()),
              "status_arm": Counter(r["status"] for r in arm.values()),
              "scored_control_not_arm": Counter(arm[k]["status"] if k in arm else "absent" for k in ctl
                                                if ctl[k]["status"] == "scored" and (k not in arm or arm[k]["status"] != "scored")),
              "all": compare(both, ctl, arm)}
    for kd in ("total", "phospho"):
        report[kd] = compare([k for k in both if kind.get(k[2]) == kd], ctl, arm, n_boot=1000)
    json.dump(report, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
