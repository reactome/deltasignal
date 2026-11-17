#!/usr/bin/env python3
"""specs/165: precision and coverage of change calls under truth-blind label
filters (robust, confident), per axis, against the unfiltered baseline.

Usage: score_165.py LABELS_TSV CONFIDENCE_JSON CURATOR_CASES EXPERIMENTAL_CASES KORKUT_CASES
"""
import csv
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))

THRESH = 0.85
N_BOOT = 10_000
SEED = 165
CODE = {"0": "down", "1": "not", "2": "up"}


def mag_bin(fold):
    if fold == 1.0:
        return "baseline"
    a = abs(math.log2(fold)) if fold > 0 else math.inf
    return "le1" if a <= 1 else ("le3" if a <= 3 else "gt3")


def confidence(table, call, fold, converged):
    """specs/146 lookup from the no_basis level on (the record has no basis)."""
    if fold == 1.0:
        return None
    call = "not" if call == "normal" else call
    railed = "True" if (fold <= 0.01 or fold >= 99.0) else "False"
    conv = "True" if converged else "False"
    mag = mag_bin(fold)
    min_cell = int(table["min_cell"])
    for level, key in (("no_basis", f"{call}|{railed}|{mag}|{conv}"),
                       ("no_magnitude", f"{call}|{railed}|{conv}"),
                       ("no_converged", f"{call}|{railed}")):
        cell = table["table"][level].get(key)
        if cell and int(cell["n"]) >= min_cell:
            return float(cell["confidence"])
    cell = table["table"]["call_only"].get(call)
    return float(cell["confidence"]) if cell else None


def evaluate(units):
    """units: {axis: [{unit, robust, confident, correct}]} (change calls on reachable cases)."""
    out = {}
    rng = random.Random(SEED)
    filters = {"F0": lambda u: True, "F1 robust": lambda u: u["robust"],
               "F2 confident": lambda u: u["confident"],
               "F3 robust+confident": lambda u: u["robust"] and u["confident"]}
    for axis, rows in units.items():
        n0 = len(rows)
        res = {}
        for name, keep in filters.items():
            kept = [u for u in rows if keep(u)]
            res[name] = {"precision": (sum(u["correct"] for u in kept) / len(kept)) if kept else None,
                         "coverage": len(kept) / n0 if n0 else None, "n": len(kept)}
        by = defaultdict(list)
        for u in rows:
            by[u["unit"]].append(u)
        keys = sorted(by)
        deltas = []
        for _ in range(N_BOOT):
            pick = [r for k in (rng.choice(keys) for _ in keys) for r in by[k]]
            f3 = [u for u in pick if u["robust"] and u["confident"]]
            if not pick or not f3:
                continue
            deltas.append(sum(u["correct"] for u in f3) / len(f3) - sum(u["correct"] for u in pick) / len(pick))
        deltas.sort()
        if deltas:
            lo = deltas[int(0.025 * len(deltas))]
            hi = deltas[int(0.975 * len(deltas)) - 1]
        else:
            lo = hi = float("nan")
        p0, p3 = res["F0"]["precision"], res["F3 robust+confident"]["precision"]
        res["gate"] = {"delta": (p3 - p0) if (p0 is not None and p3 is not None) else None,
                       "lo": lo, "hi": hi, "pass": bool(deltas) and lo > 0}
        out[axis] = res
    return out


def load(labels_tsv, conf_json, curator_tsv, experimental_tsv, korkut_tsv):
    from holdout_report import TUNING_PATHWAYS
    from engagement_161 import marker_names
    table = json.load(open(conf_json))
    marks = marker_names()
    fragile = {}
    for r in csv.DictReader(open(labels_tsv), delimiter="\t"):
        fragile[(r["axis"], r["pathway"], r["case_a"], r["case_b"], r["readout"])] = r["fragile"] == "1"
    units = defaultdict(list)
    for axis, path in (("curator held-out", curator_tsv), ("experimental", experimental_tsv)):
        lab = "curator" if axis.startswith("curator") else "experimental"
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r.get("valid") != "1":
                continue
            if lab == "curator" and r["pathway"] in TUNING_PATHWAYS:
                continue
            fold = float(r["pred_ui"])
            call, truth = CODE[r["predicted"]], CODE[r["expected"]]
            if fold == 1.0 or call == "not":
                continue
            k = (lab, r["pathway"], r["gene"], r["direction"], r["key_output"])
            if k not in fragile:
                continue
            c = confidence(table, call, fold, r.get("converged") in ("1", "True", "true"))
            units[axis].append({"unit": (r["pathway"], r["gene"], r["direction"]),
                                "robust": not fragile[k],
                                "confident": c is not None and c >= THRESH,
                                "correct": call == truth})
    for r in csv.DictReader(open(korkut_tsv), delimiter="\t"):
        if r["status"] != "scored":
            continue
        if any(r["readout"] in marks.get(d, set()) for d in r["drug_set"].split("+")):
            continue
        fold, call = float(r["pred_ui"]), r["predicted"]
        if fold == 1.0 or call == "not":
            continue
        k = ("korkut", r["pathway_id"], r["drug_set"], "", r["readout"])
        if k not in fragile:
            continue
        c = confidence(table, call, fold, r.get("converged") in ("1", "True", "true"))
        units["korkut"].append({"unit": r["drug_set"], "robust": not fragile[k],
                                "confident": c is not None and c >= THRESH,
                                "correct": call == r["measured"]})
    return units


def main():
    import dataset_roles
    dataset_roles.require_role("korkut2015", "fit")
    res = evaluate(load(*sys.argv[1:6]))
    for axis, r in res.items():
        print(f"--- {axis}")
        for name in ("F0", "F1 robust", "F2 confident", "F3 robust+confident"):
            x = r[name]
            p = "n/a" if x["precision"] is None else f"{x['precision']:.3f}"
            print(f"   {name:>20}: precision {p}  coverage {x['coverage']:.3f}  (n={x['n']})")
        g = r["gate"]
        print(f"   gate F3 - F0: {g['delta']:+.3f} [{g['lo']:+.3f}, {g['hi']:+.3f}] -> "
              f"{'PASS' if g['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
