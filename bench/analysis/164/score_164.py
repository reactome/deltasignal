#!/usr/bin/env python3
"""specs/164: leave-one-dataset-out selection and test of the specs/129 curve
grid on korkut, decryptM and CPPA.

Per dataset, per arm: a unit's predicted call (and log2 fold), against one
measured call per unit. Selection for held-out D: the arm with the highest
mean balanced accuracy over the OTHER two datasets (ties: identity, then larger
K, then larger h). Test on D, chosen vs identity: moved calls fixed > broken by
a one-sided sign test (p < 0.05) AND reachable sensitivity (identity log2 != 0)
drops by at most 0.05. Generalizes only if all three held-out tests pass.
CPPA-derived numbers go to OUT_JSON (outside the repository).

Usage: score_164.py RUNS_DIR CPPA_UNPACKED SPLIT_TSV MARKERS_TSV DRUG_TARGETS_TSV OUT_JSON
"""
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))
sys.path.insert(0, str(HERE.parents[1] / "162"))

CLIP = math.log2(100)
T = math.log2(1.5)
LO, HI = math.log2(0.85), math.log2(1.15)
GUARD = 0.05
IDENTITY = "cf_off"


def l2(ui):
    return max(-CLIP, min(CLIP, math.log2(ui))) if ui > 0 else -CLIP


def pcall(v):
    return "down" if v < LO else ("up" if v >= HI else "not")


def sign_p(fixed, broken):
    n = fixed + broken
    if n == 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) - n * math.log(2)
            for k in range(fixed, n + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs))


def arm_key(arm):
    """Tie-break order: identity first, then larger K, then larger h."""
    if arm == IDENTITY:
        return (0, 0.0, 0.0)
    k = arm.split("_K")[1].split("_h")[0]
    h = arm.split("_h")[1]
    kv = math.inf if k == "inf" else float(k)
    hv = {"05": 0.5, "075": 0.75, "1": 1.0}[h]
    return (1, -kv, -hv)


def balanced(units, calls):
    per = []
    for m in ("down", "not", "up"):
        sub = [u for u, t in units.items() if t == m]
        if sub:
            per.append(sum(calls[u][0] == m for u in sub) / len(sub))
    return sum(per) / len(per) if per else float("nan")


def card(units, calls, ident):
    reach = [u for u in units if ident[u][1] != 0.0]
    chg = [u for u in reach if units[u] != "not"]
    nc = [u for u in reach if units[u] == "not"]
    called = [u for u in reach if calls[u][0] != "not"]
    tp = sum(calls[u][0] == units[u] for u in chg)
    return {"sensitivity": tp / len(chg) if chg else float("nan"),
            "specificity": sum(calls[u][0] == "not" for u in nc) / len(nc) if nc else float("nan"),
            "precision": tp / len(called) if called else float("nan"),
            "balanced": balanced(units, calls)}


def lodo(data):
    """data: {dataset: {"truth": {unit: call}, "arms": {arm: {unit: (call, l2)}}}}."""
    out = {}
    for D in sorted(data):
        others = [x for x in data if x != D]
        arms = sorted(set.intersection(*(set(data[x]["arms"]) for x in data)), key=arm_key)
        obj = {a: sum(balanced(data[x]["truth"], data[x]["arms"][a]) for x in others) / len(others)
               for a in arms}
        best = max(obj.values())
        chosen = min((a for a in arms if obj[a] == best), key=arm_key)
        truth = data[D]["truth"]
        ident = data[D]["arms"][IDENTITY]
        res = {"chosen": chosen, "objective": obj}
        if chosen == IDENTITY:
            res.update({"tested": False, "pass": False})
        else:
            calls = data[D]["arms"][chosen]
            moved = [u for u in truth if calls[u][0] != ident[u][0]]
            fixed = sum(calls[u][0] == truth[u] for u in moved)
            broken = sum(ident[u][0] == truth[u] for u in moved)
            p = sign_p(fixed, broken)
            ci, cc = card(truth, ident, ident), card(truth, calls, ident)
            guard = (ci["sensitivity"] - cc["sensitivity"]) <= GUARD
            res.update({"tested": True, "moved": len(moved), "fixed": fixed, "broken": broken,
                        "sign_p": p, "identity": ci, "chosen_card": cc, "guard_ok": guard,
                        "pass": fixed > broken and p < 0.05 and guard})
        out[D] = res
    out["generalizes"] = all(out[D]["pass"] for D in data)
    return out


def load_korkut(runs, arms):
    from engagement_161 import marker_names
    marks = marker_names()
    truth, arm_calls = {}, defaultdict(dict)
    for a in arms:
        for r in csv.DictReader(open(runs / a / "korkut" / "korkut_cases.tsv"), delimiter="\t"):
            if r["status"] != "scored":
                continue
            if any(r["readout"] in marks.get(d, set()) for d in r["drug_set"].split("+")):
                continue
            u = (r["pathway_id"], r["drug_set"], r["readout"])
            truth[u] = r["measured"]
            v = l2(float(r["pred_ui"]))
            arm_calls[a][u] = (r["predicted"], v)
    return truth, arm_calls


def load_decryptm(runs, arms):
    truth, arm_calls = {}, defaultdict(dict)
    for a in arms:
        for r in csv.DictReader(open(runs / a / "decryptm" / "decryptm_cases.tsv"), delimiter="\t"):
            if r["status"] != "scored":
                continue
            if r["site_gene"] in set(r["targets"].split(";")):
                continue
            u = (r["pathway_id"], r["drug"], r["uniprot"], r["position"])
            truth[u] = r["measured"]
            arm_calls[a][u] = (r["predicted"], l2(float(r["pred_ui"])))
    return truth, arm_calls


def load_cppa(runs, arms, unpacked, split_tsv, markers_tsv, drug_targets):
    import cppa_data_121 as cd
    from engagement_162 import markers as load_markers
    marks = load_markers(markers_tsv, drug_targets)
    preds = {}
    for a in arms:
        acc = defaultdict(list)
        for r in csv.DictReader(open(runs / a / "cppa" / "cppa_predict_cases.tsv"), delimiter="\t"):
            if r["status"] == "scored":
                acc[(r["drug"], r["readout"])].append(l2(float(r["pred_ui"])))
        preds[a] = {k: sum(v) / len(v) for k, v in acc.items()}
    keys = set.intersection(*(set(p) for p in preds.values()))
    by_drug = defaultdict(set)
    for d, ro in keys:
        by_drug[d].add(ro)
    truth, arm_calls = {}, defaultdict(dict)
    for portion in cd.PORTIONS:
        for s in cd.load(unpacked, split_tsv, portion, final=True):
            d = s["drug"]
            for ro in by_drug.get(d, ()):
                if ro in marks.get(d, set()) or ro not in s["values"]:
                    continue
                dp = s["values"][ro]
                u = (s["set"], s["uid"], ro)
                truth[u] = "up" if dp >= T else ("down" if dp <= -T else "not")
                for a in arms:
                    v = preds[a][(d, ro)]
                    arm_calls[a][u] = (pcall(v), v)
    return truth, arm_calls


def main():
    import dataset_roles
    dataset_roles.require_role("korkut2015", "fit")
    runs = Path(sys.argv[1])
    unpacked, split_tsv, markers_tsv, drug_targets, out_json = sys.argv[2:7]
    arms = sorted(p.name for p in runs.iterdir() if (p / "cppa" / "cppa_predict_cases.tsv").exists())
    if IDENTITY not in arms:
        raise SystemExit("the identity arm cf_off is missing")
    data = {}
    for name, loader in (("korkut", lambda: load_korkut(runs, arms)),
                         ("decryptm", lambda: load_decryptm(runs, arms))):
        t, ac = loader()
        common = set.intersection(*(set(ac[a]) for a in arms)) & set(t)
        data[name] = {"truth": {u: t[u] for u in common},
                      "arms": {a: {u: ac[a][u] for u in common} for a in arms}}
    t, ac = load_cppa(runs, arms, unpacked, split_tsv, markers_tsv, drug_targets)
    data["cppa"] = {"truth": t, "arms": dict(ac)}
    res = lodo(data)
    for D in ("korkut", "decryptm", "cppa"):
        res[D]["units"] = len(data[D]["truth"])
    Path(out_json).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps(res, indent=1, default=str))
    print(f"\nARMS: {len(arms)}  GENERALIZES: {res['generalizes']}")


if __name__ == "__main__":
    main()
