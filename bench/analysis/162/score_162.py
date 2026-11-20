#!/usr/bin/env python3
"""specs/162: score engagement pins vs full-knockout pins on CPPA.

Unit: (sample, antibody) for mapped readouts, EXCLUDING the drug's own
specs/139 markers (both arms). Measured call: delta-p >= log2 1.5 up,
<= -log2 1.5 down. Predicted: the mean over pathways of the clipped log2 fold
(as specs/121), called at log2 0.85 / log2 1.15. Gate (spec.md): moved calls
fixed > broken by a one-sided sign test (p < 0.05) AND reachable sensitivity
drops by at most 0.05. Numbers go to OUT_JSON (outside the repository, CPPA
terms); the terminal prints them for the operator, and only the verdict is
committed.

Usage: score_162.py CPPA_UNPACKED SPLIT_TSV MARKERS_TSV DRUG_TARGETS_TSV \
  ENGAGEMENT_TSV READOUTS_TSV OUT_JSON ARMS_TSV [ARMS_TSV ...]
"""
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent))

CLIP = math.log2(100)
T = math.log2(1.5)
LO, HI = math.log2(0.85), math.log2(1.15)
GUARD = 0.05


def pcall(l2):
    return "down" if l2 < LO else ("up" if l2 >= HI else "not")


def mcall(dp):
    return "up" if dp >= T else ("down" if dp <= -T else "not")


def sign_test_p(fixed, broken):
    """One-sided exact binomial P(X >= fixed | n, 0.5)."""
    n = fixed + broken
    if n == 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) - n * math.log(2)
            for k in range(fixed, n + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs))


def auroc(pos, neg):
    """Mann-Whitney via ranks (ties averaged)."""
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum, i = 0.0, 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        r = (i + 1 + j) / 2.0
        rank_sum += r * sum(1 for k in range(i, j) if allv[k][1] == 1)
        i = j
    n1, n0 = len(pos), len(neg)
    return (rank_sum - n1 * (n1 + 1) / 2.0) / (n1 * n0) if n1 and n0 else float("nan")


def score(samples, pred, eng, markers):
    """samples: [{drug, cell_line, values}]; pred: {(arm, drug, line, readout): log2};
    eng: {drug: {line: e}}; markers: {drug: set}. Returns the stats dict."""
    units = []
    for s in samples:
        d, c = s["drug"], s["cell_line"]
        if c not in eng.get(d, {}):
            continue
        for ab, dp in s["values"].items():
            if ab in markers.get(d, set()):
                continue
            b = pred.get(("baseline", d, "", ab))
            e = pred.get(("engagement", d, c, ab))
            if b is None or e is None:
                continue
            units.append({"m": mcall(dp), "b": pcall(b), "e": pcall(e), "bl2": b, "el2": e,
                          "drug": d, "dp": dp})
    moved = [u for u in units if u["b"] != u["e"]]
    fixed = sum(u["e"] == u["m"] for u in moved)
    broken = sum(u["b"] == u["m"] for u in moved)
    p = sign_test_p(fixed, broken)
    reach = [u for u in units if u["bl2"] != 0.0]

    def card(col, l2col):
        chg = [u for u in reach if u["m"] != "not"]
        nc = [u for u in reach if u["m"] == "not"]
        called = [u for u in reach if u[col] != "not"]
        tp = sum(u[col] == u["m"] for u in chg)
        per = []
        for m in ("down", "not", "up"):
            sub = [u for u in units if u["m"] == m]
            if sub:
                per.append(sum(u[col] == m for u in sub) / len(sub))
        return {"sensitivity": tp / len(chg) if chg else float("nan"),
                "specificity": sum(u[col] == "not" for u in nc) / len(nc) if nc else float("nan"),
                "precision": tp / len(called) if called else float("nan"),
                "balanced": sum(per) / len(per) if per else float("nan"),
                "auroc": auroc([abs(u[l2col]) for u in chg], [abs(u[l2col]) for u in nc])}
    cb, ce = card("b", "bl2"), card("e", "el2")
    guard_ok = (cb["sensitivity"] - ce["sensitivity"]) <= GUARD
    return {"units": len(units), "moved": len(moved), "fixed": fixed, "broken": broken,
            "sign_test_p": p, "baseline": cb, "engagement": ce,
            "guard_ok": guard_ok, "pass": (fixed > broken and p < 0.05 and guard_ok)}


def main():
    import cppa_data_121 as cd
    import dataset_roles
    for ds in ("cppa_train", "cppa_test"):       # specs/051: this score decides the gate
        dataset_roles.require_role(ds, "select")
    from engagement_162 import markers as load_markers
    (unpacked, split_tsv, markers_tsv, drug_targets, engagement_tsv,
     readouts_tsv, out_json, *arm_files) = sys.argv[1:]
    from korkut_bench import read_tsv
    mapped = {r["readout"] for r in read_tsv(readouts_tsv) if r["kind"] != "exclude"}
    markers = load_markers(markers_tsv, drug_targets)
    eng = defaultdict(dict)
    for r in csv.DictReader(open(engagement_tsv), delimiter="\t"):
        eng[r["drug"]][r["cell_line"]] = float(r["engagement"])
    acc = defaultdict(list)
    for path in arm_files:
        for r in csv.DictReader(open(path), delimiter="\t"):
            ui = float(r["pred_ui"])
            l2 = max(-CLIP, min(CLIP, math.log2(ui))) if ui > 0 else -CLIP
            acc[(r["arm"], r["drug"], r["cell_line"], r["readout"])].append(l2)
    pred = {k: sum(v) / len(v) for k, v in acc.items()}
    samples = []
    for portion in cd.PORTIONS:
        for s in cd.load(unpacked, split_tsv, portion, final=True):
            if s["cell_line"] in eng.get(s["drug"], {}):
                s["values"] = {a: v for a, v in s["values"].items() if a in mapped}
                samples.append(s)
    res = score(samples, pred, eng, markers)
    Path(out_json).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    print(f"\nGATE: {'PASS' if res['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
