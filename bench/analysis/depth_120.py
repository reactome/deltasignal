#!/usr/bin/env python3
"""specs/120: the drug inhibition depth (the UI a drug target is pinned at), fitted from the other drugs.

One case log per depth (the frozen arm at 0 and three arms at 0.1, 0.25, 0.5). Cases scored in every arm. Leave one
drug out: the depth with the highest balanced accuracy of the model's call on the training drugs predicts the held-out
drug (ties: the smaller depth). H3 re-runs specs/117's R2 layer with the fitted-depth call in place of the frozen one.
95% bootstrap over held-out drugs.

Usage: depth_120.py KIND DEPTH=CASES_TSV...   (KIND korkut or decryptm; one DEPTH=path per arm, 0 among them)
"""
import collections
import csv
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dataset_roles  # noqa: E402
import train_117 as t117  # noqa: E402
from train_116 import balanced_accuracy, drugs_of, macro_f1  # noqa: E402

KEYS = {"korkut": ("pathway_id", "drug_set", "readout"),
        "decryptm": ("pathway_id", "drug", "uniprot", "position", "cell_lines")}


def load(path, kind):
    """{case key: row} for scored cases whose measured call is not a conflict."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["status"] != "scored" or r["measured"] == "conflict":
            continue
        k = KEYS[kind]
        out[tuple(r[c] for c in k)] = {"group": r[k[1]], "readout": tuple(r[c] for c in k[2:4] if c != "cell_lines"),
                                       "measured": r["measured"], "model": r["predicted"],
                                       "reach": r["reach_predicted"]}
    return out


def merge(arms):
    """Rows scored in every arm, each with calls {depth: class}; and how many each arm dropped."""
    common = set.intersection(*(set(a) for a in arms.values()))
    dropped = {d: len(a) - len(common) for d, a in arms.items()}
    base = arms[min(arms)]
    rows = []
    for k in sorted(common):
        r = dict(base[k], calls={d: a[k]["model"] for d, a in arms.items()})
        assert all(a[k]["measured"] == r["measured"] for a in arms.values()), k
        rows.append(r)
    return rows, dropped


def choose(train, depths):
    """The depth whose call scores the highest balanced accuracy on train; ties go to the smaller depth."""
    return max(depths, key=lambda d: (balanced_accuracy([(r["calls"][d], r["measured"]) for r in train]), -d))


def r2(train, test):
    """specs/117's R2 layer (responsiveness + the rows' `model` call) fitted on train, predicting test."""
    resp = t117.responsiveness(train)
    ftr, fte = t117.features(train, resp), t117.features(test, resp)
    tables = t117.fit([(f, f["measured"]) for f in ftr], t117.LAYERS["R2"])
    return [t117.predict(tables, f) for f in fte]


def lodo(rows, kind):
    depths = sorted(rows[0]["calls"])
    frozen = min(depths)
    for r in rows:
        r["drugs"] = drugs_of(r["group"], kind)
    preds, chosen = [], {}
    for d in sorted(set().union(*(r["drugs"] for r in rows))):
        train = [r for r in rows if d not in r["drugs"]]
        test = [r for r in rows if d in r["drugs"]]
        g = chosen[d] = choose(train, depths)
        at = lambda rs, x: [dict(r, model=r["calls"][x]) for r in rs]  # noqa: E731
        r2_fit = r2(at(train, g), at(test, g))
        r2_frozen = r2(at(train, frozen), at(test, frozen))
        for r, a, b in zip(test, r2_fit, r2_frozen):
            preds.append({"drug": d, "measured": r["measured"], "fitted": r["calls"][g], "frozen": r["calls"][frozen],
                          "R2_fitted": a, "R2_frozen": b, "reach": r["reach"], "always_not": "not",
                          **{f"depth_{x:g}": r["calls"][x] for x in depths}})
    return preds, chosen


PAIRS = (("fitted", "frozen"), ("R2_fitted", "R2_frozen"))


def false_change(pairs):
    """Share of change calls (down/up) made where nothing was measured."""
    calls = [m for p, m in pairs if p != "not"]
    return sum(m == "not" for m in calls) / len(calls) if calls else 0.0


def summarize(preds, seed=120, n_boot=2000):
    arms = [k for k in preds[0] if k not in ("drug", "measured")]
    out = {a: {"balanced_accuracy": balanced_accuracy([(p[a], p["measured"]) for p in preds]),
               "macro_f1": macro_f1([(p[a], p["measured"]) for p in preds]),
               "false_change": false_change([(p[a], p["measured"]) for p in preds])} for a in arms}
    by = collections.defaultdict(list)
    for p in preds:
        by[p["drug"]].append(p)
    drugs = sorted(by)
    rng = random.Random(seed)
    diffs = collections.defaultdict(list)
    for _ in range(n_boot):
        sample = [p for d in (rng.choice(drugs) for _ in drugs) for p in by[d]]
        ba = {a: balanced_accuracy([(p[a], p["measured"]) for p in sample]) for pr in PAIRS for a in pr}
        for a, b in PAIRS:
            diffs[(a, b)].append(ba[a] - ba[b])
    for (a, b), v in diffs.items():
        v.sort()
        out[f"{a}_minus_{b}"] = {"point": out[a]["balanced_accuracy"] - out[b]["balanced_accuracy"],
                                 "ci95": [v[int(0.025 * n_boot)], v[int(0.975 * n_boot) - 1]]}
    out["cases"], out["drugs"] = len(preds), len(drugs)
    return out


def main():
    kind, specs = sys.argv[1], sys.argv[2:]
    dataset_roles.require_role({"korkut": "korkut2015", "decryptm": "decryptm1"}[kind], "fit")
    arms = {float(s.split("=", 1)[0]): load(s.split("=", 1)[1], kind) for s in specs}
    assert 0.0 in arms, "the frozen arm (depth 0) is the control"
    rows, dropped = merge(arms)
    preds, chosen = lodo(rows, kind)
    print(json.dumps({"kind": kind, "dropped": dropped, "chosen_depth": chosen, **summarize(preds)}, indent=1))


if __name__ == "__main__":
    main()
