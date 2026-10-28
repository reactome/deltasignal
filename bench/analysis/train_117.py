#!/usr/bin/env python3
"""specs/117: readout responsiveness, learned from the other drugs, with and without the frozen model's call.

For each readout, the share of its distinct training measurements that changed (binned) and their majority direction.
R1 tables (bin, direction); R2 adds the frozen model's call; R3 adds reachability's call instead. Each cell predicts the
class with the largest n(class in cell) / n(class in training) (specs/116's rule). Scored leave-one-drug-out within
each dataset; 95% bootstrap intervals over held-out drugs.

Usage: train_117.py KORKUT_CASES DECRYPTM_CASES
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
from train_116 import CLASSES, MIN_CELL, balanced_accuracy, drugs_of, macro_f1  # noqa: E402

RESP_EDGES = (0.1, 0.25, 0.5)


def resp_bin(changed, total):
    if total == 0:
        return "unseen"
    s = changed / total
    if s == 0:
        return 0
    return 1 + sum(s > e for e in RESP_EDGES)


def responsiveness(train_rows):
    """{readout: (bin, direction)} from distinct (group, readout) measurements in train_rows."""
    meas = {(r["group"], r["readout"]): r["measured"] for r in train_rows}
    by = collections.defaultdict(collections.Counter)
    for (_, ro), m in meas.items():
        by[ro][m] += 1
    out = {}
    for ro, c in by.items():
        d = "down" if c["down"] > c["up"] else "up" if c["up"] > c["down"] else "tie"
        out[ro] = (resp_bin(c["down"] + c["up"], sum(c.values())), d)
    return out


def fit(train, keys):
    """Nested tables for each prefix of keys, longest first: [(keys[:n], {cell: class})]."""
    total = collections.Counter(m for _, m in train)

    def best(counts):
        return max(CLASSES, key=lambda c: (counts[c] / total[c] if total[c] else 0, c == "not"))

    tables = []
    for n in range(len(keys), 1, -1):
        cells = collections.defaultdict(collections.Counter)
        for f, m in train:
            cells[tuple(f[k] for k in keys[:n])][m] += 1
        tables.append((keys[:n], {k: best(v) for k, v in cells.items() if sum(v.values()) >= MIN_CELL}))
    return tables


def predict(tables, f):
    for keys, table in tables:
        c = table.get(tuple(f[k] for k in keys))
        if c:
            return c
    return "not"


LAYERS = {"R1": ("bin", "dir"), "R2": ("bin", "dir", "model"), "R3": ("bin", "dir", "reach")}


def features(rows, resp):
    return [dict(r, bin=resp.get(r["readout"], ("unseen", "tie"))[0], dir=resp.get(r["readout"], ("unseen", "tie"))[1])
            for r in rows]


def lodo(rows, kind):
    for r in rows:
        r["drugs"] = drugs_of(r["group"], kind)
    preds = []
    for d in sorted(set().union(*(r["drugs"] for r in rows))):
        train = [r for r in rows if d not in r["drugs"]]
        test = [r for r in rows if d in r["drugs"]]
        resp = responsiveness(train)
        ftr, fte = features(train, resp), features(test, resp)
        models = {k: fit([(f, f["measured"]) for f in ftr], keys) for k, keys in LAYERS.items()}
        for f in fte:
            p = {k: predict(m, f) for k, m in models.items()}
            preds.append(dict(p, drugs={d}, measured=f["measured"], frozen=f["model"], reach=f["reach"],
                              always_not="not"))
    return preds


ARMS = ("R2", "R1", "R3", "frozen", "reach", "always_not")
PAIRS = (("R2", "R1"), ("R2", "R3"), ("R1", "frozen"), ("R2", "frozen"))


def summarize(preds, seed=117, n_boot=2000):
    out = {a: {"balanced_accuracy": balanced_accuracy([(p[a], p["measured"]) for p in preds]),
               "macro_f1": macro_f1([(p[a], p["measured"]) for p in preds])} for a in ARMS}
    by = collections.defaultdict(list)
    for p in preds:
        by[next(iter(p["drugs"]))].append(p)
    drugs = sorted(by)
    rng = random.Random(seed)
    diffs = collections.defaultdict(list)
    for _ in range(n_boot):
        sample = [p for d in (rng.choice(drugs) for _ in drugs) for p in by[d]]
        ba = {a: balanced_accuracy([(p[a], p["measured"]) for p in sample]) for a in {x for pr in PAIRS for x in pr}}
        for a, b in PAIRS:
            diffs[(a, b)].append(ba[a] - ba[b])
    for (a, b), v in diffs.items():
        v.sort()
        out[f"{a}_minus_{b}"] = [v[int(0.025 * n_boot)], v[int(0.975 * n_boot) - 1]]
    out["cases"], out["drugs"] = len(preds), len(drugs)
    return out


def load(cases_tsv, key):
    out = []
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored" or r["measured"] == "conflict":
            continue
        out.append({"group": r[key[0]], "readout": tuple(r[k] for k in key[1:]), "measured": r["measured"],
                    "model": r["predicted"], "reach": r["reach_predicted"]})
    return out


def main():
    kc, dc = sys.argv[1:3]
    dataset_roles.require_role("korkut2015", "fit")
    dataset_roles.require_role("decryptm1", "fit")
    res = {"lodo_korkut": summarize(lodo(load(kc, ("drug_set", "readout")), "korkut")),
           "lodo_decryptm": summarize(lodo(load(dc, ("drug", "uniprot", "position")), "decryptm"))}
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
