#!/usr/bin/env python3
"""specs/116: a trained readout layer on the frozen model's case logs.

A lookup table over (fold bin of the frozen model's predicted fold, hop band); each cell predicts the class with the
largest n(class in cell) / n(class in training), i.e. the balanced-accuracy-optimal call. Cells with fewer than
MIN_CELL training cases fall back to their hop band alone, and bands below that to NOT. The topology-only control uses
the hop band alone. Scored leave-one-drug-out within each dataset and across datasets, against the frozen model,
reachability, the topology-only layer and always-NOT; 95% bootstrap intervals over held-out drugs.

Usage: train_116.py KORKUT_CASES KORKUT_HOPS DECRYPTM_CASES DECRYPTM_HOPS
"""
import collections
import csv
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dataset_roles  # noqa: E402

CLASSES = ("down", "not", "up")
MIN_CELL = 10
BANDS = ((1, 5), (6, 10), (11, 20), (21, 40), (41, 10 ** 9))
FOLD_EDGES = (-3.32, -1.0, -0.585, 0.585, 1.0, 3.32)


def fold_bin(pred_ui):
    """Bin of log2(pred_ui): 'zero' for 0, else 0..6 by FOLD_EDGES (left-closed except the middle bin)."""
    x = float(pred_ui)
    if x <= 0:
        return "zero"
    l = math.log2(x)
    if l <= -3.32:
        return 0
    if l <= -1.0:
        return 1
    if l <= -0.585:
        return 2
    if l < 0.585:
        return 3
    if l < 1.0:
        return 4
    if l < 3.32:
        return 5
    return 6


def hop_band(h):
    if h in ("", None):
        return "none"
    h = int(h)
    return next(i for i, (lo, hi) in enumerate(BANDS) if lo <= h <= hi)


def fit(train, features):
    """{feature tuple: class} and a band fallback, from training cases [(features dict, measured)]."""
    total = collections.Counter(m for _, m in train)

    def best(counts):
        return max(CLASSES, key=lambda c: (counts[c] / total[c] if total[c] else 0, c == "not"))

    cells, bands = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
    for f, m in train:
        cells[tuple(f[k] for k in features)][m] += 1
        bands[f["band"]][m] += 1
    table = {k: best(v) for k, v in cells.items() if sum(v.values()) >= MIN_CELL}
    fallback = {k: best(v) for k, v in bands.items() if sum(v.values()) >= MIN_CELL}
    return table, fallback


def predict(model, f, features):
    table, fallback = model
    return table.get(tuple(f[k] for k in features), fallback.get(f["band"], "not"))


def balanced_accuracy(pairs):
    """pairs: [(predicted, measured)]; mean recall over the measured classes present."""
    rec = []
    for c in CLASSES:
        n = [p for p, m in pairs if m == c]
        if n:
            rec.append(sum(p == c for p in n) / len(n))
    return sum(rec) / len(rec)


def macro_f1(pairs):
    out = []
    for c in CLASSES:
        tp = sum(p == m == c for p, m in pairs)
        pp, ap = sum(p == c for p, _ in pairs), sum(m == c for _, m in pairs)
        pr, rc = (tp / pp if pp else 0), (tp / ap if ap else 0)
        out.append(2 * pr * rc / (pr + rc) if pr + rc else 0)
    return sum(out) / 3


def load(cases_tsv, hops_tsv, key):
    hops = {tuple(r[k] for k in key): r["hops"] for r in csv.DictReader(open(hops_tsv), delimiter="\t")}
    out = []
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored" or r["measured"] == "conflict":
            continue
        k = tuple(r[x] for x in key)
        if k not in hops:
            raise SystemExit(f"no hop distance for {k}: rebuild the hops table")
        out.append({"fold": fold_bin(r["pred_ui"]), "band": hop_band(hops[k]), "measured": r["measured"],
                    "model": r["predicted"], "reach": r["reach_predicted"], "group": r[key[1]]})
    return out


def drugs_of(group, kind):
    return set(group.split("+")) if kind == "korkut" else {group}


def evaluate(train_rows, test_rows):
    """Predictions of every arm on test_rows, layers fitted on train_rows."""
    tr = [(r, r["measured"]) for r in train_rows]
    full, topo = fit(tr, ("fold", "band")), fit(tr, ("band",))
    return [{"drugs": r["drugs"], "measured": r["measured"], "trained": predict(full, r, ("fold", "band")),
             "topology": predict(topo, r, ("band",)), "frozen": r["model"], "reach": r["reach"], "always_not": "not"}
            for r in test_rows]


def lodo(rows, kind):
    for r in rows:
        r["drugs"] = drugs_of(r["group"], kind)
    preds = []
    for d in sorted(set().union(*(r["drugs"] for r in rows))):
        test = [dict(r, drugs={d}) for r in rows if d in r["drugs"]]
        preds += evaluate([r for r in rows if d not in r["drugs"]], test)
    return preds


def transfer(train_rows, test_rows, test_kind):
    for r in test_rows:
        r["drugs"] = drugs_of(r["group"], test_kind)
    return evaluate(train_rows, [dict(r, drugs={min(r["drugs"])}) for r in test_rows])


ARMS = ("trained", "frozen", "reach", "topology", "always_not")


def summarize(preds, seed=116, n_boot=2000):
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
        base = balanced_accuracy([(p["trained"], p["measured"]) for p in sample])
        for a in ARMS[1:]:
            diffs[a].append(base - balanced_accuracy([(p[a], p["measured"]) for p in sample]))
    for a in ARMS[1:]:
        v = sorted(diffs[a])
        out[f"trained_minus_{a}"] = [v[int(0.025 * n_boot)], v[int(0.975 * n_boot) - 1]]
    out["cases"], out["drugs"] = len(preds), len(drugs)
    return out


def main():
    kc, kh, dc, dh = sys.argv[1:5]
    dataset_roles.require_role("korkut2015", "fit")
    dataset_roles.require_role("decryptm1", "fit")
    kor = load(kc, kh, ("pathway_id", "drug_set", "readout"))
    dec = load(dc, dh, ("pathway_id", "drug", "uniprot", "position"))
    res = {"lodo_korkut": summarize(lodo(kor, "korkut")), "lodo_decryptm": summarize(lodo(dec, "decryptm")),
           "korkut_to_decryptm": summarize(transfer(kor, dec, "decryptm")),
           "decryptm_to_korkut": summarize(transfer(dec, kor, "korkut"))}
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
