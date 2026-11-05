#!/usr/bin/env python3
"""specs/139: gate a CPPA sample by its drug target's resting activity in that cell line, chosen by CV over cell lines.

A sample's activity A(s) is the mean, over its drug's marker antibodies (markers.py), of the marker's control level p0
ranked among the training cell lines of the same set (sets with >= MIN_LINES lines only; CPPA sets are normalised
separately). Under a gate at theta a sample with A(s) < theta predicts no change for every antibody. Every arm leaves
the drug's own markers out of scoring. Judge: three-class balanced accuracy, theta chosen per fold of training cell
lines on the other folds, gain bootstrapped over drugs, against a random gate of the same size per drug; the per-drug
correlation must not fall. Training portion only; CPPA values never leave this machine.

Usage: score_139.py CPPA_UNPACKED SPLIT_TSV READOUTS_TSV MARKERS_TSV CASES_TSV [CASES_TSV ...] > out.json
"""
import collections
import csv
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cppa_data_121  # noqa: E402
import score_121  # noqa: E402
import score_130  # noqa: E402

THETAS = (None, 0.2, 0.4, 0.6)                    # None: the identity (no gate)
MIN_LINES = 4
SEED = 139
FOLDS = 5
N_BOOT = 2000
N_NULL = 200
CLASSES = score_130.CLASSES
cls = score_130.cls


def read_markers(path):
    """{drug: set of marker antibodies}; drugs without markers are left out."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["markers"]:
            out[r["drug"]] = set(r["markers"].split(","))
    return out


def ranks(baselines, min_lines=MIN_LINES):
    """{(set, cell line): {antibody: rank in [0, 1]}} among the lines of the same set that measured the antibody;
    sets with fewer than min_lines lines are left out. Ties take the mid-rank."""
    by_set = collections.defaultdict(dict)
    for (s, c), vals in baselines.items():
        by_set[s][c] = vals
    out = {}
    for s, lines in by_set.items():
        if len(lines) < min_lines:
            continue
        per_ab = collections.defaultdict(list)
        for c, vals in lines.items():
            for a, v in vals.items():
                per_ab[a].append((v, c))
        for a, vc in per_ab.items():
            n = len(vc)
            if n < min_lines:
                continue
            vals = sorted(v for v, _ in vc)
            for v, c in vc:
                lo = vals.index(v)
                hi = n - 1 - vals[::-1].index(v)
                out.setdefault((s, c), {})[a] = ((lo + hi) / 2) / (n - 1)
    return out


def activity(sample, marker_set, rank_table):
    """Mean rank of the markers measured at rest in the sample's (set, cell line), or None."""
    r = rank_table.get((sample["set"], sample["cell_line"]), {})
    vals = [r[a] for a in marker_set if a in r]
    return sum(vals) / len(vals) if vals else None


def tally(pairs):
    """Confusion tallies for one sample: (total by measured class, hits by the prediction, hits if gated)."""
    tot, hit = collections.Counter(), collections.Counter()
    for pred, meas in pairs:
        m = cls(meas)
        tot[m] += 1
        hit[m] += m == cls(pred)
    return tot, hit, collections.Counter({0: tot[0]})


def balanced(tot, hit):
    present = [c for c in CLASSES if tot[c]]
    return sum(hit[c] / tot[c] for c in present) / len(present) if present else float("nan")


def build(samples, fold_pred, readouts, markers):
    """Per-sample records: drug, cell line, pairs (prediction, delta-p) with the drug's markers left out, tallies."""
    out = []
    for s in samples:
        drop = markers.get(s["drug"], set())
        pairs = [(fold_pred[(s["drug"], a)], v) for a, v in s["values"].items()
                 if a in readouts and a not in drop and (s["drug"], a) in fold_pred]
        if not pairs:
            continue
        tot, hit, gated = tally(pairs)
        out.append({"sample": s, "drug": s["drug"], "line": s["cell_line"], "pairs": pairs,
                    "tot": tot, "hit": hit, "gated": gated})
    return out


def gated_at(rec, a, theta):
    return theta is not None and a is not None and a < theta


def score(recs, acts, gate):
    """Balanced accuracy over records; gate(i) says whether record i is gated."""
    tot, hit = collections.Counter(), collections.Counter()
    for i, r in enumerate(recs):
        tot.update(r["tot"])
        hit.update(r["gated"] if gate(i) else r["hit"])
    return balanced(tot, hit)


def choose(recs, acts, idx):
    """The theta with the highest balanced accuracy over records idx; ties to the identity, then the larger theta."""
    sub = [recs[i] for i in idx]
    suba = [acts[i] for i in idx]
    def rank(t):
        return (round(score(sub, suba, lambda j: gated_at(sub[j], suba[j], t)), 12), t is None, t or 0)
    return max(THETAS, key=rank)


def cross_validate(recs, acts, seed=SEED):
    """(gated flags per record under CV, theta chosen per fold, full-data theta)."""
    folds = score_130.assign_folds({r["line"] for r in recs}, FOLDS, seed)
    gated, chosen = [False] * len(recs), []
    for f in range(FOLDS):
        t = choose(recs, acts, [i for i, r in enumerate(recs) if folds[r["line"]] != f])
        chosen.append(t)
        for i, r in enumerate(recs):
            if folds[r["line"]] == f:
                gated[i] = gated_at(r, acts[i], t)
    return gated, chosen, choose(recs, acts, range(len(recs)))


def per_drug_tallies(recs, gated):
    """{drug: (tot, hit under gate, hit under identity)}."""
    out = collections.defaultdict(lambda: (collections.Counter(), collections.Counter(), collections.Counter()))
    for r, g in zip(recs, gated):
        tot, hg, hi = out[r["drug"]]
        tot.update(r["tot"])
        hg.update(r["gated"] if g else r["hit"])
        hi.update(r["hit"])
    return out


def bootstrap_ba_gain(by_drug, n=N_BOOT, seed=SEED):
    drugs = sorted(by_drug)
    def gain(sample):
        tot, hg, hi = collections.Counter(), collections.Counter(), collections.Counter()
        for d in sample:
            t, g, i = by_drug[d]
            tot.update(t)
            hg.update(g)
            hi.update(i)
        return balanced(tot, hg) - balanced(tot, hi)
    rng = random.Random(seed)
    draws = sorted(gain([rng.choice(drugs) for _ in drugs]) for _ in range(n))
    return gain(drugs), draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


def random_null(recs, acts, gated, n=N_NULL, seed=SEED):
    """Gains of gating, per drug, as many randomly chosen eligible samples as the CV gate gated."""
    base = score(recs, acts, lambda i: False)
    k = collections.Counter(r["drug"] for r, g in zip(recs, gated) if g)
    pool = collections.defaultdict(list)
    for i, (r, a) in enumerate(zip(recs, acts)):
        if a is not None:
            pool[r["drug"]].append(i)
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        chosen = set()
        for d, kd in k.items():
            chosen |= set(rng.sample(pool[d], kd))
        out.append(score(recs, acts, lambda i: i in chosen) - base)
    return sorted(out)


def correlations(recs, gated):
    """({drug: mean per-sample correlation, gated}, {drug: identity}); a gated sample predicts a constant, so 0."""
    g, i = collections.defaultdict(list), collections.defaultdict(list)
    for r, flag in zip(recs, gated):
        c = score_121.pearson([x for x, _ in r["pairs"]], [y for _, y in r["pairs"]])
        if c is None:
            continue
        i[r["drug"]].append(c)
        g[r["drug"]].append(0.0 if flag else c)
    mean = lambda by: {d: sum(v) / len(v) for d, v in by.items()}  # noqa: E731
    return mean(g), mean(i)


def premise(recs, acts):
    """P3: share of changed (non-marker) pairs in eligible samples with A >= 0.5 and with A < 0.5."""
    acc = {True: [0, 0], False: [0, 0]}
    for r, a in zip(recs, acts):
        if a is None:
            continue
        t = r["tot"]
        acc[a >= 0.5][0] += t[-1] + t[1]
        acc[a >= 0.5][1] += sum(t.values())
    return {("high" if k else "low"): {"changed": v[0], "pairs": v[1], "share": v[0] / v[1] if v[1] else None}
            for k, v in acc.items()}


def evaluate(recs, acts):
    gated, chosen, full = cross_validate(recs, acts)
    by_drug = per_drug_tallies(recs, gated)
    gain = bootstrap_ba_gain(by_drug)
    null = random_null(recs, acts, gated)
    cg, ci = correlations(recs, gated)
    drugs = sorted(set(cg) & set(ci))
    corr = score_130.bootstrap_gain(cg, ci, drugs, N_BOOT, SEED)
    key = lambda t: "identity" if t is None else t  # noqa: E731
    return {"identity_ba": score(recs, acts, lambda i: False),
            "cv_ba": score(recs, acts, lambda i: gated[i]),
            "cv_gain": gain[0], "cv_gain_ci95": list(gain[1:]),
            "null_p95": null[int(0.95 * len(null)) - 1], "null_mean": sum(null) / len(null),
            "beats_null": gain[0] > null[int(0.95 * len(null)) - 1],
            "corr_gain": corr[0], "corr_gain_ci95": list(corr[1:]),
            "theta_per_fold": [key(t) for t in chosen], "theta_full": key(full),
            "gated_samples": sum(gated), "gated_drugs": len({r["drug"] for r, g in zip(recs, gated) if g}),
            "by_theta_full_data": {str(key(t)): score(recs, acts, lambda i, t=t: gated_at(recs[i], acts[i], t))
                                   for t in THETAS}}


def main():
    unpacked, split_tsv, readouts_tsv, markers_tsv, *logs = sys.argv[1:]
    kinds = {r["readout"]: r["kind"] for r in csv.DictReader(
        (ln for ln in open(readouts_tsv) if not ln.startswith("#")), delimiter="\t")}
    readouts = {a for a, k in kinds.items() if k != "exclude"}
    phospho = {a for a, k in kinds.items() if k == "phospho"}
    markers = read_markers(markers_tsv)
    samples = cppa_data_121.load(unpacked, split_tsv, "train")
    fold_pred = score_121.predictions(logs)["fold"]
    rank_table = ranks(score_130.control_baselines(unpacked, split_tsv))

    recs = build(samples, fold_pred, readouts, markers)
    acts = [activity(r["sample"], markers[r["drug"]], rank_table) if r["drug"] in markers else None for r in recs]
    general = [activity(r["sample"], phospho, rank_table) for r in recs]
    with_markers = [a for r, a in zip(recs, acts) if r["drug"] in markers]
    recs_all = build(samples, fold_pred, readouts, {})

    print(json.dumps({
        "portion": "train", "samples_scored": len(recs),
        "marker_drugs_scored": len({r["drug"] for r in recs if r["drug"] in markers}),
        "p1_identity_ba_markers_out": score(recs, acts, lambda i: False),
        "p1_identity_ba_markers_in": score(recs_all, [None] * len(recs_all), lambda i: False),
        "p2_marker_drug_samples": len(with_markers),
        "p2_eligible": sum(a is not None for a in with_markers),
        "p3_change_share": premise(recs, acts),
        "marker_gate": evaluate(recs, acts),
        "general_gate": evaluate(recs, general),
    }, indent=1, default=str))


if __name__ == "__main__":
    main()
