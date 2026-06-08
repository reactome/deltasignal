#!/usr/bin/env python3
"""Produce a publication-style summary table for one or more per-case dumps.

Combines case_dump_overview, classify_failures, and threshold_sweep into
one report. Designed to be copy-pasted into a manuscript appendix or run
against a future model to compare apples-to-apples.

Output sections per dump:
  1. Overall accuracy + confusion matrix
  2. Per-pathway accuracy
  3. Failure breakdown by root-cause category (requires gene cache)
  4. Continuous-output distribution per expected class

Multiple dumps print side-by-side in a delta-style table.

Usage
-----
  python report.py /tmp/ds_or_max.tsv
  python report.py /tmp/ds_or_max.tsv /tmp/ds_or_mean.tsv \
      --labels or_max,or_mean
"""
import csv, sys
from collections import defaultdict
from pathlib import Path

CLS = {"0": "DOWN", "1": "NORM", "2": "UP"}


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def confusion(rows):
    cm = defaultdict(lambda: defaultdict(int))
    for r in rows:
        cm[r["expected"]][r["predicted"]] += 1
    return cm


def per_pathway(rows):
    pw = defaultdict(lambda: {"c": 0, "t": 0})
    for r in rows:
        pw[r["pathway"]]["t"] += 1
        if r["predicted"] == r["expected"]:
            pw[r["pathway"]]["c"] += 1
    return pw


def maybe_classify(rows):
    """Attempt failure classification; return None if gene cache or networks
    aren't available."""
    try:
        sys.path.insert(0, str(Path(__file__).parent))
        from _common import (Networks, load_gene_stids, bfs_upstream,
                             resolve_gene_uuids, resolve_readout_uuids)
        gene_stids = load_gene_stids()
        nets = Networks()
    except (ImportError, SystemExit) as e:
        return None

    cats = defaultdict(int)
    for r in rows:
        if r["predicted"] == r["expected"]:
            continue
        inc, _, _, rev, proxies = nets.get(r["pathway"])
        gene_uuids = resolve_gene_uuids(r["gene"], gene_stids, rev)
        ko_uuids = resolve_readout_uuids(r["key_output"], rev, proxies)
        if not gene_uuids: cats["GENE_NOT_IN_PATHWAY"] += 1; continue
        if not ko_uuids:   cats["READOUT_NOT_IN_PATHWAY"] += 1; continue
        upstream = bfs_upstream(inc, set(ko_uuids), max_depth=30)
        if not (set(gene_uuids) & upstream):
            cats["NO_PATH_GENE_TO_READOUT"] += 1
        else:
            e, p = r["expected"], r["predicted"]
            if (e, p) in (("0", "2"), ("2", "0")):
                cats["PATH_EXISTS_DIRECTION_FLIPPED"] += 1
            elif e == "1":
                cats["PATH_EXISTS_FALSE_POSITIVE"] += 1
            elif p == "1":
                cats["PATH_EXISTS_SIGNAL_LOST"] += 1
            else:
                cats[f"PATH_EXISTS_OTHER"] += 1
    return cats


def dist_by_expected(rows):
    by = defaultdict(list)
    for r in rows:
        by[r["expected"]].append(float(r["pred_ui"]))
    out = {}
    for e, vs in by.items():
        vs.sort()
        if not vs: continue
        n = len(vs)
        out[e] = {"n": n, "min": vs[0], "p25": vs[n//4], "p50": vs[n//2],
                  "p75": vs[3*n//4], "max": vs[-1]}
    return out


def print_one(rows, label):
    n = len(rows)
    c = sum(1 for r in rows if r["predicted"] == r["expected"])
    print(f"\n=================== {label} ===================")
    print(f"Accuracy: {c}/{n} = {c*100/n:.2f}%")

    cm = confusion(rows)
    print("\nConfusion (rows=expected, cols=predicted):")
    print(f"  {'':6s} {'DOWN':>5s} {'NORM':>5s} {'UP':>5s}")
    for e in ("0", "1", "2"):
        print(f"  {CLS[e]:6s} " + " ".join(f"{cm[e][p]:>5d}" for p in ("0","1","2")))

    pw = per_pathway(rows)
    print("\nPer pathway:")
    for p in sorted(pw, key=lambda p: -pw[p]["t"]):
        d = pw[p]
        print(f"  {p[:38]:40s} {d['c']:3d}/{d['t']:3d} "
              f"({d['c']*100/d['t']:5.1f}%)")

    cats = maybe_classify(rows)
    if cats is not None:
        print("\nFailure root-cause breakdown:")
        for k, v in sorted(cats.items(), key=lambda x: -x[1]):
            print(f"  {k:35s} {v}")

    dist = dist_by_expected(rows)
    print("\npred_ui distribution by expected class:")
    for e in ("0", "1", "2"):
        if e not in dist: continue
        d = dist[e]
        print(f"  exp={CLS[e]:4s} (n={d['n']:3d})  "
              f"min={d['min']:6.2f} p25={d['p25']:6.2f} p50={d['p50']:6.2f} "
              f"p75={d['p75']:6.2f} max={d['max']:6.2f}")


def print_delta(rows_a, rows_b, label_a, label_b):
    na, nb = len(rows_a), len(rows_b)
    ca = sum(1 for r in rows_a if r["predicted"] == r["expected"])
    cb = sum(1 for r in rows_b if r["predicted"] == r["expected"])
    print(f"\n=================== {label_a} vs {label_b} ===================")
    print(f"{label_a}: {ca}/{na} = {ca*100/na:.2f}%")
    print(f"{label_b}: {cb}/{nb} = {cb*100/nb:.2f}%")
    print(f"Δ      : {cb - ca:+d}")
    a_pw = per_pathway(rows_a); b_pw = per_pathway(rows_b)
    pws = sorted(set(a_pw) | set(b_pw))
    print("\nPer-pathway delta:")
    for p in pws:
        a = a_pw.get(p, {"c":0,"t":0}); b = b_pw.get(p, {"c":0,"t":0})
        ap = a["c"]*100/a["t"] if a["t"] else 0
        bp = b["c"]*100/b["t"] if b["t"] else 0
        print(f"  {p[:38]:40s} {a['c']:3d}/{a['t']:3d} ({ap:5.1f}%) → "
              f"{b['c']:3d}/{b['t']:3d} ({bp:5.1f}%)  Δ={b['c']-a['c']:+d}")


def main(paths, labels=None):
    if not labels:
        labels = [Path(p).stem for p in paths]
    dumps = [load(p) for p in paths]
    for rows, lab in zip(dumps, labels):
        print_one(rows, lab)
    # Pairwise deltas
    for i in range(len(dumps)):
        for j in range(i + 1, len(dumps)):
            print_delta(dumps[i], dumps[j], labels[i], labels[j])


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__, file=sys.stderr); sys.exit(1)
    labels = None
    if "--labels" in args:
        idx = args.index("--labels")
        labels = args[idx + 1].split(",")
        args = args[:idx] + args[idx + 2:]
    main(args, labels)
