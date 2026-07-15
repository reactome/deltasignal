#!/usr/bin/env python3
"""Sweep up/down classification cutoffs to MAXIMIZE macro-F1 on a per-case dump.

Unlike threshold_sweep.py (which maximizes accuracy and therefore tends to
over-predict NORMAL), this optimizes macro-F1 across the DOWN/NORMAL/UP classes,
so the chosen cutoffs actually detect up/down-regulation rather than defaulting
to "unchanged". One GLOBAL (DOWN, UP) pair is chosen (applied to every pathway),
so it transfers to held-out pathways.

Reads pred_ui (continuous, UI 0-100 scale, 1=baseline) + expected class
(0=down,1=normal,2=up) from a --dump-cases TSV.

Usage:
  python f1_threshold_sweep.py dump.tsv [--down lo:hi:step] [--up lo:hi:step] [--metric macro3|macroUD]
"""
import csv, sys
from collections import defaultdict


def classify(ui, down, up):
    return "0" if ui < down else ("2" if ui >= up else "1")


def frange(lo, hi, step):
    out, v = [], lo
    while v <= hi + 1e-9:
        out.append(round(v, 6)); v += step
    return out


def f1_scores(rows, down, up):
    """Return per-class F1 dict for classes 0/1/2."""
    tp = defaultdict(int); fp = defaultdict(int); fn = defaultdict(int)
    for r in rows:
        pred = classify(r["pred_ui"], down, up)
        exp = r["expected"]
        for c in ("0", "1", "2"):
            if pred == c and exp == c: tp[c] += 1
            elif pred == c and exp != c: fp[c] += 1
            elif pred != c and exp == c: fn[c] += 1
    f1 = {}
    for c in ("0", "1", "2"):
        p = tp[c] / (tp[c] + fp[c]) if (tp[c] + fp[c]) else 0.0
        rec = tp[c] / (tp[c] + fn[c]) if (tp[c] + fn[c]) else 0.0
        f1[c] = 2 * p * rec / (p + rec) if (p + rec) else 0.0
    return f1


def main():
    args = sys.argv[1:]
    path = args[0]
    down_spec = up_spec = None
    metric = "macroUD"  # macro-F1 over DOWN+UP (the regulation calls)
    for i, a in enumerate(args):
        if a == "--down": down_spec = args[i + 1]
        if a == "--up": up_spec = args[i + 1]
        if a == "--metric": metric = args[i + 1]
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    # only cases with a valid prediction + expected class
    rows = [r for r in rows if r.get("pred_ui") not in (None, "") and r.get("expected") in ("0", "1", "2")]
    for r in rows: r["pred_ui"] = float(r["pred_ui"])

    downs = frange(*[float(x) for x in down_spec.split(":")]) if down_spec else frange(0.50, 0.99, 0.01)
    ups = frange(*[float(x) for x in up_spec.split(":")]) if up_spec else frange(1.01, 2.50, 0.01)

    def score(f1):
        return (f1["0"] + f1["2"]) / 2 if metric == "macroUD" else (f1["0"] + f1["1"] + f1["2"]) / 3

    dist = defaultdict(int)
    for r in rows: dist[r["expected"]] += 1
    print(f"N={len(rows)}  class dist (down/normal/up): {dist['0']}/{dist['1']}/{dist['2']}  metric={metric}")

    best = None
    for d in downs:
        for u in ups:
            if u <= d: continue
            f1 = f1_scores(rows, d, u)
            s = score(f1)
            if best is None or s > best[0]:
                best = (s, d, u, f1)
    s, d, u, f1 = best
    acc = sum(1 for r in rows if classify(r["pred_ui"], d, u) == r["expected"]) / len(rows)
    print(f"\nBEST macro-F1 ({metric}) = {s:.4f}  at DOWN={d}  UP={u}")
    print(f"  per-class F1: down={f1['0']:.3f} normal={f1['1']:.3f} up={f1['2']:.3f}")
    print(f"  accuracy at this cutoff: {acc*100:.2f}%")
    # per-pathway accuracy at the chosen cutoff
    bypw = defaultdict(lambda: [0, 0])
    for r in rows:
        ok = classify(r["pred_ui"], d, u) == r["expected"]
        bypw[r["pathway"]][0] += ok; bypw[r["pathway"]][1] += 1
    print("  per-pathway acc @ best cutoff:")
    for pw in sorted(bypw):
        c, n = bypw[pw]; print(f"    {pw:36} {c}/{n} ({100*c/n:.1f}%)")


if __name__ == "__main__":
    main()
