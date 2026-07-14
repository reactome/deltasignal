#!/usr/bin/env python3
"""Sweep classification thresholds against a per-case dump.

Reads pred_ui (continuous, 0–100 scale) and re-classifies each case for
a range of (DOWN, UP) cutoffs. Finds the single GLOBAL (DOWN, UP) pair
that maximizes accuracy across all pathways — and reports the per-pathway
accuracy at that pair. Same cutoff is applied to every pathway so the
result extends cleanly to held-out pathways.

For bimodal output distributions (e.g., multiplicative AND aggregation),
the sweep typically finds many cutoffs giving the same accuracy — that
indicates the threshold isn't the bottleneck. For continuous outputs
(e.g., hill_log), a clear optimum emerges.

Usage
-----
  python threshold_sweep.py /tmp/ds_best_cases.tsv
  python threshold_sweep.py /tmp/ds_best_cases.tsv --down 0.5:1.0:0.025 --up 1.05:2.5:0.025
"""
import csv, sys
from collections import defaultdict


def classify(ui, down_cut, up_cut):
    if ui < down_cut:  return "0"
    if ui >= up_cut:   return "2"
    return "1"


def parse_range(spec, default):
    if spec is None:
        return default
    lo, hi, step = (float(x) for x in spec.split(":"))
    out = []
    v = lo
    while v <= hi + 1e-9:
        out.append(round(v, 6))
        v += step
    return out


def main(path, down_spec=None, up_spec=None):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rows: r["pred_ui"] = float(r["pred_ui"])
    n = len(rows)
    baseline_acc = sum(1 for r in rows if r["predicted"] == r["expected"]) * 100 / n
    print(f"Baseline (as-dumped): {baseline_acc:.2f}%   N={n}")

    downs = parse_range(down_spec, [round(0.5 + 0.025*i, 4) for i in range(21)])  # 0.5 → 1.0
    ups   = parse_range(up_spec,   [round(1.05 + 0.025*i, 4) for i in range(59)])  # 1.05 → 2.5

    best = (0, None, None, None)
    grid = []
    for d in downs:
        for u in ups:
            if u <= d: continue
            c = sum(1 for r in rows if classify(r["pred_ui"], d, u) == r["expected"])
            grid.append((c, d, u))
            if c > best[0]:
                best = (c, d, u, None)

    print(f"\nBest GLOBAL cutoffs: DOWN<{best[1]:.3f}, UP≥{best[2]:.3f}")
    print(f"  → {best[0]}/{n} = {best[0]*100/n:.2f}%  "
          f"(Δ from dumped: {best[0]*100/n - baseline_acc:+.2f}pp)")

    # Per-pathway accuracy at best cutoffs (for sanity-check)
    print("\n=== Per-pathway accuracy at GLOBAL best cutoffs ===")
    pw = defaultdict(lambda: [0, 0])
    for r in rows:
        c = classify(r["pred_ui"], best[1], best[2])
        pw[r["pathway"]][0] += int(c == r["expected"])
        pw[r["pathway"]][1] += 1
    for p in sorted(pw, key=lambda p: -pw[p][1]):
        c, t = pw[p]
        print(f"  {p[:38]:40s} {c:3d}/{t:3d} ({c*100/t:5.1f}%)")

    # Show top-10 cutoff pairs (often a wide tie band → bimodal outputs)
    grid.sort(reverse=True)
    print("\n=== Top 10 (DOWN, UP) cutoff pairs ===")
    for c, d, u in grid[:10]:
        print(f"  DOWN<{d:.3f}  UP≥{u:.3f}  → {c:3d}/{n} = {c*100/n:.2f}%")


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__, file=sys.stderr); sys.exit(1)
    path = args[0]
    ds = us = None
    if "--down" in args: ds = args[args.index("--down") + 1]
    if "--up" in args:   us = args[args.index("--up") + 1]
    main(path, ds, us)
