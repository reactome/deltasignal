#!/usr/bin/env python3
"""specs/161: each korkut drug's measured engagement, the pin value the arm uses.

From the drug's OWN single-drug conditions (both doses), its direct markers
(specs/138 MARKERS, fixed before this spec): m = 2^(mean signed-consistent
log2 fold over markers and doses). Inhibitor markers (expected sign -):
e = min(1, m) — the target's remaining activity read from its direct
substrate. Positive-sign markers (Nutlin: MDM2 blocked -> p53/p21 up):
e = min(1, 1/m). No marker on the array (HDAC, PKC): e = 0, the recorded full
knockout. Registered in spec.md; nothing fitted.
"""
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "138"))
import korkut_bench as k  # noqa: E402
from engagement import MARKERS  # noqa: E402


def engagement():
    names = [l.strip() for l in open(k.DATA / "node_Index.csv")]
    vals = {}
    for (_, ds), row in zip(k.conditions(), csv.reader(open(k.DATA / "expr.csv"))):
        if len(ds) == 1:
            vals.setdefault(ds[0], []).append(dict(zip(names, map(float, row))))
    out = {}
    for d, marks in MARKERS.items():
        if not marks:
            out[d] = (0.0, "no marker: full knockout (recorded protocol)")
            continue
        xs = [v[m] for m, _ in marks for v in vals[d]]
        sign = marks[0][1]
        m = 2.0 ** (sum(xs) / len(xs))
        e = min(1.0, m) if sign < 0 else min(1.0, 1.0 / m)
        out[d] = (e, f"markers {','.join(mk for mk, _ in marks)} mean log2 {sum(xs)/len(xs):+.3f}")
    return out


def marker_names():
    return {d: {m for m, _ in marks} for d, marks in MARKERS.items()}


if __name__ == "__main__":
    for d, (e, why) in sorted(engagement().items()):
        print(f"{d}\t{e:.4f}\t{why}")
