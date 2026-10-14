"""specs/057: edge-type counts per pathway, canonical build vs a variant build.

Under LNG_UB_SUBSTRATE_CHANGED only depletion edges may differ (fewer), and
the variant's own census (ub_census.py) must find none left.

    python3 bench/analysis/057/edge_diff.py BUILD_A BUILD_B
"""
import collections
import csv
import glob
import os
import sys


def counts(build):
    out = {}
    for d in sorted(glob.glob(os.path.join(build, "R-HSA-*"))):
        c = collections.Counter()
        for f in ("logic_network.csv", "boundary_edges.csv"):
            p = os.path.join(d, f)
            if os.path.exists(p):
                for r in csv.DictReader(open(p)):
                    c[r["edge_type"]] += 1
        out[os.path.basename(d)] = c
    return out


def main(a, b):
    ca, cb = counts(a), counts(b)
    total = collections.Counter()
    for pid in sorted(set(ca) | set(cb)):
        x, y = ca.get(pid, collections.Counter()), cb.get(pid, collections.Counter())
        diff = {t: y[t] - x[t] for t in set(x) | set(y) if y[t] != x[t]}
        if diff:
            print(pid, diff)
            total.update(diff)
    print("total", dict(total))


if __name__ == "__main__":
    main(*sys.argv[1:3])
