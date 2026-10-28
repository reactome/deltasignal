#!/usr/bin/env python3
"""specs/115: where the frozen model's errors on a case log fall, overall and by hop distance from the pinned targets.

Classes over scored cases whose measured call is not `conflict`:
- changed / no route: measured DOWN or UP, reachability NOT (no route from a pinned target);
- changed / right; changed / model NOT; changed / wrong direction (a route exists);
- unchanged / right; unchanged / false change call.
By hop band (needs a hops table from korkut_hops.py): cases, measured-changed share, model change calls, the share of
those right, and direction agreement where both measured and model call a change.

Usage: error_anatomy_115.py korkut|decryptm CASES_TSV [HOPS_TSV]
"""
import collections
import csv
import json
import sys

KEY = {"korkut": ("pathway_id", "drug_set", "readout")}
BANDS = ((1, 5), (6, 10), (11, 20), (21, 40), (41, 10 ** 9))


def kind(r):
    m, p, reach = r["measured"], r["predicted"], r["reach_predicted"]
    if m == "not":
        return "unchanged / right" if p == "not" else "unchanged / false change call"
    if reach == "not":
        return "changed / no route"
    if p == m:
        return "changed / right"
    return "changed / model NOT" if p == "not" else "changed / wrong direction"


def band(h):
    if h in ("", None):
        return "no route"
    h = int(h)
    return next(f"{lo}-{hi}" if hi < 10 ** 9 else f">{lo - 1}" for lo, hi in BANDS if lo <= h <= hi)


def anatomy(rows, hops=None, key=None):
    rows = [r for r in rows if r["status"] == "scored" and r["measured"] != "conflict"]
    out = {"classes": dict(collections.Counter(map(kind, rows)))}
    if hops is not None:
        groups = collections.defaultdict(list)
        for r in rows:
            groups[band(hops.get(tuple(r[k] for k in key), ""))].append(r)
        out["by_hops"] = {}
        for b, g in groups.items():
            calls = [r for r in g if r["predicted"] != "not"]
            both = [r for r in calls if r["measured"] != "not"]
            out["by_hops"][b] = {"cases": len(g), "changed": sum(r["measured"] != "not" for r in g),
                                 "model_calls": len(calls), "model_calls_right": sum(r["predicted"] == r["measured"] for r in calls),
                                 "both_changed": len(both), "direction_agrees": sum(r["predicted"] == r["measured"] for r in both)}
    return out


def main():
    kind_, path = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    hops = None
    if len(sys.argv) > 3:
        hops = {tuple(r[k] for k in KEY[kind_]): r["hops"] for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")}
    print(json.dumps(anatomy(rows, hops, KEY.get(kind_)), indent=1))


if __name__ == "__main__":
    main()
