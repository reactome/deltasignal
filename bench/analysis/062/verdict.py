#!/usr/bin/env python3
"""specs/062 step 2: the registered verdict from two korkut_cases.tsv files.

Usage: verdict.py CONTROL_CASES ARM_CASES [measured|measured_secondary]
"""
import collections
import csv
import math
import sys


def mcnemar_p(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def macro_f1(pairs):
    f1 = []
    for lab in ("down", "not", "up"):
        tp = sum(p == lab and m == lab for p, m in pairs)
        fp = sum(p == lab and m != lab for p, m in pairs)
        fn = sum(p != lab and m == lab for p, m in pairs)
        f1.append(0.0 if tp == 0 else 2 * tp / (2 * tp + fp + fn))
    return sum(f1) / 3


def main():
    ctrl, arm = sys.argv[1:3]
    col = sys.argv[3] if len(sys.argv) > 3 else "measured"
    load = lambda p: {(r["pathway_id"], r["drug_set"], r["readout"]): r
                      for r in csv.DictReader(open(p), delimiter="\t")}
    a, b = load(ctrl), load(arm)
    keys = sorted(k for k in a.keys() & b.keys()
                  if a[k]["status"] == b[k]["status"] == "scored" and a[k][col] != "conflict")
    one = sum(1 for k in a.keys() | b.keys()
              if (a.get(k, {}).get("status") == "scored") != (b.get(k, {}).get("status") == "scored"))
    print(f"measured column: {col}; cases scored in both arms (conflicts left out): {len(keys)}; "
          f"scored in one arm only (left out): {one}")
    for name, d in (("control", a), ("arm", b)):
        pairs = [(d[k]["predicted"], d[k][col]) for k in keys]
        acc = sum(p == m for p, m in pairs) / len(pairs)
        print(f"  {name}: accuracy {acc:.4f}, macro-F1 {macro_f1(pairs):.4f}")
    down = [("down", a[k][col]) for k in keys]
    print(f"  always-DOWN: accuracy {sum(p == m for p, m in down) / len(down):.4f}, macro-F1 {macro_f1(down):.4f}")
    for name, d in (("control", a), ("arm", b)):
        if all(d[k].get("reach_predicted") for k in keys):
            rp = [(d[k]["reach_predicted"], d[k][col]) for k in keys]
            print(f"  reachability ({name} pins): accuracy {sum(p == m for p, m in rp) / len(rp):.4f}, "
                  f"macro-F1 {macro_f1(rp):.4f}")
    nots = [("not", a[k][col]) for k in keys]
    print(f"  always-NOT: accuracy {sum(p == m for p, m in nots) / len(nots):.4f}, macro-F1 {macro_f1(nots):.4f}")
    print(f"  measured: {dict(collections.Counter(a[k][col] for k in keys))}")
    disc = [k for k in keys if a[k]["predicted"] != b[k]["predicted"]]
    fixed = [k for k in disc if b[k]["predicted"] == b[k][col] and a[k]["predicted"] != a[k][col]]
    broke = [k for k in disc if a[k]["predicted"] == a[k][col] and b[k]["predicted"] != b[k][col]]
    p = mcnemar_p(len(fixed), len(broke))
    net = len(fixed) - len(broke)
    verdict = ("confirms" if net > 0 and p < 0.05 else "contradicts" if net < 0 and p < 0.05 else "inconclusive")
    print(f"discordant predictions {len(disc)}: arm right / control wrong {len(fixed)}, "
          f"control right / arm wrong {len(broke)}, both wrong {len(disc) - len(fixed) - len(broke)}; "
          f"net {net:+d}, McNemar p {p:.4g} -> VERDICT: {verdict}")
    # clustered by drug set: one vote per drug set (its net sign)
    by = collections.Counter()
    for k in fixed:
        by[k[1]] += 1
    for k in broke:
        by[k[1]] -= 1
    pos, neg = sum(v > 0 for v in by.values()), sum(v < 0 for v in by.values())
    print(f"clustered by drug set: {len(by)} drug sets with a fix or break; net positive {pos}, "
          f"negative {neg}, zero {len(by) - pos - neg}; sign-test p {mcnemar_p(pos, neg):.4g}")
    units = collections.Counter()
    for k in disc:
        units[(k[0], k[2], a[k]["predicted"], b[k]["predicted"], a[k][col])] += 1
    print("by (pathway, readout, control -> arm, measured):")
    for u, n in sorted(units.items()):
        print(f"  {n:3d} {u}")


if __name__ == "__main__":
    main()
