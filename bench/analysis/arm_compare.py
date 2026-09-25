#!/usr/bin/env python3
"""Paired comparison of two case dumps on the cases valid in BOTH.

Reports, per ground-truth axis and split (held-out / tuning / all): accuracy and
macro-F1 of each arm on the shared valid set, fixed / broke / net, exact
McNemar p, how many pathways moved (and in which direction), how many distinct
perturbations the discordant cases span, and the valid-set sizes (a case that
becomes invalid in one arm is not silently scored as a loss or a gain).

Usage: arm_compare.py BASE_DIR ARM_DIR   (each holding {curator,experimental}_cases.tsv)
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import mcnemar_exact  # noqa: E402
from holdout_report import TUNING_PATHWAYS  # noqa: E402


def load(path: Path) -> dict:
    with open(path, newline="") as fh:
        return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
                for r in csv.DictReader(fh, delimiter="\t")}


def macro_f1(pairs) -> float:
    f1s = []
    for c in ("0", "1", "2"):
        tp = sum(1 for p, e in pairs if p == c and e == c)
        fp = sum(1 for p, e in pairs if p == c and e != c)
        fn = sum(1 for p, e in pairs if p != c and e == c)
        f1s.append(2 * tp / (2 * tp + fp + fn) if tp else 0.0)
    return sum(f1s) / 3


def compare(a: dict, b: dict, keep) -> dict:
    va = {k for k, r in a.items() if r["valid"] == "1" and keep(k)}
    vb = {k for k, r in b.items() if r["valid"] == "1" and keep(k)}
    shared = va & vb
    per = defaultdict(lambda: [0, 0])
    genes = set()
    for k in shared:
        ra = a[k]["predicted"] == a[k]["expected"]
        rb = b[k]["predicted"] == b[k]["expected"]
        if ra != rb:
            per[k[0]][0 if rb else 1] += 1
            genes.add((k[0], k[1], k[2]))
    fx = sum(v[0] for v in per.values())
    br = sum(v[1] for v in per.values())
    acc = lambda d: sum(d[k]["predicted"] == d[k]["expected"] for k in shared) / max(1, len(shared))
    mf = lambda d: macro_f1([(d[k]["predicted"], d[k]["expected"]) for k in shared])
    return {"valid_base": len(va), "valid_arm": len(vb), "shared": len(shared),
            "acc_base": acc(a), "acc_arm": acc(b), "mf1_base": mf(a), "mf1_arm": mf(b),
            "fixed": fx, "broke": br, "p": mcnemar_exact(fx, br),
            "moved": len(per), "up": sum(v[0] > v[1] for v in per.values()),
            "down": sum(v[0] < v[1] for v in per.values()), "perturbations": len(genes),
            "per_pathway": {p: v[0] - v[1] for p, v in per.items()}}


def main() -> int:
    base, arm = Path(sys.argv[1]), Path(sys.argv[2])
    for ax in ("curator", "experimental"):
        a, b = load(base / f"{ax}_cases.tsv"), load(arm / f"{ax}_cases.tsv")
        splits = (("held-out", lambda k: k[0] not in TUNING_PATHWAYS),
                  ("tuning", lambda k: k[0] in TUNING_PATHWAYS), ("all", lambda k: True))
        for name, keep in splits:
            r = compare(a, b, keep)
            print(f"{ax:12} {name:8} valid {r['valid_base']:>5}->{r['valid_arm']:<5} shared {r['shared']:>5}  "
                  f"acc {r['acc_base']:.4f}->{r['acc_arm']:.4f}  mF1 {r['mf1_base']:.4f}->{r['mf1_arm']:.4f}  "
                  f"net {r['fixed'] - r['broke']:+5d} ({r['fixed']}/{r['broke']}) p={r['p']:.2g}  "
                  f"pathways {r['moved']} (+{r['up']}/-{r['down']}) perturbations {r['perturbations']}")
            if ax == "curator" and name == "all" and r["per_pathway"]:
                pp = sorted(r["per_pathway"].items(), key=lambda kv: kv[1])
                print("   worst:", ", ".join(f"{p[:28]} {n:+d}" for p, n in pp[:4]))
                print("   best: ", ", ".join(f"{p[:28]} {n:+d}" for p, n in pp[-4:][::-1]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
