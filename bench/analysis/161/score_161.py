#!/usr/bin/env python3
"""specs/161: score engagement-scaled pins against the full-knockout canonical
run on the same rows.

Truth: the recorded 98 korkut TSV's `measured`. Baseline: the specs/157
canonical run's korkut calls (same build, same /api/readouts valuation).
Excluded from BOTH arms: any readout that is a direct marker of a drug in the
row's drug set (the markers set the pin values). Gate: balanced accuracy,
95% bootstrap over drug sets (10k, seed 161), interval excluding zero.

Usage: score_161.py CANONICAL_CALLS[.gz] ARM_CALLS KORKUT_CASES_TSV
"""
import csv
import gzip
import random
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent))
import dataset_roles  # noqa: E402
from engagement_161 import engagement, marker_names  # noqa: E402

N_BOOT = 10_000
SEED = 161
RAF = "R-HSA-5673001"


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_calls(path):
    out = {}
    with opener(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["axis"] == "korkut" and r["call"]:
                out[(r["pathway"], r["case_a"], r["readout"])] = (r["call"], float(r["fold"]))
    return out


def ba(rows, col):
    per = []
    for m in ("down", "not", "up"):
        sub = [r for r in rows if r["measured"] == m]
        if sub:
            per.append(sum(r[col] == m for r in sub) / len(sub))
    return sum(per) / len(per)


def scorecard(rows, col, fold_col):
    reach = [r for r in rows if r["base_fold"] != 1.0]
    chg = [r for r in reach if r["measured"] != "not"]
    nc = [r for r in reach if r["measured"] == "not"]
    called = [r for r in reach if r[col] != "not"]
    tp = sum(r[col] == r["measured"] for r in chg)
    import math
    def s(r):
        f = r[fold_col]
        return abs(math.log2(f)) if f > 0 else 99.0
    pos = [s(r) for r in chg]; neg = [s(r) for r in nc]
    auc = sum((1 if p > q else 0.5 if p == q else 0) for p in pos for q in neg) / (len(pos) * len(neg))
    return (tp / len(chg), sum(r[col] == "not" for r in nc) / len(nc),
            tp / len(called) if called else float("nan"), auc, len(reach))


def main():
    dataset_roles.require_role("korkut2015", "fit")
    base = load_calls(sys.argv[1])
    arm = load_calls(sys.argv[2])
    marks = marker_names()
    eng = engagement()
    rows, excluded = [], 0
    for r in csv.DictReader(open(sys.argv[3]), delimiter="\t"):
        if r["status"] != "scored":
            continue
        k = (r["pathway_id"], r["drug_set"], r["readout"])
        if k not in base or k not in arm:
            continue
        drugs = r["drug_set"].split("+")
        if any(r["readout"] in marks.get(d, set()) for d in drugs):
            excluded += 1
            continue
        engaged_all = all(eng[d][0] <= 0.85 for d in drugs)  # anatomy only
        rows.append({"drug_set": r["drug_set"], "pathway": r["pathway_id"],
                     "measured": r["measured"],
                     "base": base[k][0], "base_fold": base[k][1],
                     "arm": arm[k][0], "arm_fold": arm[k][1],
                     "group": "strongly engaged" if engaged_all else "weak/none engaged"})
    print(f"{len(rows)} scored rows ({excluded} marker readouts excluded from both arms)")

    by_set = defaultdict(list)
    for r in rows:
        by_set[r["drug_set"]].append(r)
    sets = sorted(by_set)
    rng = random.Random(SEED)
    deltas = []
    for _ in range(N_BOOT):
        flat = [x for s in (rng.choice(sets) for _ in sets) for x in by_set[s]]
        deltas.append(ba(flat, "arm") - ba(flat, "base"))
    deltas.sort()
    lo, hi = deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1]
    d = ba(rows, "arm") - ba(rows, "base")
    acc = lambda col: sum(r[col] == r["measured"] for r in rows) / len(rows)  # noqa: E731
    print(f"full-knockout pins: acc {acc('base'):.4f}  balanced {ba(rows, 'base'):.4f}")
    print(f"engagement pins:    acc {acc('arm'):.4f}  balanced {ba(rows, 'arm'):.4f}")
    print(f"GATE: delta {d:+.4f} [{lo:+.4f}, {hi:+.4f}] over {len(sets)} drug sets -> "
          f"{'PASS' if lo > 0 else 'FAIL'}")
    nr = [r for r in rows if r["pathway"] != RAF]
    print(f"  (RAF set aside: delta {ba(nr, 'arm') - ba(nr, 'base'):+.4f})")

    print("\nreachable-case scorecard (reachable = full-knockout fold != 1):")
    for col, fc, name in (("base", "base_fold", "full-knockout"), ("arm", "arm_fold", "engagement")):
        se, sp, pr, auc, n = scorecard(rows, col, fc)
        print(f"  {name:>14}: sensitivity {se:.3f}  specificity {sp:.3f}  precision {pr:.3f}  "
              f"AUROC {auc:.3f}  (n={n})")

    moved = [r for r in rows if r["arm"] != r["base"]]
    print(f"\ncalls moved {len(moved)}: fixed {sum(r['arm'] == r['measured'] for r in moved)}, "
          f"broken {sum(r['base'] == r['measured'] for r in moved)}")
    for g in ("weak/none engaged", "strongly engaged"):
        mg = [r for r in moved if r["group"] == g]
        print(f"  {g:>18}: moved {len(mg)}, fixed {sum(r['arm'] == r['measured'] for r in mg)}, "
              f"broken {sum(r['base'] == r['measured'] for r in mg)}")


if __name__ == "__main__":
    main()
