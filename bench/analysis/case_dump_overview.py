#!/usr/bin/env python3
"""High-level overview of a per-case dump (no network access required).

Reports:
  - overall accuracy and confusion matrix (expected × predicted)
  - per-pathway accuracy and which class label dominates the errors
  - per-(pathway, key_output) error concentration — picks out readouts that
    fail on most/all perturbations (often indicates a structural network
    issue — missing edges or wrong sign at that specific readout)
  - continuous-output distribution per expected class — sanity check that
    UP-cases are predicted with higher UI than DOWN-cases on average

Usage
-----
  python case_dump_overview.py /tmp/ds_best_cases.tsv
"""
import csv, sys, statistics
from collections import Counter, defaultdict


CLS = {"0": "DOWN", "1": "NORM", "2": "UP"}


def main(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    n = len(rows)
    correct = sum(1 for r in rows if r["predicted"] == r["expected"])
    print(f"Total: {n}  Correct: {correct}  Acc: {correct*100/n:.2f}%")

    # Confusion matrix
    print("\n=== Confusion matrix (rows=expected, cols=predicted) ===")
    cm = defaultdict(lambda: defaultdict(int))
    for r in rows:
        cm[r["expected"]][r["predicted"]] += 1
    print(f"  {'':6s} {'DOWN':>6s} {'NORM':>6s} {'UP':>6s}")
    for e in ("0", "1", "2"):
        bits = " ".join(f"{cm[e][p]:>6d}" for p in ("0", "1", "2"))
        print(f"  {CLS[e]:6s} {bits}")

    # Per-pathway accuracy
    print("\n=== Per-pathway ===")
    pw = defaultdict(lambda: {"correct": 0, "total": 0,
                              "errs": Counter()})
    for r in rows:
        pw[r["pathway"]]["total"] += 1
        if r["predicted"] == r["expected"]:
            pw[r["pathway"]]["correct"] += 1
        else:
            pw[r["pathway"]]["errs"][f"E{CLS[r['expected']]}_P{CLS[r['predicted']]}"] += 1
    for p in sorted(pw, key=lambda p: -pw[p]["total"]):
        d = pw[p]
        acc = d["correct"] * 100 / d["total"]
        worst = d["errs"].most_common(2)
        worst_s = ", ".join(f"{k}={v}" for k, v in worst)
        print(f"  {p[:38]:40s} {d['correct']:3d}/{d['total']:3d} "
              f"({acc:5.1f}%)  top errs: {worst_s}")

    # Readout concentration: which (pathway, key_output) are mostly wrong?
    print("\n=== Readouts with ≥4 wrong cases (structural candidates) ===")
    by_ko = defaultdict(lambda: {"wrong": 0, "total": 0})
    for r in rows:
        k = (r["pathway"], r["key_output"])
        by_ko[k]["total"] += 1
        if r["predicted"] != r["expected"]:
            by_ko[k]["wrong"] += 1
    bad = [(k, v) for k, v in by_ko.items() if v["wrong"] >= 4]
    bad.sort(key=lambda x: -x[1]["wrong"])
    for (pw_, ko), v in bad[:20]:
        print(f"  {pw_[:30]:32s} → {ko:8s}  "
              f"wrong={v['wrong']:2d}/{v['total']:2d}")

    # Continuous-output distribution per expected class
    print("\n=== pred_ui distribution by expected class ===")
    by_e = defaultdict(list)
    for r in rows:
        by_e[r["expected"]].append(float(r["pred_ui"]))
    for e in ("0", "1", "2"):
        vs = by_e[e]
        if not vs: continue
        vs.sort()
        p25 = vs[len(vs)//4]; p50 = vs[len(vs)//2]; p75 = vs[3*len(vs)//4]
        print(f"  exp={CLS[e]:4s} (n={len(vs):3d})  min={vs[0]:6.2f} "
              f"p25={p25:6.2f} p50={p50:6.2f} p75={p75:6.2f} max={vs[-1]:6.2f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr); sys.exit(1)
    main(sys.argv[1])
