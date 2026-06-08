#!/usr/bin/env python3
"""Diff two per-case dumps to see what flipped between configurations.

Each dump is a TSV with one row per case:
    pathway gene direction key_output predicted expected method pred_ui

This tool joins on (pathway, gene, direction, key_output) and reports:
  - aggregate accuracy delta
  - per-pathway accuracy delta
  - cases that flipped: WRONG→RIGHT (gain), RIGHT→WRONG (loss), and
    WRONG→WRONG with a changed prediction (sideways)

Usage
-----
  python compare_configs.py baseline.tsv candidate.tsv [--label-a NAME --label-b NAME]
"""
import csv, sys
from collections import defaultdict


CLS = {"0": "DOWN", "1": "NORM", "2": "UP"}


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return {(r["pathway"], r["gene"], r["direction"], r["key_output"]): r
            for r in rows}


def main(a_path, b_path, label_a="A", label_b="B"):
    A = load(a_path)
    B = load(b_path)
    keys = set(A) & set(B)
    only_a = set(A) - set(B)
    only_b = set(B) - set(A)
    if only_a or only_b:
        print(f"  [warn] {len(only_a)} cases only in A, {len(only_b)} only in B",
              file=sys.stderr)

    a_correct = sum(1 for k in A if A[k]["predicted"] == A[k]["expected"])
    b_correct = sum(1 for k in B if B[k]["predicted"] == B[k]["expected"])
    print(f"{label_a}: {a_correct}/{len(A)} = {a_correct*100/len(A):.2f}%")
    print(f"{label_b}: {b_correct}/{len(B)} = {b_correct*100/len(B):.2f}%")
    print(f"Δ      : {b_correct - a_correct:+d}")

    # Per-pathway delta
    print("\n=== Per-pathway accuracy delta ===")
    pw_acc = defaultdict(lambda: {"a": [0, 0], "b": [0, 0]})
    for k in keys:
        pw = A[k]["pathway"]
        pw_acc[pw]["a"][0] += int(A[k]["predicted"] == A[k]["expected"])
        pw_acc[pw]["a"][1] += 1
        pw_acc[pw]["b"][0] += int(B[k]["predicted"] == B[k]["expected"])
        pw_acc[pw]["b"][1] += 1
    for pw in sorted(pw_acc):
        a, b = pw_acc[pw]["a"], pw_acc[pw]["b"]
        ap, bp = a[0]*100/a[1] if a[1] else 0, b[0]*100/b[1] if b[1] else 0
        print(f"  {pw[:38]:40s} {a[0]:3d}/{a[1]:3d} ({ap:5.1f}%) → "
              f"{b[0]:3d}/{b[1]:3d} ({bp:5.1f}%)  Δ={b[0]-a[0]:+d}")

    # Flip detail
    gain, loss, sideways = [], [], []
    for k in keys:
        ra, rb = A[k], B[k]
        a_ok = ra["predicted"] == ra["expected"]
        b_ok = rb["predicted"] == rb["expected"]
        if not a_ok and b_ok:    gain.append((ra, rb))
        elif a_ok and not b_ok:  loss.append((ra, rb))
        elif not a_ok and not b_ok and ra["predicted"] != rb["predicted"]:
            sideways.append((ra, rb))

    print(f"\n=== Flips ===")
    print(f"  Gains  (wrong→right): {len(gain)}")
    print(f"  Losses (right→wrong): {len(loss)}")
    print(f"  Sideways  (wrong→wrong, different): {len(sideways)}")

    def fmt(r):
        d = "KO" if r["direction"] == "0" else "OE"
        return (f"{r['pathway'][:30]:32s} {r['gene']:9s} {d}  →{r['key_output']:8s}"
                f"  exp={CLS[r['expected']]:4s} pred={CLS[r['predicted']]:4s} "
                f"ui={float(r['pred_ui']):.2f}")

    if gain:
        print(f"\n--- Gains ({label_a} → {label_b}) ---")
        for ra, rb in gain[:50]:
            print(f"  {label_a}: {fmt(ra)}\n  {label_b}: {fmt(rb)}\n")
    if loss:
        print(f"\n--- Losses ({label_a} → {label_b}) ---")
        for ra, rb in loss[:50]:
            print(f"  {label_a}: {fmt(ra)}\n  {label_b}: {fmt(rb)}\n")


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__, file=sys.stderr); sys.exit(1)
    a, b = args[0], args[1]
    la, lb = "A", "B"
    if "--label-a" in args:
        la = args[args.index("--label-a") + 1]
    if "--label-b" in args:
        lb = args[args.index("--label-b") + 1]
    main(a, b, la, lb)
