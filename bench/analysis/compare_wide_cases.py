#!/usr/bin/env python3
"""Compare two `benchmark_vs_mpbiopath.py --dump-cases` arms, conditioned.

The wide curator set is the only evidence that has held up at scale, but a
naive diff of two arms answers the wrong question whenever a change alters the
network: the benchmark perturbs root inputs, so any edge-removing arm promotes
new nodes into the perturbed set and the two arms are then scoring different
experiments. `n_gene_uuids` and `n_ko_uuids` are that experiment's fingerprint,
so the default report is restricted to cases where both are unchanged.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def key(row: dict) -> tuple:
    return (row["pathway"], row["gene"], row["direction"], row["key_output"])


def load(path: Path) -> dict:
    with path.open() as fh:
        return {key(r): r for r in csv.DictReader(fh, delimiter="\t")}


def macro_f1(rows, labels=None) -> float:
    # Derive the classes from the data rather than hard-coding them: this dump
    # encodes DOWN/NORMAL/UP as 0/1/2 while the ten-pathway one uses -1/0/1,
    # and a hard-coded set silently drops a whole class from the average.
    #
    # `labels` must be passed when comparing two arms. Derived per-arm, an arm
    # that emits one stray prediction of a class the other never emits divides
    # its F1 sum by a larger denominator and is penalised for it — in a tool
    # whose only job is a paired comparison, that is the wrong denominator.
    if labels is None:
        labels = {r["expected"] for r in rows} | {r["predicted"] for r in rows}
    if not labels:
        return 0.0
    total = 0.0
    for label in labels:
        tp = sum(1 for r in rows if r["predicted"] == label and r["expected"] == label)
        fp = sum(1 for r in rows if r["predicted"] == label and r["expected"] != label)
        fn = sum(1 for r in rows if r["predicted"] != label and r["expected"] == label)
        denom = 2 * tp + fp + fn
        total += (2 * tp / denom) if denom else 0.0
    return total / len(labels)


def scored(rows):
    return [r for r in rows if r.get("valid") == "1"]


def summarise(name: str, rows, labels=None) -> None:
    ok = sum(1 for r in rows if r["predicted"] == r["expected"])
    print(f"{name:<14} scored {len(rows):>6}  correct {ok:>6}  "
          f"acc {ok / len(rows) if rows else 0:.4f}  "
          f"macro-F1 {macro_f1(rows, labels):.4f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", type=Path, required=True)
    ap.add_argument("--arm", type=Path, required=True)
    ap.add_argument("--all-cases", action="store_true",
                    help="report the unconditioned diff as well; the experiment "
                         "may differ between arms, so read it as an upper bound")
    args = ap.parse_args()

    base, arm = load(args.baseline), load(args.arm)
    shared = [k for k in base if k in arm]

    same = [k for k in shared
            if base[k]["n_gene_uuids"] == arm[k]["n_gene_uuids"]
            and base[k]["n_ko_uuids"] == arm[k]["n_ko_uuids"]]
    moved = len(shared) - len(same)

    base_rows = scored([base[k] for k in same])
    arm_rows = scored([arm[k] for k in same])
    valid_keys = {key(r) for r in base_rows} & {key(r) for r in arm_rows}
    base_rows = [r for r in base_rows if key(r) in valid_keys]
    arm_rows = [r for r in arm_rows if key(r) in valid_keys]

    print(f"cases: baseline {len(base)}  arm {len(arm)}  shared {len(shared)}")
    print(f"SAME EXPERIMENT: {len(same)} of {len(shared)}"
          + (f"   (experiment MOVED in {moved})" if moved else ""))
    print(f"scored in both: {len(valid_keys)}\n")
    labels = ({r["expected"] for r in base_rows} | {r["predicted"] for r in base_rows}
              | {r["expected"] for r in arm_rows} | {r["predicted"] for r in arm_rows})
    summarise("baseline", base_rows, labels)
    summarise("arm", arm_rows, labels)

    by_key = {key(r): r for r in arm_rows}
    gained = lost = 0
    per_pathway: Counter = Counter()
    for r in base_rows:
        a = by_key[key(r)]
        if r["predicted"] == a["predicted"]:
            continue
        was = r["predicted"] == r["expected"]
        now = a["predicted"] == a["expected"]
        if now and not was:
            gained += 1
            per_pathway[r["pathway"]] += 1
        elif was and not now:
            lost += 1
            per_pathway[r["pathway"]] -= 1
    changed = sum(1 for r in base_rows if r["predicted"] != by_key[key(r)]["predicted"])
    print(f"\nchanged predictions: {changed}   net {gained - lost:+d} "
          f"(+{gained} / -{lost})")
    if per_pathway:
        print("per-pathway net (nonzero):")
        for pathway, net in sorted(per_pathway.items(), key=lambda kv: kv[1]):
            if net:
                print(f"  {net:+4d}  {pathway}")

    if args.all_cases:
        print("\n--- unconditioned (experiment may differ) ---")
        b_all, a_all = scored([base[k] for k in shared]), scored([arm[k] for k in shared])
        all_labels = ({r["expected"] for r in b_all} | {r["predicted"] for r in b_all}
                      | {r["expected"] for r in a_all} | {r["predicted"] for r in a_all})
        summarise("baseline", b_all, all_labels)
        summarise("arm", a_all, all_labels)


if __name__ == "__main__":
    main()
