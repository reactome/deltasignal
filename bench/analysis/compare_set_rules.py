#!/usr/bin/env python3
"""Compare set-value combining rules on identical cases.

Every arm reports macro-F1, per-pathway net change, the both-arms-converged
count and the coverage delta — the last first, because a rule that scores
fewer cases can look better on accuracy while being worse.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

NAME = {"0": "DOWN", "1": "NO_CHANGE", "2": "UP"}


def load(directory: Path) -> list[dict]:
    with (directory / "benchmark_cases.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def key(row: dict) -> tuple:
    return (row["pathway_id"], row["gene"], row["key_output_dbid"], row["direction"])


def score(rows: list[dict], column: str = "prediction") -> tuple:
    rows = [r for r in rows if r[column]]
    correct = sum(1 for r in rows if r[column] == r["expected"])
    f1s = {}
    for label in "012":
        tp = sum(1 for r in rows if r[column] == label and r["expected"] == label)
        fp = sum(1 for r in rows if r[column] == label and r["expected"] != label)
        fn = sum(1 for r in rows if r[column] != label and r["expected"] == label)
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1s[NAME[label]] = (2 * precision * recall / (precision + recall)
                            if precision + recall else 0.0)
    macro = sum(f1s.values()) / len(f1s)
    return correct, len(rows), macro, f1s


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--arms", type=Path, nargs="+", required=True)
    args = parser.parse_args()

    base = load(args.baseline)
    base_by_key = {key(r): r for r in base}
    arms = [(p.name, load(p)) for p in args.arms]

    print(f"{'arm':18s} {'scored':>7} {'correct':>8} {'acc':>7} {'macro-F1':>9}"
          f" {'F1_DOWN':>8} {'F1_NC':>7} {'F1_UP':>7}")
    b = score(base)
    print(f"{args.baseline.name:18s} {b[1]:7d} {b[0]:8d} {b[0]/b[1]:7.4f} {b[2]:9.4f}"
          f" {b[3]['DOWN']:8.3f} {b[3]['NO_CHANGE']:7.3f} {b[3]['UP']:7.3f}")
    for name, rows in arms:
        s = score(rows)
        print(f"{name:18s} {s[1]:7d} {s[0]:8d} {s[0]/s[1]:7.4f} {s[2]:9.4f}"
              f" {s[3]['DOWN']:8.3f} {s[3]['NO_CHANGE']:7.3f} {s[3]['UP']:7.3f}")

    for name, rows in arms:
        print(f"\n## {name} vs {args.baseline.name}")
        shared = [r for r in rows
                  if r["prediction"] and base_by_key.get(key(r))
                  and base_by_key[key(r)]["prediction"]]
        # Coverage first: a rule scoring fewer cases can flatter itself.
        only_arm = sum(1 for r in rows if r["prediction"] and base_by_key.get(key(r))
                       and not base_by_key[key(r)]["prediction"])
        only_base = sum(1 for r in base if r["prediction"]
                        and not next((x for x in rows if key(x) == key(r)
                                      and x["prediction"]), None))
        print(f"  coverage: +{only_arm} scored only here, -{only_base} scored only in baseline")

        changed = [r for r in shared if r["prediction"] != base_by_key[key(r)]["prediction"]]
        both_conv = sum(1 for r in changed
                        if r["converged"] == "True"
                        and base_by_key[key(r)]["converged"] == "True")
        gained = sum(1 for r in changed if r["prediction"] == r["expected"])
        lost = sum(1 for r in changed
                   if base_by_key[key(r)]["prediction"] == r["expected"])
        print(f"  changed predictions: {len(changed)}  (both arms converged: {both_conv})")
        print(f"  net on changed: +{gained} -{lost} = {gained-lost:+d}")
        if changed:
            per = defaultdict(lambda: [0, 0])
            for r in changed:
                p = r["pathway_name"]
                if r["prediction"] == r["expected"]:
                    per[p][0] += 1
                elif base_by_key[key(r)]["prediction"] == r["expected"]:
                    per[p][1] += 1
            print(f"  {'pathway':34s} {'+':>4} {'-':>4} {'net':>5}")
            for p, (g, l) in sorted(per.items(), key=lambda kv: -(kv[1][0] - kv[1][1])):
                print(f"  {p[:34]:34s} {g:4d} {l:4d} {g-l:5d}")
        # Only set-resolved cases can differ; anything else is noise.
        modes = Counter(r["output_mapping_mode"] for r in changed)
        print(f"  changed by mapping mode: {dict(modes)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
