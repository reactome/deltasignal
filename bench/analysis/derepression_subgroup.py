#!/usr/bin/env python3
"""Report the targeted subgroup: knockouts predicted to increase.

FR-009. An overall score can improve while the change does the wrong thing —
a bound tight enough to suppress every de-repression call would remove 50
wrong answers and 28 right ones, and look fine on the total. This reports
both halves so that cannot pass unnoticed.

Baseline (ceiling 11x): 78 such calls, 28 correct.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

NAME = {"0": "DOWN", "1": "NO_CHANGE", "2": "UP"}


def load(directory: Path) -> list[dict]:
    with (directory / "benchmark_cases.tsv").open(newline="") as handle:
        return [r for r in csv.DictReader(handle, delimiter="\t") if r["prediction"]]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--arms", type=Path, nargs="+", required=True)
    args = parser.parse_args()

    print(f"{'arm':16s} {'KO→UP':>7} {'right':>6} {'acc':>7} "
          f"{'wrong removed':>14} {'right lost':>11} {'overall':>9}")
    base_calls: set | None = None
    base_right: set | None = None
    for directory in args.arms:
        rows = load(directory)
        key = lambda r: (r["pathway_id"], r["gene"], r["key_output_dbid"], r["direction"])
        calls = {key(r) for r in rows if r["direction"] == "0" and r["prediction"] == "2"}
        right = {key(r) for r in rows
                 if r["direction"] == "0" and r["prediction"] == "2" and r["expected"] == "2"}
        overall = sum(1 for r in rows if r["prediction"] == r["expected"])
        if base_calls is None:
            base_calls, base_right = calls, right
            removed = lost = "-"
        else:
            # Of the baseline's calls, how many are gone — and were they the
            # wrong ones or the right ones?
            gone = base_calls - calls
            removed = str(len(gone - base_right))
            lost = str(len(gone & base_right))
        acc = len(right) / len(calls) if calls else float("nan")
        print(f"{directory.name:16s} {len(calls):7d} {len(right):6d} {acc:7.3f} "
              f"{removed:>14} {lost:>11} {overall:9d}")
    print("\n  'wrong removed' is the point; 'right lost' is the cost. An arm "
          "that\n  removes all 78 calls scores 22 on this subgroup and should "
          "not be\n  mistaken for an improvement.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
