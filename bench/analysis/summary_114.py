#!/usr/bin/env python3
"""specs/114 part 2: the registered korkut figures, and the multi-network count for both benches.

For a case log (korkut_cases.tsv or decryptm_cases.tsv) over scored cases whose measured call is not `conflict`:
- accuracy of the model (`predicted`), of reachability (`reach_predicted`) and of always-NOT;
- regulated cases (measured down or up) right, for the model and reachability;
- exact McNemar of model against reachability;
- distinct (perturbation, readout) units, and how many get different model calls from different networks
  (a reaction under two top-level pathways is in both networks).

Usage: summary_114.py korkut|decryptm CASES_TSV
"""
import collections
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from decryptm_bench import mcnemar_p  # noqa: E402

UNIT = {"korkut": ("drug_set", "readout"), "decryptm": ("drug", "uniprot", "position")}


def summarize(kind, rows):
    rows = [r for r in rows if r["status"] == "scored" and r["measured"] != "conflict"]
    reg = [r for r in rows if r["measured"] in ("down", "up")]
    model_only = sum(r["predicted"] == r["measured"] != r["reach_predicted"] for r in rows)
    reach_only = sum(r["reach_predicted"] == r["measured"] != r["predicted"] for r in rows)
    calls = collections.defaultdict(set)
    for r in rows:
        calls[tuple(r[k] for k in UNIT[kind])].add(r["predicted"])
    return {
        "scored": len(rows),
        "accuracy": {"model": sum(r["predicted"] == r["measured"] for r in rows) / len(rows),
                     "reach": sum(r["reach_predicted"] == r["measured"] for r in rows) / len(rows),
                     "always_not": sum(r["measured"] == "not" for r in rows) / len(rows)},
        "regulated": {"n": len(reg), "model": sum(r["predicted"] == r["measured"] for r in reg),
                      "reach": sum(r["reach_predicted"] == r["measured"] for r in reg)},
        "mcnemar": {"model_only": model_only, "reach_only": reach_only, "p": mcnemar_p(model_only, reach_only)},
        "units": len(calls),
        "units_with_disagreeing_networks": sum(len(c) > 1 for c in calls.values()),
    }


def main():
    kind, path = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    print(json.dumps(summarize(kind, rows), indent=1))


if __name__ == "__main__":
    main()
