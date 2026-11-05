#!/usr/bin/env python3
"""specs/138 amendment 1: korkut reported by target engagement, beside the raw score (Adam, 2026-11-05).

The raw score stays the headline. Beside it, the same metrics for the model, always-NOT and reachability on drug sets
in which every drug, some, or none engaged its target (engagement.py, its marker list frozen 2026-11-05). Each drug
set's own markers are left out of every split row (they are measured changed by definition when engaged); the `all`
row is the raw score. The split describes; it never chooses a rule.

Usage: engagement_report.py CASES > report.json
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import engagement  # noqa: E402
import score_134 as s134  # noqa: E402

SPLITS = ("every drug engaged", "some engaged", "none engaged")


def split_of(drug_set, engaged):
    states = {engaged[d] for d in drug_set.split("+")}
    return SPLITS[0] if states == {"yes"} else SPLITS[2] if "yes" not in states else SPLITS[1]


def markers_of(drug_set, markers):
    return {m for d in drug_set.split("+") for m, _ in markers[d]}


def methods(keys, cases):
    return {"model": s134.metrics([(cases[k]["predicted"], cases[k]["measured"]) for k in keys]),
            "always_not": s134.metrics([("not", cases[k]["measured"]) for k in keys]),
            "reachability": s134.metrics([(cases[k]["reach_predicted"], cases[k]["measured"]) for k in keys])}


def report(cases, engaged, markers):
    keys = sorted(k for k, r in cases.items() if r["status"] == "scored" and r["measured"] in s134.CLASSES)
    out = {"all (raw)": methods(keys, cases)}
    rest = [k for k in keys if k[2] not in markers_of(k[1], markers)]
    for s in SPLITS:
        sub = [k for k in rest if split_of(k[1], engaged) == s]
        out[s] = dict(methods(sub, cases), drug_sets=len({k[1] for k in sub}))
    return out


def main():
    engaged = {r["drug"]: r["engaged"] for r in engagement.table()}
    json.dump(report(s134.load(sys.argv[1]), engaged, engagement.MARKERS), sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
