#!/usr/bin/env python3
"""specs/178: how the Lun predictions differ between two releases' builds. Values-blind (prediction files only).

Usage: compare_178.py OLD_PREDICTIONS NEW_PREDICTIONS OUT_JSON
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def pairs(path):
    """{(gene, readout): (model call, reach call)} over the networks that score the pair, as the scorers combine them."""
    import lun_score as L
    frozen, reach = defaultdict(list), defaultdict(list)
    n_rows = n_scored = 0
    for r in csv.DictReader(open(path), delimiter="\t"):
        n_rows += 1
        if r["self"] == "1" or r["status"] != "scored":
            continue
        n_scored += 1
        frozen[(r["gene"], r["readout"])].append(r["predicted"])
        reach[(r["gene"], r["readout"])].append(r["reach_predicted"])
    return ({k: (L.combine_networks(frozen[k]), L.combine_networks(reach[k])) for k in frozen}, n_rows, n_scored)


def compare(old, new):
    both = old.keys() & new.keys()
    moved = [k for k in both if old[k][0] != new[k][0]]
    reach_moved = [k for k in both if old[k][1] != new[k][1]]
    return {
        "pairs_old": len(old), "pairs_new": len(new), "pairs_both": len(both),
        "only_old": len(old.keys() - new.keys()), "only_new": len(new.keys() - old.keys()),
        "proteins_old": len({g for g, _ in old}), "proteins_new": len({g for g, _ in new}),
        "model_calls_moved": len(moved),
        "model_transitions": dict(Counter(f"{old[k][0]}->{new[k][0]}" for k in moved).most_common()),
        "model_moved_by_readout": dict(Counter(k[1] for k in moved).most_common()),
        "reach_calls_moved": len(reach_moved),
        "reach_transitions": dict(Counter(f"{old[k][1]}->{new[k][1]}" for k in reach_moved).most_common()),
        "only_new_by_readout": dict(Counter(k[1] for k in new.keys() - old.keys()).most_common()),
        "only_old_by_readout": dict(Counter(k[1] for k in old.keys() - new.keys()).most_common()),
    }


def main():
    old_f, new_f, out = sys.argv[1:4]
    (old, ro, so), (new, rn, sn) = pairs(old_f), pairs(new_f)
    res = {"rows_old": ro, "rows_new": rn, "scored_rows_old": so, "scored_rows_new": sn, **compare(old, new)}
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
