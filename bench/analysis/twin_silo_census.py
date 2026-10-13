#!/usr/bin/env python3
"""Exact identity twins on opposite sides of a perturbation (specs/053 trace).

Two nodes of one pathway are identity twins when their `nodes.csv` rows are
equal in every column but the uuid (node_kind, diagram_entity_id,
compartment, member_leaves, source_sets, chosen_members). The NF1 trace
found SPRED dimer:NF1 split this way: the copy made from NF1 goes only to
ubiquitination, while the copy that acts as a RAS GAP is fed only by its own
release, so an NF1 knockout never reaches the GAP.

For every valid case of a results directory, this reads the pinned uuids the
bench logged (`gene_uuids`), takes everything the pins reach (any edge), and
counts the case as SPLIT when some reached node has a twin the pins do not
reach. Graph only; nothing is solved, so this bounds what a join could move.

Usage: twin_silo_census.py BUILD_DIR RESULTS_DIR
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from network_files import open_network  # noqa: E402

KEY = ("node_kind", "diagram_entity_id", "compartment", "member_leaves", "source_sets", "chosen_members")


def load(d):
    key = {}
    with open(d / "nodes.csv", newline="") as f:
        for r in csv.DictReader(f):
            key[r["uuid"]] = tuple(r.get(k, "") for k in KEY)
    groups = collections.defaultdict(list)
    for u, k in key.items():
        groups[k].append(u)
    twins = {u: [v for v in groups[key[u]] if v != u] for u in key if len(groups[key[u]]) > 1}
    adj = collections.defaultdict(set)
    with open_network(d) as f:
        for r in csv.DictReader(f):
            adj[r["source_id"]].add(r["target_id"])
    return key, twins, adj


def reach(adj, src):
    seen, stack = set(src), list(src)
    while stack:
        for w in adj.get(stack.pop(), ()):
            if w not in seen:
                seen.add(w); stack.append(w)
    return seen


def main(build, results):
    build, results = Path(build), Path(results)
    owner, cache = {}, {}
    n_groups = n_twin_nodes = 0
    for d in sorted(p for p in build.iterdir() if (p / "nodes.csv").exists()):
        with open(d / "nodes.csv", newline="") as f:
            for r in csv.DictReader(f):
                owner[r["uuid"]] = d
    for axis in ("curator", "experimental"):
        tally, by_path, genes = collections.Counter(), collections.Counter(), collections.Counter()
        with open(results / f"{axis}_cases.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["valid"] != "1":
                    continue
                tally["valid"] += 1
                pins = [u for u in r["gene_uuids"].split("|") if u in owner]
                if not pins:
                    continue
                d = owner[pins[0]]
                if d not in cache:
                    cache[d] = load(d)
                key, twins, adj = cache[d]
                seen = reach(adj, pins)
                split = any(t not in seen for u in seen for t in twins.get(u, ()))
                if split:
                    ok = r["predicted"] == r["expected"]
                    tally["split"] += 1
                    tally["split, correct"] += ok
                    tally["split, no change predicted"] += r["predicted"] == "1"
                    tally["split, no change predicted, wrong"] += (r["predicted"] == "1") and not ok
                    by_path[r["pathway"]] += 1
                    genes[(r["pathway"], r["gene"])] += 1
        print(f"== {axis}")
        for k in ("valid", "split", "split, correct", "split, no change predicted",
                  "split, no change predicted, wrong"):
            print(f"  {k:40s} {tally[k]}")
        print("  pathways:", ", ".join(f"{p} {n}" for p, n in by_path.most_common(8)))
        print(f"  distinct (pathway, gene): {len(genes)}")
    for d, (key, twins, adj) in cache.items():
        n_twin_nodes += len(twins)
    print(f"twin nodes in the {len(cache)} pathways scored: {n_twin_nodes}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
