#!/usr/bin/env python3
"""Scope of the specs/053 closed-twin join, measured DS-side on a built catalog.

A node X is CLOSED when it has an in-edge, every in-edge is an `output` edge,
and every producing reaction has an `input` from a node whose member leaves
strictly contain X's (X is released from a complex it was bound into). A FED
TWIN is a node with X's `nodes.csv` identity (every column but the uuid) that
is not closed and has an in-edge. specs/053 joins X to its fed twin when it has
exactly one (amendment 1).

For each valid case of a results directory, the case is TOUCHED when its pins
reach a fed twin but not the closed copy joined to it. Graph only.

Usage: closed_twin_census.py BUILD_DIR RESULTS_DIR
"""
import collections
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import twin_silo_census as t  # noqa: E402
from network_files import open_network  # noqa: E402


def joins(d):
    """(adjacency, {closed: [fed twins]}) for one pathway, every closed copy
    with at least one fed twin."""
    leaves, kind = {}, {}
    with open(d / "nodes.csv", newline="") as f:
        for r in csv.DictReader(f):
            leaves[r["uuid"]] = set(filter(None, (r["member_leaves"] or "").split("|"))) or {r["diagram_entity_id"]}
            kind[r["uuid"]] = r["node_kind"]
    key, twins, adj = t.load(d)
    inn = collections.defaultdict(set)
    for s, ts in adj.items():
        for x in ts:
            inn[x].add(s)
    etype = collections.defaultdict(set)
    with open_network(d) as f:
        for r in csv.DictReader(f):
            etype[(r["source_id"], r["target_id"])].add(r["edge_type"])

    def releases(p, x):
        return any("input" in etype[(s, p)] and s != x and leaves.get(x, set()) < leaves.get(s, set())
                   for s in inn[p])

    closed = set()
    for x in key:
        if kind.get(x) == "reaction" or not inn[x]:
            continue
        if all(etype[(p, x)] == {"output"} for p in inn[x]) and all(releases(p, x) for p in inn[x]):
            closed.add(x)
    out = {x: [w for w in twins.get(x, []) if w not in closed and inn[w]] for x in closed}
    return adj, {x: ws for x, ws in out.items() if ws}, len(closed)


def main(build, results):
    build, results = Path(build), Path(results)
    owner, info = {}, {}
    n_closed = 0
    per_path = collections.Counter()
    fan = collections.Counter()
    for d in sorted(p for p in build.iterdir() if (p / "nodes.csv").exists()):
        with open(d / "nodes.csv", newline="") as f:
            for r in csv.DictReader(f):
                owner[r["uuid"]] = d
        adj, j, nc = joins(d)
        n_closed += nc
        for ws in j.values():
            fan[len(ws)] += 1
        single = {x: ws[0] for x, ws in j.items() if len(ws) == 1}
        per_path[d.name] = len(single)
        info[d] = (adj, single)
    print(f"closed copies {n_closed}; with fed twins {sum(fan.values())}; "
          f"joined (exactly one fed twin) {fan[1]}")
    print("fed twins per closed copy:", sorted(fan.items()))
    print("joins by pathway:", ", ".join(f"{p} {n}" for p, n in per_path.most_common(12) if n))
    for axis in ("curator", "experimental"):
        tally, by_path = collections.Counter(), collections.Counter()
        with open(results / f"{axis}_cases.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["valid"] != "1":
                    continue
                pins = [u for u in r["gene_uuids"].split("|") if u in owner]
                if not pins:
                    continue
                adj, single = info[owner[pins[0]]]
                if not single:
                    continue
                seen = t.reach(adj, pins)
                if any(x not in seen and w in seen for x, w in single.items()):
                    tally["touched"] += 1
                    tally["touched, correct"] += r["predicted"] == r["expected"]
                    by_path[r["pathway"]] += 1
        print(f"== {axis}: touched {tally['touched']}, correct today {tally['touched, correct']}")
        print("   ", ", ".join(f"{p} {n}" for p, n in by_path.most_common(10)))


if __name__ == "__main__":
    main(*sys.argv[1:3])
