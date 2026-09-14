#!/usr/bin/env python3
"""Classify each benchmark case by whether it touches a cycle.

Feature 004 must report DeltaSignal against MP-BioPath separately on cyclic
and acyclic cases, because our networks are 27.1% cycle-resident and
MP-BioPath's published ones are 0.8%. Pathway-level cyclicity is not a
substitute: a pathway can be 36% cyclic while a given case never touches a
cycle.

Runs against a completed benchmark output directory, so arms already scored
do not have to be re-run.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict, deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cycle_structure import read_lng_network, tarjan_sccs  # noqa: E402


def scc_membership(pathway_dir: Path) -> dict[str, int]:
    edges = read_lng_network(pathway_dir)
    adjacency: dict[str, list[str]] = defaultdict(list)
    nodes: set[str] = set()
    for e in edges:
        nodes.add(e["source"])
        nodes.add(e["target"])
        adjacency[e["source"]].append(e["target"])
    membership: dict[str, int] = {}
    for i, component in enumerate(tarjan_sccs(adjacency, nodes)):
        if len(component) > 1:
            for node in component:
                membership[node] = i
    return membership


def component_sizes(membership: dict[str, int]) -> dict[int, int]:
    sizes: dict[int, int] = defaultdict(int)
    for cid in membership.values():
        sizes[cid] += 1
    return sizes


def largest_scc_on_path(adjacency, sources, targets, membership, sizes) -> int:
    """Largest component touched by any shortest path from source to target.

    Plain BFS over nodes, unlike the signed traversal the baselines use: we
    want whether a cycle is on the route, not what the route's polarity is.
    """
    if not sources or not targets:
        return 0
    target_set = set(targets)
    parents: dict[str, list[str]] = defaultdict(list)
    depth = {s: 0 for s in sources}
    queue = deque(sources)
    reached_at = None
    while queue:
        node = queue.popleft()
        if reached_at is not None and depth[node] >= reached_at:
            continue
        if node in target_set:
            reached_at = depth[node]
            continue
        for nxt in adjacency.get(node, ()):
            if nxt not in depth:
                depth[nxt] = depth[node] + 1
                parents[nxt].append(node)
                queue.append(nxt)
            elif depth[nxt] == depth[node] + 1:
                parents[nxt].append(node)
    if reached_at is None:
        return 0
    # Walk the shortest-path DAG backwards from every reached target.
    seen: set[str] = set()
    stack = [t for t in target_set if depth.get(t) == reached_at]
    largest = 0
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        cid = membership.get(node)
        if cid is not None:
            largest = max(largest, sizes[cid])
        stack.extend(parents.get(node, ()))
    return largest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True,
                        help="A benchmark output directory holding benchmark_cases.tsv")
    parser.add_argument("--catalog", type=Path, required=True,
                        help="The catalog that arm was solved on")
    parser.add_argument("--out", type=Path,
                        help="Write the annotated table here (default: alongside the input)")
    args = parser.parse_args()

    cases_path = args.results / "benchmark_cases.tsv"
    with cases_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    dirs = {p.name.rsplit("_", 1)[-1]: p for p in args.catalog.iterdir() if p.is_dir()}
    cache: dict[str, tuple] = {}
    for row in rows:
        pid = row["pathway_id"]
        if pid not in cache:
            directory = dirs.get(pid)
            if directory is None:
                cache[pid] = ({}, {}, {})
            else:
                membership = scc_membership(directory)
                sizes = component_sizes(membership)
                adjacency: dict[str, list[str]] = defaultdict(list)
                for e in read_lng_network(directory):
                    adjacency[e["source"]].append(e["target"])
                cache[pid] = (membership, sizes, adjacency)
        membership, sizes, adjacency = cache[pid]
        # The benchmark records counts, not the uuids themselves, so cyclicity
        # is derived from the path columns it does record plus the network.
        gene_uuids = [u for u in (row.get("gene_uuids") or "").split("|") if u]
        out_uuids = [u for u in (row.get("output_uuids") or "").split("|") if u]
        row["readout_in_scc"] = str(any(u in membership for u in out_uuids))
        row["perturbation_in_scc"] = str(any(u in membership for u in gene_uuids))
        row["largest_scc_on_path"] = str(
            largest_scc_on_path(adjacency, gene_uuids, out_uuids, membership, sizes))
        row["case_is_cyclic"] = str(
            row["readout_in_scc"] == "True" or row["largest_scc_on_path"] != "0")

    out = args.out or (args.results / "benchmark_cases_cyclicity.tsv")
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    scored = [r for r in rows if r["prediction"]]
    cyclic = [r for r in scored if r["case_is_cyclic"] == "True"]
    print(f"wrote {out}")
    print(f"scored cases {len(scored)}: cyclic {len(cyclic)}, acyclic {len(scored)-len(cyclic)}")
    if not any(r.get("gene_uuids") for r in rows):
        print("WARNING: no gene_uuids column in this results file — the "
              "classification above is vacuous. Re-run the benchmark with the "
              "uuid columns emitted before trusting any cyclic/acyclic split.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
