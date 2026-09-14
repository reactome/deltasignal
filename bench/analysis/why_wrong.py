#!/usr/bin/env python3
"""Attribute DeltaSignal's remaining errors to structural properties.

The control settled the big question: given MP-BioPath's own acyclic networks
the same propagator scores 538 of 740 against their 544, while on our
networks it scores 490. The propagator is not the problem, so something about
our network structure is. This asks what.

For every case it computes, on OUR network: whether the readout is
cycle-resident, the size of the largest strongly connected component on the
path, the shortest path length, and how many nodes the perturbation reaches.
Then it splits accuracy by each, and contrasts the cases we get wrong that
the control gets right against those we both get right.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cycle_structure import read_lng_network, tarjan_sccs  # noqa: E402


def network_facts(pathway_dir: Path):
    edges = read_lng_network(pathway_dir)
    adjacency: dict[str, list[str]] = defaultdict(list)
    nodes: set[str] = set()
    for e in edges:
        nodes.add(e["source"])
        nodes.add(e["target"])
        adjacency[e["source"]].append(e["target"])
    membership: dict[str, int] = {}
    sizes: dict[int, int] = defaultdict(int)
    for i, component in enumerate(tarjan_sccs(adjacency, nodes)):
        if len(component) > 1:
            for node in component:
                membership[node] = i
            sizes[i] = len(component)
    return adjacency, membership, dict(sizes)


def bfs(adjacency, sources, targets):
    """Shortest hop count and the largest SCC met on the way."""
    target_set = set(targets)
    depth = {s: 0 for s in sources}
    queue = deque(sources)
    reached: list[str] = []
    order: list[str] = list(sources)
    while queue:
        node = queue.popleft()
        if node in target_set:
            reached.append(node)
            continue
        for nxt in adjacency.get(node, ()):
            if nxt not in depth:
                depth[nxt] = depth[node] + 1
                order.append(nxt)
                queue.append(nxt)
    if not reached:
        return None, len(depth), order
    return min(depth[t] for t in reached), len(depth), order


def load_dbid_map(pathway_dir: Path) -> dict[str, list[str]]:
    """Same cascade as the benchmark: nodes.csv, else stid_to_uuid_mapping.csv.

    The adapted MP-BioPath catalog carries only the two files the server
    requires, so the fallback is what makes the control analysable.
    """
    mapping: dict[str, list[str]] = defaultdict(list)
    nodes_file = pathway_dir / "nodes.csv"
    if nodes_file.exists():
        with nodes_file.open(newline="") as handle:
            for row in csv.DictReader(handle):
                ids = {(row.get("diagram_entity_id") or "").strip()}
                ids |= {x.strip() for x in (row.get("member_leaves") or "").split("|")}
                for stable_id in ids - {""}:
                    mapping[stable_id.rsplit("-", 1)[-1]].append(row["uuid"])
    else:
        with (pathway_dir / "stid_to_uuid_mapping.csv").open(newline="") as handle:
            for row in csv.DictReader(handle):
                mapping[row["stable_id"].rsplit("-", 1)[-1]].append(row["uuid"])
    # A set-valued readout has no node of its own; without its member nodes
    # these cases read as "readout absent" and the no-path rate is inflated.
    resolution = pathway_dir / "node_resolution.csv"
    if resolution.exists():
        with resolution.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["relation"] == "set_member":
                    dbid = row["stable_id"].rsplit("-", 1)[-1]
                    if row["uuid"] not in mapping[dbid]:
                        mapping[dbid].append(row["uuid"])
    return mapping


def bucket_report(title, rows, bucket_fn):
    buckets: dict[object, list[dict]] = defaultdict(list)
    for r in rows:
        buckets[bucket_fn(r)].append(r)
    print(f"\n## {title}")
    print(f"  {'bucket':>18} {'n':>6} {'correct':>8} {'acc':>7}")
    for key in sorted(buckets, key=lambda k: (k is None, k)):
        group = buckets[key]
        correct = sum(1 for r in group if r["prediction"] == r["expected"])
        print(f"  {str(key):>18} {len(group):6d} {correct:8d} {correct/len(group):7.3f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True,
                        help="db_id_to_name_mapping.txt, to resolve gene names")
    parser.add_argument("--control", type=Path,
                        help="DS-on-MP-BioPath-networks results, for the contrast")
    args = parser.parse_args()

    with (args.results / "benchmark_cases.tsv").open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["prediction"]]

    dirs = {p.name.rsplit("_R-HSA-", 1)[-1]: p for p in args.catalog.iterdir() if p.is_dir()}

    from benchmark_mpbiopath_cases import load_gene_dbids, load_root_input_nodes
    gene_dbids = load_gene_dbids(args.id_map, {r["gene"] for r in rows})
    roots_by_pid = {pid: load_root_input_nodes(d) for pid, d in dirs.items()}
    cache: dict[str, tuple] = {}
    for row in rows:
        pid = row["pathway_id"].replace("R-HSA-", "")
        if pid not in cache:
            directory = dirs.get(pid)
            cache[pid] = ((*network_facts(directory), load_dbid_map(directory))
                          if directory else (None, None, None, None))
        adjacency, membership, sizes, dbid_map = cache[pid]
        if adjacency is None:
            row["_scc"] = row["_hops"] = row["_reach"] = None
            continue
        out_uuids = dbid_map.get(row["key_output_dbid"], [])
        # The benchmark records uuid COUNTS, not uuids, so recompute the
        # perturbed set the same way it does: gene -> dbids -> uuids,
        # restricted to root inputs.
        gene_uuids = sorted({
            u for dbid in gene_dbids.get(row["gene"], ())
            for u in dbid_map.get(dbid, ())
        } & roots_by_pid.get(pid, set()))
        row["_readout_cyclic"] = any(u in membership for u in out_uuids)
        row["_scc"] = max((sizes.get(membership.get(u), 0) for u in out_uuids), default=0)
        if gene_uuids:
            hops, reach, _ = bfs(adjacency, gene_uuids, out_uuids)
            row["_hops"], row["_reach"] = hops, reach
        else:
            row["_hops"] = row["_reach"] = None

    print(f"cases analysed: {len(rows)}")
    overall = sum(1 for r in rows if r["prediction"] == r["expected"])
    print(f"overall: {overall}/{len(rows)} = {overall/len(rows):.4f}")

    bucket_report("accuracy by whether the READOUT is cycle-resident", rows,
                  lambda r: r.get("_readout_cyclic"))

    def scc_bucket(r):
        n = r.get("_scc") or 0
        if n == 0:
            return "0 (acyclic)"
        if n < 50:
            return "1-49"
        if n < 200:
            return "50-199"
        if n < 600:
            return "200-599"
        return "600+"
    bucket_report("accuracy by the size of the readout's component", rows, scc_bucket)

    def reach_bucket(r):
        n = r.get("_reach")
        if n is None:
            return "unknown"
        for hi, label in ((50, "1  <50"), (150, "2  50-149"), (300, "3  150-299"),
                          (600, "4  300-599"), (1200, "5  600-1199")):
            if n < hi:
                return label
        return "6  1200+"
    bucket_report("accuracy by how many nodes the perturbation REACHES", rows, reach_bucket)

    def hop_bucket(r):
        n = r.get("_hops")
        if n is None:
            return "no path"
        for hi, label in ((5, "1  <5"), (9, "2  5-8"), (13, "3  9-12"), (19, "4  13-18")):
            if n < hi:
                return label
        return "5  19+"
    bucket_report("accuracy by shortest path length (hops)", rows, hop_bucket)

    bucket_report("accuracy by convergence", rows,
                  lambda r: "converged" if r["converged"] == "True" else "NOT converged")

    if args.control:
        with (args.control / "benchmark_cases.tsv").open(newline="") as handle:
            control = {(r["pathway_id"], r["gene"], r["key_output_dbid"], r["direction"]): r
                       for r in csv.DictReader(handle, delimiter="\t") if r["prediction"]}
        ours_wrong_theirs_right = []
        both_right = []
        for r in rows:
            k = (r["pathway_id"], r["gene"], r["key_output_dbid"], r["direction"])
            c = control.get(k)
            if not c:
                continue
            if r["prediction"] != r["expected"] and c["prediction"] == c["expected"]:
                ours_wrong_theirs_right.append(r)
            elif r["prediction"] == r["expected"] and c["prediction"] == c["expected"]:
                both_right.append(r)
        print(f"\n## the contrast: {len(ours_wrong_theirs_right)} cases our networks get "
              f"wrong that MP-BioPath's networks get right")
        print(f"   (against {len(both_right)} both get right)")
        for label, group in (("ours wrong / theirs right", ours_wrong_theirs_right),
                             ("both right", both_right)):
            cyc = sum(1 for r in group if r.get("_readout_cyclic"))
            nonconv = sum(1 for r in group if r["converged"] != "True")
            hops = [r["_hops"] for r in group if r.get("_hops") is not None]
            reach = [r["_reach"] for r in group if r.get("_reach") is not None]
            print(f"   {label:26s} n={len(group):4d}  readout cyclic {100*cyc/max(1,len(group)):5.1f}%"
                  f"  non-converged {100*nonconv/max(1,len(group)):5.1f}%"
                  f"  median hops {sorted(hops)[len(hops)//2] if hops else '-'}"
                  f"  median reach {sorted(reach)[len(reach)//2] if reach else '-'}")
        print("   by pathway:",
              dict(Counter(r["pathway_name"][:20] for r in ours_wrong_theirs_right).most_common(6)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
