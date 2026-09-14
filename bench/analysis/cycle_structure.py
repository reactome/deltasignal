#!/usr/bin/env python3
"""Report the cycle structure of a network catalog.

The comparison against MP-BioPath was never like-for-like: its published
networks are effectively acyclic (zero self-loops, largest strongly connected
component eleven nodes across the ten benchmark pathways) while ours hold
components of 836, 643, 465 and 443 nodes. This tool is the reproducible
measurement behind that claim, for both formats.

It also classifies components. A genuine feedback loop is a chain of distinct
curated reactions, so its node count and its reaction count are of the same
order. A recycling artifact is one or two curated reactions instantiated as
hundreds of virtual reactions, so nodes-per-reaction is large. Polarity is
reported but is deliberately NOT the rule: 21 of 34 components in the current
catalog contain no negative edge at all, yet ERBB2's 465-node artifact
contains 121, so polarity alone misclassifies.

See specs/004-loop-handling/contracts/cycle-structure-report.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

DEFAULT_RATIO_THRESHOLD = 15.0


def tarjan_sccs(adjacency: dict[str, list[str]], nodes: set[str]) -> list[list[str]]:
    """Iterative Tarjan. Recursive descent overflows on an 836-node component."""
    index: dict[str, int] = {}
    low: dict[str, int] = {}
    on_stack: dict[str, bool] = {}
    stack: list[str] = []
    counter = 0
    components: list[list[str]] = []

    for root in sorted(nodes):
        if root in index:
            continue
        work: list[tuple[str, int]] = [(root, 0)]
        while work:
            node, child_i = work[-1]
            if child_i == 0:
                index[node] = low[node] = counter
                counter += 1
                stack.append(node)
                on_stack[node] = True
            recursed = False
            children = adjacency.get(node, ())
            for i in range(child_i, len(children)):
                child = children[i]
                if child not in index:
                    work[-1] = (node, i + 1)
                    work.append((child, 0))
                    recursed = True
                    break
                if on_stack.get(child):
                    low[node] = min(low[node], index[child])
            if recursed:
                continue
            if low[node] == index[node]:
                component = []
                while True:
                    member = stack.pop()
                    on_stack[member] = False
                    component.append(member)
                    if member == node:
                        break
                components.append(sorted(component))
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[node])
    return components


def read_lng_network(path: Path) -> list[dict[str, str]]:
    with (path / "logic_network.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"source_id", "target_id"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path.name}/logic_network.csv missing columns: {missing}")
        return [
            {
                "source": (row["source_id"] or "").strip(),
                "target": (row["target_id"] or "").strip(),
                "sign": (row.get("pos_neg") or "").strip(),
                "edge_type": (row.get("edge_type") or "").strip(),
                "reaction": (row.get("edge_reaction_id") or "").strip(),
            }
            for row in reader
        ]


def read_mpbiopath_network(path: Path) -> list[dict[str, str]]:
    """Four columns: parent dbid, child dbid, polarity (1/-1), conjunction (0=AND).

    No header. Reaction provenance is not carried in this format, so the
    nodes-per-reaction ratio is undefined for these networks and the
    classification is skipped rather than faked.
    """
    edges = []
    with path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"{path.name}:{lineno} has fewer than two columns")
            polarity = parts[2].strip() if len(parts) > 2 else "1"
            edges.append(
                {
                    "source": parts[0].strip(),
                    "target": parts[1].strip(),
                    "sign": "neg" if polarity == "-1" else "pos",
                    "edge_type": "",
                    "reaction": "",
                }
            )
    return edges


def analyse(name: str, edges: list[dict[str, str]], threshold: float) -> dict:
    adjacency: dict[str, list[str]] = defaultdict(list)
    nodes: set[str] = set()
    self_loops = 0
    for edge in edges:
        source, target = edge["source"], edge["target"]
        nodes.add(source)
        nodes.add(target)
        if source == target:
            self_loops += 1
        adjacency[source].append(target)

    components = [c for c in tarjan_sccs(adjacency, nodes) if len(c) > 1]
    membership: dict[str, int] = {}
    for i, component in enumerate(components):
        for node in component:
            membership[node] = i

    detail = []
    for i, component in enumerate(components):
        intra = [
            e for e in edges
            if membership.get(e["source"]) == i and membership.get(e["target"]) == i
        ]
        reactions = {e["reaction"] for e in intra}
        # An empty reaction id is a real distinct value here: it marks bridges
        # with no owning reaction, and collapsing it into "unknown" would
        # understate the reaction count and inflate the ratio.
        ratio = len(component) / len(reactions) if reactions else None
        detail.append(
            {
                "nodes": len(component),
                "intra_edges": len(intra),
                "negative_edges": sum(1 for e in intra if e["sign"] == "neg"),
                "reaction_stids": len(reactions),
                "nodes_per_reaction": ratio,
                "edge_types": dict(Counter(e["edge_type"] for e in intra)),
                "classification": (
                    None if ratio is None
                    else "recycling_artifact" if ratio >= threshold
                    else "candidate_feedback"
                ),
            }
        )
    detail.sort(key=lambda d: -d["nodes"])

    return {
        "pathway": name,
        "nodes": len(nodes),
        "edges": len(edges),
        "self_loops": self_loops,
        "components": len(components),
        "cycle_resident_nodes": sum(len(c) for c in components),
        "largest_component": max((len(c) for c in components), default=0),
        "component_detail": detail,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog", type=Path,
                        help="Directory of LNG pathway directories")
    parser.add_argument("--mpbiopath", type=Path,
                        help="Directory of MP-BioPath four-column TSVs")
    parser.add_argument("--classify", action="store_true",
                        help="Report per-component classification")
    parser.add_argument("--ratio-threshold", type=float, default=DEFAULT_RATIO_THRESHOLD,
                        help=f"nodes-per-reaction cut (default {DEFAULT_RATIO_THRESHOLD})")
    parser.add_argument("--json", type=Path, help="Also write machine-readable output")
    args = parser.parse_args()

    # Exactly one source. Defaulting here would silently measure the wrong
    # network set, which is the specific mistake this feature exists to undo.
    if bool(args.catalog) == bool(args.mpbiopath):
        parser.error("give exactly one of --catalog or --mpbiopath")

    results = []
    if args.catalog:
        directories = sorted(p for p in args.catalog.iterdir()
                             if p.is_dir() and (p / "logic_network.csv").exists())
        if not directories:
            print(f"No pathway directories under {args.catalog}", file=sys.stderr)
            return 1
        for directory in directories:
            results.append(analyse(directory.name, read_lng_network(directory),
                                   args.ratio_threshold))
        source = f"LNG catalog {args.catalog}"
    else:
        files = sorted(args.mpbiopath.glob("*.tsv"))
        if not files:
            print(f"No .tsv files under {args.mpbiopath}", file=sys.stderr)
            return 1
        for path in files:
            results.append(analyse(path.stem, read_mpbiopath_network(path),
                                   args.ratio_threshold))
        source = f"MP-BioPath networks {args.mpbiopath}"

    print(f"# {source}")
    print(f"# nodes-per-reaction threshold in force: {args.ratio_threshold}")
    header = (f"{'pathway':52s} {'nodes':>7} {'edges':>7} {'self':>5} "
              f"{'comps':>6} {'inSCC':>7} {'%':>6} {'largest':>8}")
    print(header)
    print("-" * len(header))
    for r in results:
        pct = 100 * r["cycle_resident_nodes"] / r["nodes"] if r["nodes"] else 0.0
        print(f"{r['pathway'][:52]:52s} {r['nodes']:7d} {r['edges']:7d} "
              f"{r['self_loops']:5d} {r['components']:6d} "
              f"{r['cycle_resident_nodes']:7d} {pct:5.1f}% {r['largest_component']:8d}")
    print("-" * len(header))
    total_nodes = sum(r["nodes"] for r in results)
    total_scc = sum(r["cycle_resident_nodes"] for r in results)
    print(f"{'TOTAL':52s} {total_nodes:7d} {sum(r['edges'] for r in results):7d} "
          f"{sum(r['self_loops'] for r in results):5d} "
          f"{sum(r['components'] for r in results):6d} {total_scc:7d} "
          f"{100 * total_scc / total_nodes if total_nodes else 0:5.1f}% "
          f"{max((r['largest_component'] for r in results), default=0):8d}")

    if args.classify:
        print(f"\n# component classification at threshold {args.ratio_threshold}")
        cls_header = (f"{'pathway':30s} {'nodes':>6} {'intra':>6} {'neg':>5} "
                      f"{'rxns':>5} {'n/rxn':>7}  class")
        print(cls_header)
        print("-" * len(cls_header))
        counts: Counter = Counter()
        for r in results:
            for d in r["component_detail"]:
                ratio = "n/a" if d["nodes_per_reaction"] is None else f"{d['nodes_per_reaction']:.1f}"
                cls = d["classification"] or "unclassified (no reaction provenance)"
                counts[cls] += 1
                print(f"{r['pathway'][:30]:30s} {d['nodes']:6d} {d['intra_edges']:6d} "
                      f"{d['negative_edges']:5d} {d['reaction_stids']:5d} {ratio:>7}  {cls}")
        print("-" * len(cls_header))
        for cls, n in counts.most_common():
            held = sum(d["nodes"] for r in results for d in r["component_detail"]
                       if (d["classification"] or "unclassified (no reaction provenance)") == cls)
            print(f"  {cls}: {n} components holding {held} nodes")

    if args.json:
        args.json.write_text(json.dumps(
            {"source": source, "ratio_threshold": args.ratio_threshold, "pathways": results},
            indent=2))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
