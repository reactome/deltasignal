#!/usr/bin/env python3
"""SCC edge-type breakdown for the logic-network catalog.

WHY
---
The giant SCCs in the generated logic networks (TP53 ~1250 nodes, WNT ~622) are
catalytic recycling cycles, not feedback. The loop-closing edges are `catalyst`
edges: a downstream reaction regenerates a catalyst that feeds an upstream one
(see project_loop_taxonomy_finding memory). If that hypothesis holds, then
removing `catalyst` edges (modeling a regenerated catalyst as a conserved-pool
modulator rather than a loop-carrying signal) should dissolve the giant SCCs.

This tool quantifies, per pathway:
  1. The nontrivial SCC structure (count, largest, % of nodes in loops).
  2. Edge-type composition of edges INSIDE the largest SCC.
  3. How the largest SCC shrinks when each edge_type is removed from the graph
     (the lever test): catalyst, regulator, assembly/dissociation.

This works at the logic-network CSV level (where `edge_type`/`pos_neg` live),
unlike loop_diagnostics.jl which works at the parsed-reaction (solver) level.

USAGE
-----
  python3 bench/analysis/scc_edge_breakdown.py              # 9 exp pathways
  python3 bench/analysis/scc_edge_breakdown.py TP53 WNT     # substring filter
  PATHWAY_CATALOG=/path/to/output python3 .../scc_edge_breakdown.py

Reads the catalog from $PATHWAY_CATALOG (default
~/gitroot/logic-network-generator/output).
"""
import csv, os, sys
from collections import defaultdict
from pathlib import Path

CATALOG = Path(os.environ.get(
    "PATHWAY_CATALOG", str(Path.home() / "gitroot" / "logic-network-generator" / "output")))

# Default 9-pathway experimental benchmark dirs (R-HSA-suffixed).
EXP_DIRS = [
    "Transcriptional_Regulation_by_TP53_R-HSA-3700989",
    "Mitotic_G1-G1_S_phases_R-HSA-453279",
    "S_Phase_R-HSA-69242",
    "Signaling_by_ERBB2_R-HSA-1227986",
    "Cell_Cycle_Checkpoints_R-HSA-69620",
    "Signaling_by_WNT_R-HSA-195721",
    "PIP3_activates_AKT_signaling_R-HSA-1257604",
    "Mitotic_Prophase_R-HSA-68875",
    "HDR_through_Homologous_Recombination_HRR_or_Single_Strand_Annealing_SSA__R-HSA-5693567",
]


def load_edges(dirname):
    """Return (edges, stid) where edges is a list of (src,tgt,edge_type,pos_neg)
    and stid maps uuid -> stable_id."""
    d = CATALOG / dirname
    edges = []
    with open(d / "logic_network.csv") as f:
        for r in csv.DictReader(f):
            edges.append((r["source_id"], r["target_id"],
                          r["edge_type"], r["pos_neg"]))
    stid = {}
    with open(d / "stid_to_uuid_mapping.csv") as f:
        for r in csv.DictReader(f):
            stid[r["uuid"]] = r["stable_id"]
    return edges, stid


def tarjan(adj, nodes):
    """Iterative Tarjan SCC. adj: dict node->list of nodes. Returns list of
    components (each a set)."""
    index = {}
    low = {}
    onstack = set()
    stack = []
    comps = []
    counter = [0]
    for s in nodes:
        if s in index:
            continue
        work = [(s, iter(adj.get(s, ())))]
        index[s] = low[s] = counter[0]; counter[0] += 1
        stack.append(s); onstack.add(s)
        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = counter[0]; counter[0] += 1
                    stack.append(w); onstack.add(w)
                    work.append((w, iter(adj.get(w, ()))))
                    advanced = True
                    break
                elif w in onstack:
                    low[v] = min(low[v], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                p = work[-1][0]
                low[p] = min(low[p], low[v])
            if low[v] == index[v]:
                comp = set()
                while True:
                    w = stack.pop(); onstack.discard(w); comp.add(w)
                    if w == v:
                        break
                comps.append(comp)
    return comps


def build_adj(edges, drop_types=frozenset()):
    adj = defaultdict(list)
    nodes = set()
    for s, t, et, pn in edges:
        nodes.add(s); nodes.add(t)
        if et in drop_types:
            continue
        adj[s].append(t)
    return adj, nodes


def largest_scc_size(edges, drop_types=frozenset()):
    adj, nodes = build_adj(edges, drop_types)
    comps = tarjan(adj, nodes)
    nontrivial = [c for c in comps if len(c) > 1]
    largest = max((len(c) for c in nontrivial), default=0)
    in_loops = sum(len(c) for c in nontrivial)
    return largest, len(nontrivial), in_loops, nontrivial


def analyze(dirname):
    edges, stid = load_edges(dirname)
    nnodes = len({s for s, *_ in edges} | {t for _, t, *_ in edges})
    largest, ncomp, in_loops, nontrivial = largest_scc_size(edges)
    print(f"\n### {dirname}")
    print(f"  edges={len(edges)} nodes={nnodes} "
          f"nontrivial_SCCs={ncomp} largest={largest} in_loops={in_loops}")

    if largest == 0:
        return
    big = max(nontrivial, key=len)
    # internal-edge edge_type breakdown
    bt = defaultdict(int)
    for s, t, et, pn in edges:
        if s in big and t in big:
            bt[(et, pn)] += 1
    print("  largest-SCC internal edges by (edge_type, sign):")
    for (et, pn), c in sorted(bt.items(), key=lambda x: -x[1]):
        print(f"      {et:14}/{pn:3} {c}")

    # The lever test: drop each edge_type set and see how the largest SCC shrinks
    scenarios = [
        ("drop catalyst", {"catalyst"}),
        ("drop catalyst+regulator", {"catalyst", "regulator"}),
        ("drop assembly+dissociation", {"assembly", "dissociation"}),
        ("drop catalyst+assembly+dissociation", {"catalyst", "assembly", "dissociation"}),
    ]
    print("  largest-SCC after removing edge types (the modeling lever):")
    for name, drop in scenarios:
        lg, nc, il, nontriv = largest_scc_size(edges, drop)
        # characterize the residual largest SCC: stid domination tells us
        # whether what survives is still a single-reaction recycling tangle
        # (artifact) or a multi-entity structure (candidate real feedback).
        resid = ""
        if lg > 1:
            rbig = max(nontriv, key=len)
            per = defaultdict(int)
            for n in rbig:
                per[stid.get(n, "?")] += 1
            top = sorted(per.items(), key=lambda x: -x[1])[:3]
            ratio = lg / len(per)
            resid = (f"  [{len(per)} stids, {ratio:.1f} nodes/stid; top: "
                     + ", ".join(f"{s}×{c}" for s, c in top) + "]")
        print(f"      {name:38} largest={lg:6}  total_in_loops={il}{resid}")


def main():
    args = sys.argv[1:]
    dirs = EXP_DIRS
    if args:
        dirs = [d for d in EXP_DIRS if any(a.lower() in d.lower() for a in args)]
        # also allow exact dir names not in EXP_DIRS
        for a in args:
            if (CATALOG / a).is_dir() and a not in dirs:
                dirs.append(a)
    for d in dirs:
        if not (CATALOG / d).is_dir():
            print(f"\n### {d}\n  MISSING")
            continue
        analyze(d)


if __name__ == "__main__":
    main()
