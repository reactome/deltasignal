#!/usr/bin/env python3
"""Audit generator connectivity: for each (perturbation gene -> key-output) pair
that the ground truth expects to CHANGE, is the key-output directed-reachable
from the perturbed gene in the generated network?

A structurally sound network should connect most expected-change pairs. A high
"unreachable" fraction means the generator is dropping perturbation->readout
connectivity, which caps accuracy regardless of the propagator (every unreachable
expected-change pair is a forced NORMAL misprediction).

For unreachable pairs it classifies WHY, which says whether connectivity work
could recover them:
  reverse        readout is UPSTREAM of the perturbation (directionality)
  shared_ancestor readout co-regulated via a common driver (BN-style coupling)
  diff_component readout in a different weakly-connected component (missing edges
                 -- the most directly fixable by connectivity work)
  zigzag         same component but only via alternating-direction path (no clean
                 causal story; connectivity within the pathway won't cleanly help)
  unresolved     gene or readout not in the network at all

Usage:
  python bench/analysis/connectivity_audit.py [--ground-truth curator|experimental]
      [--pathway-list FILE] [--max-edges N] [--only NAME[,NAME...]]
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as b  # noqa: E402


def load_graph(pathway_dir: Path):
    fwd = defaultdict(list)
    rev = defaultdict(list)
    un = defaultdict(set)
    for r in csv.DictReader(open(pathway_dir / "logic_network.csv")):
        s, t = r["source_id"], r["target_id"]
        fwd[s].append(t)
        rev[t].append(s)
        un[s].add(t)
        un[t].add(s)
    return fwd, rev, un


def reach(adj, src):
    seen = set(src)
    fr = list(src)
    while fr:
        nx = []
        for u in fr:
            for v in adj.get(u, ()):
                if v not in seen:
                    seen.add(v)
                    nx.append(v)
        fr = nx
    return seen


def components(un, nodes):
    seen = set()
    cid = {}
    c = 0
    for n in nodes:
        if n in seen:
            continue
        stack = [n]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            cid[u] = c
            stack.extend(un[u] - seen)
        c += 1
    return cid


def classify(gu, ku, fwd, rev, cid):
    """Why is ku not forward-reachable from gu? gu/ku are uuid lists."""
    gu, ku = set(gu), set(ku)
    if reach(rev, gu) & ku:
        return "reverse"
    anc = reach(rev, gu)
    if reach(fwd, anc) & ku:
        return "shared_ancestor"
    if not ({cid.get(u) for u in gu} & {cid.get(u) for u in ku}):
        return "diff_component"
    return "zigzag"


def audit_pathway(pid, pname, gt):
    pdir = b.find_pathway_dir(pid)
    if pdir is None:
        return None
    cur = b.load_curator(pname, ground_truth=gt)
    if cur is None:
        return None
    header, rows = cur
    s2u = b.load_stid_to_uuids(pdir)
    proxies = b.load_entity_reaction_proxies(pdir)
    perts = b.parse_perturbation_columns(header)
    genes = sorted({g for g, _, _ in perts})
    g2s = b.neo4j_gene_to_stids(genes)
    g2u = {g: [u for sid in g2s.get(g, []) for u in s2u.get(sid, [])] for g in genes}
    d2s = b.neo4j_dbid_to_stid({str(r["key_output"]) for r in rows})

    def ko_uuids(ko):
        sid = d2s.get(ko) or f"R-HSA-{ko}"
        return s2u.get(sid, []) or proxies.get(sid, [])

    fwd, rev, un = load_graph(pdir)
    nodes = set(fwd) | set(rev)
    cid = components(un, nodes)

    tally = defaultdict(int)
    n_change = 0
    reachable = 0
    for g, direction, col in perts:
        gu = g2u.get(g, [])
        for r in rows:
            try:
                exp = int(r[col])
            except (KeyError, ValueError):
                continue
            if exp == 1 or exp == -999:  # only pairs expected to CHANGE
                continue
            n_change += 1
            ku = ko_uuids(str(r["key_output"]))
            if not gu or not ku:
                tally["unresolved"] += 1
                continue
            if reach(fwd, set(gu)) & set(ku):
                reachable += 1
            else:
                tally[classify(gu, ku, fwd, rev, cid)] += 1
    return {
        "name": pname, "n_change": n_change, "reachable": reachable,
        "n_components": len(set(cid.values())), "n_nodes": len(nodes),
        "tally": dict(tally),
    }


def main():
    args = sys.argv[1:]
    gt = "curator"
    pathway_list = str(b.MPBIO_ROOT / "pathway_list.tsv")
    max_edges = 200000
    only = None
    for i, a in enumerate(args):
        if a == "--ground-truth":
            gt = args[i + 1]
        elif a == "--pathway-list":
            pathway_list = args[i + 1]
        elif a == "--max-edges":
            max_edges = int(args[i + 1])
        elif a == "--only":
            only = set(args[i + 1].split(","))

    with open(pathway_list) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ii, ni = hdr.index("pathway_id"), hdr.index("pathway_name")
        pathways = []
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) > max(ii, ni):
                pathways.append((c[ii], c[ni]))

    tot_change = tot_reach = 0
    agg = defaultdict(int)
    print(f"{'pathway':52} {'chg-pairs':>9} {'reach%':>7} {'diffcomp':>8} "
          f"{'shAnc':>6} {'rev':>5} {'zig':>5} {'unres':>6} {'comps':>6}")
    for pid, pname in pathways:
        if only and pname not in only:
            continue
        n = b.network_edge_count(pid)
        if n < 0 or n > max_edges:
            continue
        res = audit_pathway(pid, pname, gt)
        if res is None or res["n_change"] == 0:
            continue
        t = res["tally"]
        pct = 100 * res["reachable"] / res["n_change"]
        print(f"{res['name'][:52]:52} {res['n_change']:>9} {pct:>6.1f}% "
              f"{t.get('diff_component',0):>8} {t.get('shared_ancestor',0):>6} "
              f"{t.get('reverse',0):>5} {t.get('zigzag',0):>5} "
              f"{t.get('unresolved',0):>6} {res['n_components']:>6}")
        tot_change += res["n_change"]
        tot_reach += res["reachable"]
        for k, v in t.items():
            agg[k] += v
    if tot_change:
        print(f"\nTOTAL expected-change pairs: {tot_change}")
        print(f"  directed-reachable (connectivity OK): {tot_reach} = {100*tot_reach/tot_change:.1f}%")
        print(f"  UNREACHABLE breakdown: {dict(agg)}")
        fixable = agg.get('diff_component', 0) + agg.get('reverse', 0)
        print(f"  potentially connectivity-recoverable (diff_component+reverse): "
              f"{fixable} = {100*fixable/tot_change:.1f}% of all change-pairs")


if __name__ == "__main__":
    main()
