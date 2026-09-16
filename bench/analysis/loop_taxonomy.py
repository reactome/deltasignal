#!/usr/bin/env python3
"""Classify every SCC in a generated network into the four curator loop types.

From "Reactome Pathway Loops Revisited". Loops are not one thing, and the four
types want opposite treatment — two are artifacts to break, two are real
biology to model:

  Type I   recycling, obscure   an artifact of decomposing a root complex or
                                set into its parts, so a component appears on
                                both sides of the same event. Break.
  Type II  recycling, explicit  the input of a preceding event is the output of
                                a following one in the same chain, with NO
                                curator `precedingEvent` backing the closing
                                direction. Break.
  Type III negative feedback    real, and closed by a NEGATIVE regulator.
                                Model, do not break.
  Type IV  positive feedback    real, and closed by curator-annotated
                                precedingEvent links. Model, do not break.

The discriminator that makes this possible is `event_status` in
`cache/reaction_connections.csv`: 4,580 connections across the catalog are
"Has Preceding Event" (a curator said so) and 1,862 are "No Preceding Event"
(the pipeline inferred it). A cycle closed only by inferred links is an
artifact; one closed by an annotated link is biology.

Counting them is the prerequisite for handling them. Nothing here changes a
solve — this reports what is in the networks.
"""

from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path


def tarjan(nodes, adj):
    """Iterative Tarjan; returns node -> component id and component members."""
    index = {}
    low = {}
    on = set()
    stack = []
    comp = {}
    members = collections.defaultdict(list)
    counter = [0]
    cid = [0]
    for root in nodes:
        if root in index:
            continue
        work = [(root, 0)]
        while work:
            v, pi = work.pop()
            if pi == 0:
                index[v] = low[v] = counter[0]
                counter[0] += 1
                stack.append(v)
                on.add(v)
            recurse = False
            neighbours = adj.get(v, ())
            for i in range(pi, len(neighbours)):
                w = neighbours[i]
                if w not in index:
                    work.append((v, i + 1))
                    work.append((w, 0))
                    recurse = True
                    break
                if w in on:
                    low[v] = min(low[v], index[w])
            if recurse:
                continue
            if low[v] == index[v]:
                while True:
                    w = stack.pop()
                    on.discard(w)
                    comp[w] = cid[0]
                    members[cid[0]].append(w)
                    if w == v:
                        break
                cid[0] += 1
            if work:
                u = work[-1][0]
                low[u] = min(low[u], low[v])
    return comp, members


def classify_pathway(d: Path):
    logic = d / "logic_network.csv"
    ctx = d / "node_reaction_context.csv"
    rconn = d / "cache" / "reaction_connections.csv"
    if not logic.exists():
        return None

    edges = []
    adj = collections.defaultdict(list)
    nodes = set()
    with logic.open() as fh:
        for e in csv.DictReader(fh):
            edges.append(e)
            adj[e["source_id"]].append(e["target_id"])
            nodes.update((e["source_id"], e["target_id"]))

    comp, members = tarjan(nodes, adj)
    sccs = {c: m for c, m in members.items() if len(m) > 1}
    # self-loops are cycles too
    for e in edges:
        if e["source_id"] == e["target_id"] and comp[e["source_id"]] not in sccs:
            sccs[comp[e["source_id"]]] = [e["source_id"]]
    if not sccs:
        return collections.Counter(), {}

    # which reaction each node belongs to, and whether a connection is annotated
    node_rxn = collections.defaultdict(set)
    if ctx.exists():
        with ctx.open() as fh:
            for r in csv.DictReader(fh):
                node_rxn[r["context_node"]].add(r["reaction_id"])
    annotated = set()
    inferred = set()
    if rconn.exists():
        with rconn.open() as fh:
            for r in csv.DictReader(fh):
                pair = (r["preceding_reaction_id"], r["following_reaction_id"])
                if r.get("event_status") == "Has Preceding Event":
                    annotated.add(pair)
                else:
                    inferred.add(pair)

    counts: collections.Counter = collections.Counter()
    detail = {}
    for cid, mem in sccs.items():
        inside = set(mem)
        internal = [e for e in edges
                    if e["source_id"] in inside and e["target_id"] in inside]
        if not internal:
            continue
        neg = [e for e in internal if e["pos_neg"] == "neg"]
        regulator = [e for e in internal if e["edge_type"] == "regulator"]
        boundary = [e for e in internal
                    if e["edge_type"] in ("assembly", "dissociation")]

        # ORDER-INDEPENDENT TEST. Which edge "closes" a cycle depends on
        # where a DFS starts, and that ambiguity gave the same 1,127-node
        # component two different types in two pathways. So do not pick a
        # closing edge at all. Instead delete every edge that is NOT backed by
        # a curator `precedingEvent` and ask whether the component is STILL
        # strongly connected:
        #
        #   still cyclic -> the loop stands on curator-annotated links. Real.
        #   falls apart  -> the loop exists only because of inferred links.
        #                   Artifact.
        #
        # This is the deck's rule stated as a property of the component rather
        # than of an arbitrarily chosen edge.
        def edge_is_annotated(e):
            for a in node_rxn.get(e["source_id"], ()):
                for b in node_rxn.get(e["target_id"], ()):
                    if a != b and ((a, b) in annotated or (b, a) in annotated):
                        return True
            return False

        kept = [e for e in internal if edge_is_annotated(e)]
        kept_adj = collections.defaultdict(list)
        for e in kept:
            kept_adj[e["source_id"]].append(e["target_id"])
        _, kept_members = tarjan(inside, kept_adj)
        survives = [m for m in kept_members.values() if len(m) > 1]
        # a self-loop on an annotated edge also survives
        if not survives and any(e["source_id"] == e["target_id"] for e in kept):
            survives = [[kept[0]["source_id"]]]

        backed = bool(survives)
        surviving_nodes = set()
        for m in survives:
            surviving_nodes.update(m)
        neg_in_loop = [e for e in kept
                       if e["pos_neg"] == "neg"
                       and e["source_id"] in surviving_nodes
                       and e["target_id"] in surviving_nodes]

        if backed and neg_in_loop:
            t = "III_negative_feedback"
        elif backed:
            t = "IV_positive_feedback"
        elif boundary:
            t = "I_recycling_obscure"
        else:
            t = "II_recycling_explicit"
        counts[t] += 1
        counts[t + "_nodes"] += len(mem)
        detail[cid] = (t, len(mem))
    return counts, detail


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()

    total: collections.Counter = collections.Counter()
    biggest = []
    for d in sorted(p for p in args.catalog.iterdir() if p.is_dir()):
        got = classify_pathway(d)
        if not got:
            continue
        counts, detail = got
        total.update(counts)
        for cid, (t, n) in detail.items():
            biggest.append((n, t, d.name))

    types = ["I_recycling_obscure", "II_recycling_explicit",
             "III_negative_feedback", "IV_positive_feedback"]
    n_scc = sum(total[t] for t in types)
    print(f"SCCs classified: {n_scc}\n")
    print(f"{'type':<26}{'SCCs':>8}{'share':>8}{'nodes':>9}")
    for t in types:
        if not total[t]:
            continue
        print(f"{t:<26}{total[t]:>8}{100*total[t]/n_scc:>7.1f}%{total[t+'_nodes']:>9}")
    artifact = total["I_recycling_obscure"] + total["II_recycling_explicit"]
    real = total["III_negative_feedback"] + total["IV_positive_feedback"]
    print(f"\n  ARTIFACT (break):  {artifact:>6} SCCs")
    print(f"  REAL (model):      {real:>6} SCCs")
    biggest.sort(reverse=True)
    print(f"\nlargest {args.top} SCCs:")
    for n, t, pw in biggest[:args.top]:
        print(f"   {n:>6} nodes  {t:<24} {pw[:40]}")


if __name__ == "__main__":
    main()
