#!/usr/bin/env python3
"""Classify loops at the REACTION level, where `precedingEvent` actually lives.

Three earlier formulations of this classifier disagreed wildly — 49% of SCCs
"real", then 18%, then 0.4% — because all three tried to test a REACTION-level
annotation (`precedingEvent`) against an ENTITY-level graph. The logic network
is bipartite (entity -> reaction -> entity), so most of its edges connect a
node to its own reaction and can never carry a reaction-to-reaction annotation.
The mapping is lossy in a way that silently biases the answer.

So build the reaction graph directly: reaction A -> reaction B when A outputs
an entity that B consumes. That is the graph the deck's rule is written about
(Rn, R_{n-i}), and `cache/reaction_connections.csv` labels exactly those pairs
as "Has Preceding Event" or "No Preceding Event".

Then, per reaction-level SCC:
    delete every link NOT curator-annotated; is it still strongly connected?
      yes -> the loop stands on curator annotation             REAL
      no  -> it exists only because of inferred links          ARTIFACT
and split REAL by whether a surviving link is a negative regulation (III) or
not (IV), and ARTIFACT by whether the entity doing the recycling is a
decomposition product (I) or a plain shared participant (II).
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

from loop_taxonomy import tarjan  # reuse the iterative Tarjan


def classify(d: Path):
    ctx = d / "node_reaction_context.csv"
    rconn = d / "cache" / "reaction_connections.csv"
    logic = d / "logic_network.csv"
    if not (ctx.exists() and rconn.exists() and logic.exists()):
        return None

    produced = collections.defaultdict(set)   # reaction -> entities it outputs
    consumed = collections.defaultdict(set)   # reaction -> entities it takes
    roles = collections.defaultdict(set)
    with ctx.open() as fh:
        for r in csv.DictReader(fh):
            node, rxn, role = r["context_node"], r["reaction_id"], r.get("role", "")
            roles[(node, rxn)].add(role)
            if role == "output":
                produced[rxn].add(node)
            else:
                consumed[rxn].add(node)

    annotated = set()
    with rconn.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("event_status") == "Has Preceding Event":
                annotated.add((r["preceding_reaction_id"], r["following_reaction_id"]))

    # negative edges, by the entity that carries them
    neg_nodes = set()
    node_kind = {}
    with logic.open() as fh:
        for e in csv.DictReader(fh):
            if e["pos_neg"] == "neg":
                neg_nodes.add(e["source_id"])
            node_kind[e["source_id"]] = e["edge_type"]

    by_entity = collections.defaultdict(set)
    for rxn, ents in produced.items():
        for e in ents:
            by_entity[e].add(rxn)

    adj = collections.defaultdict(list)
    link_entity = collections.defaultdict(set)
    reactions = set(produced) | set(consumed)
    for rxn, ents in consumed.items():
        for e in ents:
            for src in by_entity.get(e, ()):
                if src != rxn:
                    adj[src].append(rxn)
                    link_entity[(src, rxn)].add(e)

    comp, members = tarjan(reactions, adj)
    sccs = {c: m for c, m in members.items() if len(m) > 1}
    counts: collections.Counter = collections.Counter()
    for cid, mem in sccs.items():
        inside = set(mem)
        links = [(a, b) for a in inside for b in adj.get(a, ()) if b in inside]
        links = list(set(links))
        kept = [(a, b) for a, b in links if (a, b) in annotated]
        kept_adj = collections.defaultdict(list)
        for a, b in kept:
            kept_adj[a].append(b)
        _, km = tarjan(inside, kept_adj)
        survives = [m for m in km.values() if len(m) > 1]
        real = bool(survives)
        surv = {n for m in survives for n in m}

        if real:
            has_neg = any(any(x in neg_nodes for x in link_entity[(a, b)])
                          for a, b in kept if a in surv and b in surv)
            counts["III_negative_feedback" if has_neg else "IV_positive_feedback"] += 1
            counts[("III" if has_neg else "IV") + "_rxns"] += len(mem)
        else:
            # Type I: the recycling entity is a decomposition product, so it is
            # carried on a synthetic assembly/dissociation edge.
            ents = {x for a, b in links for x in link_entity[(a, b)]}
            obscure = any(node_kind.get(x) in ("assembly", "dissociation") for x in ents)
            counts["I_recycling_obscure" if obscure else "II_recycling_explicit"] += 1
            counts[("I" if obscure else "II") + "_rxns"] += len(mem)
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalog", type=Path, required=True)
    args = ap.parse_args()
    total: collections.Counter = collections.Counter()
    for d in sorted(p for p in args.catalog.iterdir() if p.is_dir()):
        got = classify(d)
        if got:
            total.update(got)
    types = ["I_recycling_obscure", "II_recycling_explicit",
             "III_negative_feedback", "IV_positive_feedback"]
    n = sum(total[t] for t in types)
    if not n:
        print("no reaction-level SCCs found")
        return
    print(f"reaction-level SCCs: {n}\n")
    print(f"{'type':<26}{'SCCs':>8}{'share':>8}{'reactions':>11}")
    for t in types:
        key = t.split("_")[0] + "_rxns"
        print(f"{t:<26}{total[t]:>8}{100*total[t]/n:>7.1f}%{total[key]:>11}")
    art = total[types[0]] + total[types[1]]
    real = total[types[2]] + total[types[3]]
    print(f"\n  ARTIFACT (break): {art:>5} ({100*art/n:.1f}%)")
    print(f"  REAL (model):     {real:>5} ({100*real/n:.1f}%)")


if __name__ == "__main__":
    main()
