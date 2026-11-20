#!/usr/bin/env python3
"""specs/167 FR-004: how much of an unselected genetic-perturbation transcript panel could the served networks score?

A structural count; no external data. Per pathway of a build:
  perturbable gene  a node with no in-edge (depletion ignored, as the bench's root rule) carries the gene
  target gene       an output of a curated gene-expression-regulated event is a node of the network
  pair              (perturbable gene, target gene), gene != target; reached when a root copy reaches a target node
Approximate on purpose: the bench's pin protocol (root_cycle) and readout rules are not applied.

Usage: transcript_feasibility.py BUILD_DIR OUT_JSON      (needs Neo4j for gene names and expression events)
"""
import csv
import json
import sys
from collections import defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))

GENES = """MATCH (e:EntityWithAccessionedSequence)-[:referenceEntity]->(re:ReferenceGeneProduct)
WHERE e.speciesName = 'Homo sapiens' AND re.geneName IS NOT NULL RETURN e.stId AS s, re.geneName[0] AS g"""
TARGETS = """MATCH (r:ReactionLikeEvent)-[:regulatedBy]->(reg)
WHERE r.speciesName = 'Homo sapiens' AND (reg:PositiveGeneExpressionRegulation OR reg:NegativeGeneExpressionRegulation)
MATCH (r)-[:output]->(o)-[:referenceEntity]->(tre) WHERE tre.geneName IS NOT NULL
RETURN DISTINCT o.stId AS s, tre.geneName[0] AS g"""


def pathway_pairs(edges, nodes, gene_of, target_nodes):
    """edges: (source, target, edge_type); nodes: (uuid, member leaf stIds); gene_of: leaf stId -> gene;
    target_nodes: gene -> uuids. Returns (all pairs, reached pairs, perturbable genes)."""
    adj, fed = defaultdict(list), set()
    for s, t, et in edges:
        adj[s].append(t)
        if et != "depletion":
            fed.add(t)
    roots = defaultdict(set)
    for uuid, leaves in nodes:
        if uuid in fed:
            continue
        for leaf in leaves:
            if leaf in gene_of:
                roots[gene_of[leaf]].add(uuid)
    pairs, reached = set(), set()
    for gene, rs in roots.items():
        seen, q = set(rs), deque(rs)
        while q:
            for w in adj.get(q.popleft(), ()):
                if w not in seen:
                    seen.add(w)
                    q.append(w)
        for t, us in target_nodes.items():
            if t == gene or not us:
                continue
            pairs.add((gene, t))
            if us & seen:
                reached.add((gene, t))
    return pairs, reached, set(roots)


def main():
    import benchmark_vs_mpbiopath as bm
    from py2neo import Graph
    build, out_json = Path(sys.argv[1]), Path(sys.argv[2])
    g = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    gene_of = {r["s"]: r["g"] for r in g.run(GENES)}
    targets = defaultdict(set)
    for r in g.run(TARGETS):
        targets[r["s"]].add(r["g"])
    release = g.run("MATCH (d:DBInfo) RETURN d.releaseNumber AS v").evaluate()
    all_pairs, all_reached, pert, tgt, rows = set(), set(), set(), set(), []
    for pdir in sorted(p for p in build.iterdir() if (p / "logic_network.csv").exists()):
        s2u = bm.load_stid_to_uuids(pdir)
        tnodes = defaultdict(set)
        for s, gs in targets.items():
            for u in s2u.get(s, []):
                for x in gs:
                    tnodes[x].add(u)
        edges = [(r["source_id"], r["target_id"], r["edge_type"])
                 for r in csv.DictReader(open(pdir / "logic_network.csv"))]
        nodes = [(r["uuid"], r["member_leaves"].split("|")) for r in csv.DictReader(open(pdir / "nodes.csv"))]
        pairs, reached, roots = pathway_pairs(edges, nodes, gene_of, tnodes)
        all_pairs |= pairs
        all_reached |= reached
        pert |= roots
        tgt |= set(tnodes)
        rows.append({"pathway": pdir.name, "perturbable_genes": len(roots), "target_genes": len(tnodes),
                     "pairs": len(pairs), "reached": len(reached)})
    res = {"build": build.name, "reactome_release": release,
           "expression_event_outputs": len(targets), "expression_target_genes": len({x for v in targets.values() for x in v}),
           "pathways": len(rows), "pathways_with_a_pair": sum(1 for r in rows if r["pairs"]),
           "pathways_with_a_reached_pair": sum(1 for r in rows if r["reached"]),
           "perturbable_genes": len(pert), "target_genes_in_networks": len(tgt),
           "unique_pairs": len(all_pairs), "unique_pairs_reached": len(all_reached),
           "by_pathway": sorted(rows, key=lambda r: -r["reached"])}
    out_json.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "by_pathway"}, indent=1))
    for r in res["by_pathway"][:10]:
        print(r)


if __name__ == "__main__":
    main()
