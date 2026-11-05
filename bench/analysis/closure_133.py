#!/usr/bin/env python3
"""specs/133 step 1: follow Reactome's curated cross-pathway handoffs back from the scored pathways, as far as they go.

Round 0 is the scored pathways (bench/catalog_pathways.tsv). In each round, every pathway X added in the previous round
receives handoffs: pairs (r in X, s outside X, r precededBy s) with an entity that s outputs and r takes as input,
catalyst or regulator (specs/132's L1 query). s's source pathway is its lowest ancestor with a diagram. Handoffs of an
entity in the extracellular region (or a set with any member there) are kept in the table, flagged, and add no
pathway. A source pathway not seen before joins the next round. Stops when a round adds nothing.

This is the superset the builder can draw from: the builder (context_build.py --transitive) keeps, for each scored
pathway, only the tree of handoffs that reach roots it keeps. Reads Neo4j only, no case data.

Writes OUT/links.tsv (pathway, r, s, q, entity, role, l2, extracellular), OUT/s_q.tsv (s, q_diagram, name) and
OUT/closure.tsv (pathway, round, name). A preceding reaction with two diagram pathways at the same least depth takes the
smaller stId, and the tie is printed. Usage: closure_133.py OUT_DIR
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as bm  # noqa: E402

LINKS = """
MATCH (P:Pathway {stId: $pid})
MATCH (P)-[:hasEvent*1..12]->(r:ReactionLikeEvent)-[:precedingEvent]->(s:ReactionLikeEvent)
WHERE NOT (P)-[:hasEvent*1..12]->(s) AND s.speciesName = 'Homo sapiens'
WITH DISTINCT P, r, s
MATCH (s)-[:output]->(e:PhysicalEntity)
WITH P, r, s, e
OPTIONAL MATCH (r)-[:input]->(ei) WHERE ei = e
OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(ec) WHERE ec = e
OPTIONAL MATCH (r)-[:regulatedBy]->(reg:Regulation)-[:regulator]->(er) WHERE er = e
WITH P, r, s, e,
     CASE WHEN ei IS NOT NULL THEN 'input' WHEN ec IS NOT NULL THEN 'catalyst'
          WHEN er IS NOT NULL THEN CASE WHEN reg:NegativeRegulation THEN 'regulator_neg' ELSE 'regulator_pos' END
          ELSE null END AS role
WHERE role IS NOT NULL
MATCH (q:Pathway)-[:hasEvent]->(s)
OPTIONAL MATCH (P)-[:hasEvent*0..12]->(:Pathway)-[:hasEncapsulatedEvent]->(E:Pathway)-[:hasEvent*1..12]->(s)
RETURN P.stId AS pathway, r.stId AS r, s.stId AS s, q.stId AS q, e.stId AS entity, role, count(E) > 0 AS l2
"""

# every diagram pathway above s with its least depth; pick_diagram chooses the lowest, ties to the smaller stId (review
# 2026-11-05: `collect(d)[0]` after `ORDER BY L` left ties to Neo4j's row order)
DIAGRAM = """
UNWIND $ids AS sid MATCH (s {stId: sid})
MATCH path = (d:Pathway)-[:hasEvent*1..12]->(s) WHERE d.hasDiagram = true
WITH sid, d, min(length(path)) AS L
RETURN sid AS s, d.stId AS q_diagram, d.displayName AS name, L
"""

# GO:0005576; a set counts as extracellular when it or any member (hasMember|hasCandidate, three levels) is
EXTRACELLULAR = """
UNWIND $ids AS eid MATCH (e:PhysicalEntity {stId: eid})
OPTIONAL MATCH (e)-[:hasMember|hasCandidate*0..3]->(m:PhysicalEntity)-[:compartment]->(c:Compartment)
WHERE c.accession = '0005576'
RETURN eid AS entity, count(c) > 0 AS extracellular
"""


def scored_pathways():
    return [line.split("\t")[0] for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv")
            if line.startswith("R-HSA")]


def pick_diagram(rows):
    """From DIAGRAM rows, each s's source pathway: the diagram pathway at the least depth above it, ties broken by
    stId. Returns ({s: q_diagram}, {q_diagram: name}, {s: [tied q_diagram, ...]} for every s with a tie)."""
    best, names, ties = {}, {}, {}
    by_s = {}
    for r in rows:
        by_s.setdefault(r["s"], []).append(r)
        names[r["q_diagram"]] = r["name"]
    for sid, rs in by_s.items():
        rs = sorted(rs, key=lambda r: (r["L"], r["q_diagram"]))
        best[sid] = rs[0]["q_diagram"]
        tied = [r["q_diagram"] for r in rs if r["L"] == rs[0]["L"]]
        if len(tied) > 1:
            ties[sid] = tied
    return best, names, ties


def closure(graph, seeds):
    """(link rows, {s: q_diagram}, {pathway: name}, {pathway: round}, {s: tied diagram pathways})."""
    rnd = {p: 0 for p in seeds}
    frontier, links, s_q, extra, names, ties, k = list(seeds), [], {}, {}, {}, {}, 0
    while frontier:
        new_rows = []
        for p in frontier:             # one receiving pathway per query keeps each bounded
            new_rows += graph.run(LINKS, pid=p).data()
        new_s = sorted({r["s"] for r in new_rows} - set(s_q))
        best, nm, tied = pick_diagram(graph.run(DIAGRAM, ids=new_s).data())
        s_q.update(best); names.update(nm); ties.update(tied)
        new_e = sorted({r["entity"] for r in new_rows} - set(extra))
        for r in graph.run(EXTRACELLULAR, ids=new_e).data():
            extra[r["entity"]] = r["extracellular"]
        for r in new_rows:
            r["extracellular"] = extra.get(r["entity"], False)
        links += new_rows
        k += 1
        added = sorted({s_q[r["s"]] for r in new_rows if r["s"] in s_q and not r["extracellular"]} - set(rnd))
        for q in added:
            rnd[q] = k
        frontier = added
    return links, s_q, names, rnd, ties


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    links, s_q, names, rnd, ties = closure(graph, scored_pathways())
    cols = ("pathway", "r", "s", "q", "entity", "role", "l2", "extracellular")
    rows = sorted({tuple(r[c] for c in cols) for r in links})
    with open(out / "links.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(cols); w.writerows(rows)
    with open(out / "s_q.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(("s", "q_diagram", "name"))
        w.writerows(sorted((s, q, names.get(q, "")) for s, q in s_q.items()))
    with open(out / "closure.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(("pathway", "round", "name"))
        w.writerows(sorted(((p, k, names.get(p, "")) for p, k in rnd.items()), key=lambda t: (t[1], t[0])))
    by = {}
    for k in rnd.values():
        by[k] = by.get(k, 0) + 1
    print(f"{len(rows)} handoff rows ({sum(1 for r in rows if r[7])} extracellular); pathways by round: "
          + ", ".join(f"{k}: {n}" for k, n in sorted(by.items())))
    print(f"{len(ties)} preceding reactions with tied source pathways (lowest diagram at equal depth; the smaller "
          "stId is taken)" + "".join(f"\n  {s}: {', '.join(t)}" for s, t in sorted(ties.items())))


if __name__ == "__main__":
    main()
