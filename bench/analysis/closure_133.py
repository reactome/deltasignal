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
OUT/closure.tsv (pathway, round, name). Usage: closure_133.py OUT_DIR
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

DIAGRAM = """
UNWIND $ids AS sid MATCH (s {stId: sid})
MATCH path = (d:Pathway)-[:hasEvent*1..12]->(s) WHERE d.hasDiagram = true
WITH sid, d, length(path) AS L ORDER BY L
WITH sid, collect(d)[0] AS d
RETURN sid AS s, d.stId AS q_diagram, d.displayName AS name
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


def closure(graph, seeds):
    """(link rows, {s: (q_diagram, name)}, {pathway: round})."""
    rnd = {p: 0 for p in seeds}
    frontier, links, s_q, extra, names, k = list(seeds), [], {}, {}, {}, 0
    while frontier:
        new_rows = []
        for p in frontier:             # one receiving pathway per query keeps each bounded
            new_rows += graph.run(LINKS, pid=p).data()
        new_s = sorted({r["s"] for r in new_rows} - set(s_q))
        for r in graph.run(DIAGRAM, ids=new_s).data():
            s_q[r["s"]] = r["q_diagram"]
            names[r["q_diagram"]] = r["name"]
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
    return links, s_q, names, rnd


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    links, s_q, names, rnd = closure(graph, scored_pathways())
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


if __name__ == "__main__":
    main()
