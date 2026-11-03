#!/usr/bin/env python3
"""specs/132 step 1a: Reactome's curated cross-pathway handoffs, from Neo4j only.

For every catalog pathway P (bench/catalog_pathways.tsv): each pair (r, s) with r inside P (hasEvent*), s outside P,
r precededBy s (`(r)-[:precedingEvent]->(s)`), and an entity that is an output of s and an input, catalyst or regulator
of r (same stId). L2 marks pairs whose s lies inside a pathway that P, or a pathway inside P, encapsulates
(`hasEncapsulatedEvent`, the diagram's process nodes). Q is each direct parent pathway of s.

Writes TSV: pathway, r, s, q, entity, role, l2. No case data is read.
Usage: curated_links_132.py OUT.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as bm  # noqa: E402

QUERY = """
UNWIND $pathways AS pid
MATCH (P:Pathway {stId: pid})
MATCH (P)-[:hasEvent*]->(r:ReactionLikeEvent)-[:precedingEvent]->(s:ReactionLikeEvent)
WHERE NOT (P)-[:hasEvent*]->(s)
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
OPTIONAL MATCH (P)-[:hasEvent*0..]->(:Pathway)-[:hasEncapsulatedEvent]->(E:Pathway)-[:hasEvent*]->(s)
RETURN P.stId AS pathway, r.stId AS r, s.stId AS s, q.stId AS q, e.stId AS entity, role, count(E) > 0 AS l2
"""


def main():
    out = sys.argv[1]
    pathways = [line.split("\t")[0] for line in open(Path(__file__).resolve().parents[1] / "catalog_pathways.tsv")
                if line.startswith("R-HSA")]
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = []
    for p in pathways:      # one pathway per query keeps each bounded
        rows += graph.run(QUERY, pathways=[p]).data()
    rows = sorted({tuple(r[k] for k in ("pathway", "r", "s", "q", "entity", "role", "l2")) for r in rows})
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pathway", "r", "s", "q", "entity", "role", "l2"])
        w.writerows(rows)
    print(f"{len(pathways)} pathways; {len(rows)} rows; "
          f"{len({(r[0], r[1], r[2]) for r in rows})} handoff pairs; "
          f"{len({r[0] for r in rows})} pathways with any; {len({r[3] for r in rows})} source pathways Q")


if __name__ == "__main__":
    main()
