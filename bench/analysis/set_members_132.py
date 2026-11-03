#!/usr/bin/env python3
"""specs/132: the members of each handed-off entity that is a set (hasMember / hasCandidate, up to 3 levels), from
Neo4j only. The generator expands a set into its members, so the set's own id has no node; a handoff of a set is a
handoff of each member present on both sides.

Writes TSV: entity, member, entity_class. No case data is read.
Usage: set_members_132.py LINKS.tsv OUT.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as bm  # noqa: E402

QUERY = """
UNWIND $ids AS id
MATCH (e:PhysicalEntity {stId: id})
OPTIONAL MATCH (e)-[:hasMember|hasCandidate*1..3]->(m:PhysicalEntity)
RETURN e.stId AS entity, e.schemaClass AS cls, collect(DISTINCT m.stId) AS members
"""


def main():
    links, out = sys.argv[1:]
    from py2neo import Graph
    ids = sorted({r["entity"] for r in csv.DictReader(open(links), delimiter="\t")})
    g = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = g.run(QUERY, ids=ids).data()
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["entity", "member", "entity_class"])
        for r in sorted(rows, key=lambda r: r["entity"]):
            for m in sorted(r["members"]):
                w.writerow([r["entity"], m, r["cls"]])
    print(f"{len(rows)} entities; {sum(bool(r['members']) for r in rows)} sets; "
          f"{sum(len(r['members']) for r in rows)} members")


if __name__ == "__main__":
    main()
