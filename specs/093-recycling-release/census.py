"""specs/093 census: open carrier-loop releases that move the carrier across
compartments (a receptor recycled from an endosome to the plasma membrane),
against those within one compartment (a partner letting go).

For every open loop row (free form F, release reaction R) in a build, the
compartments of F and of R's inputs, read from Neo4j (the build's nodes.csv
carries none). A release is a translocation when F's compartment is not among
its inputs' compartments. Counted per distinct (pathway, F entity, R reaction).
Run from the generator's poetry env with NEO4J_PASSWORD set.
Usage: python3 census.py <build dir>
"""
import csv, os, sys
from collections import Counter
sys.path.insert(0, "/home/awright/gitroot/logic-network-generator")
from src.neo4j_connector import get_graph

build = sys.argv[1]
base = lambda s: s.split("::")[0]
pairs = set()
for pw in sorted(p for p in os.listdir(build) if p.startswith("R-HSA-")):
    d = os.path.join(build, pw)
    if not os.path.exists(os.path.join(d, "carrier_loops.csv")):
        continue
    ent = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(open(os.path.join(d, "nodes.csv")))}
    for r in csv.DictReader(open(os.path.join(d, "carrier_loops.csv"))):
        if r["closed"] != "True":
            pairs.add((pw, base(r["free_stid"]), base(ent.get(r["loop_source_uuid"], ""))))
g = get_graph()
rows = g.run("""UNWIND $p AS p
MATCH (f:PhysicalEntity {stId: p[1]})-[:compartment]->(cf)
MATCH (r:ReactionLikeEvent {stId: p[2]})-[:input]->(i)-[:compartment]->(ci)
RETURN p[0] AS pw, p[1] AS f, p[2] AS r, f.displayName AS fn, r.displayName AS rn,
       collect(DISTINCT cf.name) AS cf, collect(DISTINCT ci.name) AS ci""",
             p=[list(x) for x in sorted(pairs)]).data()
tot = Counter(); per = Counter()
for x in rows:
    kind = "translocation" if not set(x["cf"]) & set(x["ci"]) else "same compartment"
    tot[kind] += 1
    if kind == "translocation":
        per[x["pw"]] += 1
        print(x["pw"], x["f"], x["fn"][:60], "|", x["r"], x["rn"][:70], "|", x["ci"], "->", x["cf"], sep="\t")
print(dict(tot), "of", len(pairs), "distinct (pathway, free form, release) pairs")
print(dict(per))
