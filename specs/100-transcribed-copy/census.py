"""specs/100 census: entities with a ROOT copy (no in-edge) and a PRODUCED copy (an
`output` in-edge) of the same stable id in one pathway, split by molecule type
(Neo4j referenceEntity class: RNA vs protein vs other), and how many root copies
are consumed by a reaction while their produced twin is not.
Run from the generator's poetry env with NEO4J_PASSWORD set.
Usage: python3 census.py BUILD_DIR"""
import csv, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, "/home/awright/gitroot/logic-network-generator")
from src.neo4j_connector import get_graph
build = sys.argv[1]
pairs = []   # (pathway, stid, n_root, n_produced)
for pw in sorted(p for p in os.listdir(build) if p.startswith("R-HSA-")):
    d = os.path.join(build, pw)
    if not os.path.exists(os.path.join(d, "nodes.csv")):
        continue
    ent = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(open(os.path.join(d, "nodes.csv")))
           if r["node_kind"] in ("simple_entity", "simple_complex")}
    has_in = set(); produced = set()
    for f in ("logic_network.csv", "boundary_edges.csv"):
        p = os.path.join(d, f)
        if not os.path.exists(p):
            continue
        for e in csv.DictReader(open(p)):
            has_in.add(e["target_id"])
            e["edge_type"] == "output" and produced.add(e["target_id"])
    by = defaultdict(lambda: [0, 0])
    for u, s in ent.items():
        if u not in has_in: by[s][0] += 1
        elif u in produced: by[s][1] += 1
    for s, (r, p) in by.items():
        if r and p:
            pairs.append((pw, s, r, p))
stids = sorted({s for _, s, _, _ in pairs})
cls = {x["s"]: x["c"] for x in get_graph().run(
    "UNWIND $s AS s MATCH (e:PhysicalEntity {stId: s}) OPTIONAL MATCH (e)-[:referenceEntity]->(r) "
    "RETURN s, coalesce(r.schemaClass, e.schemaClass) AS c", s=stids).data()}
c = Counter(); pw = defaultdict(set)
for p, s, r, n in pairs:
    k = cls.get(s, "?")
    c[k] += 1; pw[k].add(p)
for k, v in c.most_common():
    print(f"{k}: {v} entity pairs in {len(pw[k])} pathways")
rna = [(p, s) for p, s, _, _ in pairs if cls.get(s) == "ReferenceRNASequence"]
print("RNA pairs:", rna)
