"""specs/089 census: which ubiquitin-pass substrate leaves are the ligase.

For every ubiquitination (Ub input) in the catalog pathways, the generator's
substrate set is the input leaves whose gene is also in the output. This lists,
per reaction, the leaves the ligase-unit rule would drop:
  A  a leaf of the catalyst's curated activeUnit;
  B  (no activeUnit) a leaf of an output entity that is also a direct
     component of an input complex: the ligase released unchanged.
Either way a leaf is dropped only if its own stId is also an output leaf: a
ligase that is itself modified (autoubiquitination, Ub-BRAP) stays a substrate.
Run from the generator's poetry env with NEO4J_PASSWORD set.
"""
import csv, sys
sys.path.insert(0, "/home/awright/gitroot/logic-network-generator")
from src.neo4j_connector import get_graph
from src.logic_network_generator import _UBIQUITIN_STIDS

UB = list(_UBIQUITIN_STIDS)
pw = [r["id"] for r in csv.DictReader(
    (l for l in open(sys.argv[1]) if not l.startswith("#")), delimiter="\t")]
g = get_graph()
rows = g.run("""
UNWIND $pw AS p MATCH (:Pathway {stId:p})-[:hasEvent*]->(r:ReactionLikeEvent)-[:input]->(u)
WHERE u.stId IN $ub RETURN DISTINCT r.stId AS r, r.displayName AS name, collect(DISTINCT p) AS pws""",
    pw=pw, ub=UB).data()
LEAF = "OPTIONAL MATCH (x)-[:hasComponent|hasMember|hasCandidate*0..3]->(l:PhysicalEntity) OPTIONAL MATCH (l)-[:referenceEntity]->(re)"
out = csv.writer(sys.stdout, delimiter="\t")
out.writerow(["reaction", "name", "pathways", "substrate_leaves", "drop_A", "drop_B"])
for row in rows:
    r = row["r"]
    ins = g.run(f"MATCH (:ReactionLikeEvent {{stId:$r}})-[:input]->(x) WHERE NOT x.stId IN $ub {LEAF} "
                "RETURN l.stId AS s, l.displayName AS n, re.geneName AS gn", r=r, ub=UB).data()
    outs = g.run(f"MATCH (:ReactionLikeEvent {{stId:$r}})-[:output]->(x) WHERE NOT x.stId IN $ub {LEAF} "
                 "RETURN l.stId AS s, re.geneName AS gn", r=r, ub=UB).data()
    og = {gn for o in outs for gn in (o["gn"] or [])}
    ol = {o["s"] for o in outs if o["s"]}
    subst = {(i["s"], i["n"]) for i in ins if i["s"] and i["gn"] and og & set(i["gn"])}
    if not subst:
        continue
    au = {x["s"] for x in g.run(f"MATCH (:ReactionLikeEvent {{stId:$r}})-[:catalystActivity]->(ca)-[:activeUnit]->(x) {LEAF} "
                                "RETURN l.stId AS s", r=r).data() if x["s"]}
    rel = set()
    if not au:
        rel = {x["s"] for x in g.run(
            f"MATCH (rx:ReactionLikeEvent {{stId:$r}})-[:output]->(x)<-[:hasComponent]-(c)<-[:input]-(rx) {LEAF} "
            "RETURN l.stId AS s", r=r).data() if x["s"]}
    a = sorted(n for s, n in subst if s in au and s in ol)
    b = sorted(n for s, n in subst if s in rel and s in ol)
    out.writerow([r, row["name"], ",".join(row["pws"]), len(subst), "; ".join(a), "; ".join(b)])
