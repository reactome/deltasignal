# specs/086: an entity E whose STARTING copy S (no in-edge) is an input of a binding
# (or catalysed) reaction, and whose RELEASED copy P (same entity, made by a reaction) carries on
# downstream (has an out-edge), joined by a route S -> ... -> P of at most MAXR reactions
# through forms that contain E (containment table). The only way from S's level to P,
# and to whatever P feeds, is through the complex: a knockout of a complex partner cuts it.
# "Contain E" is decided by reference entity (Neo4j), because the route passes modified
# forms (p-S45-CTNNB1) that do not contain E's own stId.
# Usage (py2neo, NEO4J_PASSWORD): python3 count_released.py <build dir>
import csv, os, sys
from collections import Counter, defaultdict, deque
from py2neo import Graph
G = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"), auth=("neo4j", os.environ["NEO4J_PASSWORD"]))
_ref = {}
def refs(stids):
    need = [x for x in stids if x and x not in _ref]
    for x in need: _ref[x] = frozenset()
    if need:
        for r in G.run("""UNWIND $ids AS s MATCH (p:PhysicalEntity {stId:s})
            MATCH (p)-[:hasComponent|hasMember|hasCandidate*0..8]->(:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
            RETURN s, collect(DISTINCT re.dbId) AS r""", ids=need).data():
            _ref[r["s"]] = frozenset(r["r"])
    return _ref
MAXR = 10
b = sys.argv[1]
tot = Counter(); per = Counter(); ex = []
for pw in sorted(p for p in os.listdir(b) if p.startswith("R-HSA-")):
    d = f"{b}/{pw}"
    nodes = {r["uuid"]: r for r in csv.DictReader(open(d + "/nodes.csv"))}
    full = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
    cont = defaultdict(set)
    for r in csv.DictReader(open(d + "/containment.csv")): cont[r["stable_id"]].add(r["contains_stable_id"])
    cof = set()
    if os.path.exists(d + "/cofactors.csv"):
        cof = {list(r.values())[0] for r in csv.DictReader(open(d + "/cofactors.csv"))}
    inputs_of = defaultdict(set); rout = defaultdict(set); hasin = set(); hasout = set()
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if not os.path.exists(f"{d}/{f}"): continue
        for e in csv.DictReader(open(f"{d}/{f}")):
            s, t, et = e["source_id"], e["target_id"], e["edge_type"]
            if et == "recycle": continue
            hasin.add(t); hasout.add(s)
            if et in ("input", "catalyst"): inputs_of[s].add(t)
            if et == "output": rout[s].add(t)
    ent = lambda u: nodes.get(u, {}).get("diagram_entity_id", "")
    R = refs(sorted({ent(u) for u in nodes}))
    by_ent = defaultdict(list)
    for u in nodes: by_ent[ent(u)].append(u)
    found = set()
    for S, n in nodes.items():
        e = ent(S)
        if S in hasin or n["node_kind"] != "simple_entity" or not inputs_of.get(S) or e.startswith("R-ALL-") or e in cof:
            continue
        if len(R.get(e, ())) != 1: continue            # a single protein
        Ps = {u for u in by_ent[e] if u != S and u in hasout and any(u in rout[r] for r in rout)}
        if not Ps: continue
        # BFS S -> reaction -> output form (containing e) ... -> P
        q = deque([(S, 0)]); seen = {S}; hit = None
        while q and hit is None:
            x, k = q.popleft()
            if k >= MAXR: continue
            for rx in inputs_of.get(x, ()):
                for y in rout.get(rx, ()):
                    if y in Ps: hit = (y, k + 1); break
                    if y not in seen and R.get(e) and R.get(e) <= R.get(ent(y), frozenset()):
                        seen.add(y); q.append((y, k + 1))
                if hit: break
        if hit and hit[1] >= 2:
            found.add((e, hit[1]))
    if found:
        per[pw] = len(found); tot["entities"] += len(found)
        ex += [(pw, e, k) for e, k in sorted(found)]
print(dict(tot), "pathways", len(per))
for pw, n in per.most_common(12): print("  ", pw, n)
for x in ex[:30]: print("  ", x)
