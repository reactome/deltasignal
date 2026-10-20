# specs/081: count recycle candidates in a catalog build.
# A candidate is a protein copy P that a reaction produces (an `output` in-edge),
# that leads nowhere (no out-edge) and is not a dissociation sink, together with
# a starting copy B of the same entity (diagram_entity_id) that feeds a root complex
# C by an `assembly` edge, and a reaction route C -> ... -> P (inputs to outputs,
# at most MAXR reactions) through forms that carry P's protein leaves. That route
# plus the recycle edge P -> B is the cycle; its first form after C is the modified state M.
# "Carries P's protein" is decided by reference entity (the gene product), read from
# Neo4j, because a modified form (p-RBL2) is a different entity from the bare one.
# Usage (needs py2neo, NEO4J_PASSWORD): python3 count_recycle.py <build dir> [pathway ...]
import csv, os, sys
from py2neo import Graph
G = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"), auth=("neo4j", os.environ["NEO4J_PASSWORD"]))
_ref = {}
def refs(stids):
    need = [x for x in stids if x not in _ref]
    if need:
        for x in need: _ref[x] = frozenset()
        for r in G.run("""UNWIND $ids AS s MATCH (p:PhysicalEntity {stId:s})
            MATCH (p)-[:hasComponent|hasMember|hasCandidate*0..8]->(:PhysicalEntity)-[:referenceEntity]->(re:ReferenceEntity)
            RETURN s, collect(DISTINCT re.dbId) AS r""", ids=need).data():
            _ref[r["s"]] = frozenset(r["r"])
    return _ref
from collections import defaultdict, deque
MAXR = 4
build = sys.argv[1]
pws = sys.argv[2:] or sorted(p for p in os.listdir(build) if p.startswith("R-HSA-"))
total = 0; per = {}
for pw in pws:
    d = f"{build}/{pw}"
    nodes = {r["uuid"]: r for r in csv.DictReader(open(d + "/nodes.csv"))}
    E = []
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if os.path.exists(f"{d}/{f}"): E += list(csv.DictReader(open(f"{d}/{f}")))
    out_of = defaultdict(list); inputs_of = defaultdict(list); outputs_of = defaultdict(list)
    produced = set(); hasout = set(); assembly = defaultdict(set)
    for e in E:
        s, t, et = e["source_id"], e["target_id"], e["edge_type"]
        hasout.add(s)
        if et == "output": produced.add(t); outputs_of[s].append(t)
        if et == "input": inputs_of[s].append(t)   # entity -> reactions it is an input of
        if et == "assembly": assembly[s].add(t)
    ent_ids = {n["diagram_entity_id"] for n in nodes.values()}
    R = refs(sorted(ent_ids))
    leaves = lambda u: R.get(nodes.get(u, {}).get("diagram_entity_id", ""), frozenset())
    cof = set()
    if os.path.exists(d + "/cofactors.csv"):
        cof = {list(r.values())[0] for r in csv.DictReader(open(d + "/cofactors.csv"))}
    by_ent = defaultdict(list)
    for u, n in nodes.items(): by_ent[n["diagram_entity_id"]].append(u)
    found = []
    for P, n in nodes.items():
        if P not in produced or P in hasout or n["node_kind"] == "dissociation_sink": continue
        if n["diagram_entity_id"].startswith("R-ALL-") or n["diagram_entity_id"] in cof: continue  # proteins only
        lp = leaves(P)
        if not lp: continue
        for B in by_ent[n["diagram_entity_id"]]:
            if B == P or B not in assembly: continue
            for C in assembly[B]:
                # BFS over reactions from C to P through forms carrying P's leaves
                q = deque([(C, [])]); seen = {C}; hit = None
                while q and hit is None:
                    x, route = q.popleft()
                    if len(route) >= MAXR: continue
                    for r in inputs_of[x]:
                        for y in outputs_of[r]:
                            if y == P: hit = route + [(r, y)]; break
                            if y not in seen and lp <= leaves(y):
                                seen.add(y); q.append((y, route + [(r, y)]))
                        if hit: break
                # a modification cycle: at least two reactions (one is a plain release), and
                # on the way P's protein is present only in a modified form (a form that
                # carries P's gene product but not P's own entity) - a released partner
                # that rides along unchanged (a carrier, specs/077) does not qualify
                pe = n["diagram_entity_id"]
                lv = lambda u: set((nodes.get(u, {}).get("member_leaves") or "").split("|")) - {""}
                modified = hit and any(pe not in lv(y) and pe != nodes.get(y, {}).get("diagram_entity_id")
                                       for _, y in hit[:-1])
                if hit and len(hit) >= 2 and modified:
                    found.append((P, B, C, hit))
    e = lambda u: nodes.get(u, {}).get("diagram_entity_id", "?")
    uniq = {}
    for P, B, C, hit in found:
        uniq.setdefault((e(P), e(C)), (P, B, C, hit, []))[4].append((P, B, C))
    if uniq:
        per[pw] = uniq; total += len(uniq)
        for (eP, eC), (P, B, C, hit, copies) in sorted(uniq.items()):
            print(pw, "P", eP, "C", eC, f"({len(copies)} node triples)", "route", " -> ".join(f"{e(r)}:{e(y)}" for r, y in hit))
print("entity-level candidates", total, "pathways", len(per))
