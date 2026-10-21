# specs/080 breaks.py: for every case the database reaches (v5) and our network does not,
# where the shortest database route first leaves what our network reaches, and by
# which kind of step (reaction edge, member -> set, component -> root complex).
# Usage: python3 breaks.py reach_<sha>.tsv <build dir> <mp-biopath dir> <out.tsv>
# (shared reachability rules below)
# specs/080: for every benchmark case, is the readout reachable from the gene in the
# DATABASE (Neo4j, the pathway's own reactions), in OUR network, and (10 pathways) in
# MP-BioPath's published distance table?
#
# Database reachability (v1), directed and scoped to one pathway's ReactionLikeEvents
# (hasEvent*):
#   start   every participant of those reactions that contains the gene's product
#           (hasComponent|hasMember|hasCandidate, any depth), as the benchmark resolves
#           a gene to every network copy that is or contains it;
#   steps   entity -> reaction when it is an input, a catalyst or a regulator (either
#           sign); reaction -> each of its outputs;
#   reached the readout when a reached participant is, or contains, the readout entity
#           (the benchmark's readout forms).
#   v2 = v1 + set membership both ways (a set is "any of" its members);
#   v3 = v2 + a component reaches a complex that contains it when no reaction in the
#        pathway produces that complex (a root complex, as the generator decomposes);
#   v4 = v2 + a component reaches every complex that contains it (loose upper bound);
#   v5 = v1 + a member reaches its set (one way: a member's change is not its siblings')
#        + the v3 root-complex step. The strict reading.
# Our network: forward over every edge of logic_network.csv + boundary_edges.csv, from
# the case's gene uuids to its readout uuids (the benchmark's own no_path test, run on
# every case, not only wrong ones).
#
# Usage (needs py2neo; NEO4J_PASSWORD in the environment):
#   python3 db_reach.py <results dir> <catalog build dir> <mp-biopath-pathways dir> <out.tsv>
import csv, os, sys
from collections import defaultdict
from py2neo import Graph

REACH, BUILD, MPB, OUT = sys.argv[1:5]; R = None
g = Graph(os.environ.get("NEO4J_URL", "bolt://localhost:7687"),
          auth=("neo4j", os.environ["NEO4J_PASSWORD"]))
name2id = {}
for line in open(os.path.join(os.path.dirname(__file__), "../../bench/catalog_pathways.tsv")):
    if line.startswith(("#", "id\t")) or not line.strip(): continue
    i, n = line.rstrip("\n").split("\t")[:2]; name2id[n] = i

mpb = set()
p = f"{MPB}/analysis_results/input_to_keyoutput_distance.tsv"
for r in csv.DictReader(open(p), delimiter="\t"):
    mpb.add((r["pathway_id"], r["input_name"], r["keyoutput_id"]))
mpb_pathways = {k[0] for k in mpb}

def db_graph(pw):
    rows = g.run("""
    MATCH (:Pathway {stId:$pw})-[:hasEvent*]->(r:ReactionLikeEvent)
    OPTIONAL MATCH (r)-[:input]->(i:PhysicalEntity)
    OPTIONAL MATCH (r)-[:output]->(o:PhysicalEntity)
    OPTIONAL MATCH (r)-[:catalystActivity]->(:CatalystActivity)-[:physicalEntity]->(c:PhysicalEntity)
    OPTIONAL MATCH (r)-[:regulatedBy]->(:Regulation)-[:regulator]->(g2:PhysicalEntity)
    RETURN r.stId AS r, collect(DISTINCT i.stId) AS i, collect(DISTINCT o.stId) AS o,
           collect(DISTINCT c.stId) AS c, collect(DISTINCT g2.stId) AS g""", pw=pw).data()
    fwd = defaultdict(set); parts = set(); produced = set()
    for x in rows:
        for e in x["i"] + x["c"] + x["g"]:
            fwd[e].add(x["r"]); parts.add(e)
        for e in x["o"]:
            fwd[x["r"]].add(e); parts.add(e); produced.add(e)
    # every containment edge below a participant, with its kind
    edges = g.run("""
    UNWIND $parts AS s MATCH (p:PhysicalEntity {stId:s})
    MATCH (p)-[:hasComponent|hasMember|hasCandidate*0..8]->(a:PhysicalEntity)-[rel:hasComponent|hasMember|hasCandidate]->(b:PhysicalEntity)
    RETURN DISTINCT a.stId AS a, type(rel) AS t, b.stId AS b""", parts=sorted(parts)).data()
    leaves = g.run("""
    UNWIND $parts AS s MATCH (p:PhysicalEntity {stId:s})
    MATCH (p)-[:hasComponent|hasMember|hasCandidate*0..8]->(x:PhysicalEntity)
    OPTIONAL MATCH (x)-[:referenceEntity]->(re:ReferenceEntity)
    RETURN DISTINCT p.stId AS p, x.stId AS x, x.dbId AS xdb, re.geneName AS genes""", parts=sorted(parts)).data()
    has_gene = defaultdict(set); gene_leaf = defaultdict(set); has_db = defaultdict(set)
    for x in leaves:
        has_db[str(x["xdb"])].update((x["p"], x["x"]))
        for gn in x["genes"] or []:
            has_gene[gn].add(x["p"]); gene_leaf[gn].add(x["x"])
    sets = defaultdict(set); set_up = defaultdict(set); up_root = defaultdict(set); up_all = defaultdict(set)
    for e in edges:
        if e["t"] in ("hasMember", "hasCandidate"):
            sets[e["a"]].add(e["b"]); sets[e["b"]].add(e["a"]); set_up[e["b"]].add(e["a"])
        else:
            up_all[e["b"]].add(e["a"])
            if e["a"] not in produced: up_root[e["b"]].add(e["a"])
    def merge(*ds):
        m = defaultdict(set)
        for d in ds:
            for k, v in d.items(): m[k] |= v
        return m
    global KIND
    KIND = {}
    for a, bs in fwd.items():
        for b in bs: KIND[(a, b)] = "reaction"
    for a, bs in set_up.items():
        for b in bs: KIND.setdefault((a, b), "member->set")
    for a, bs in up_root.items():
        for b in bs: KIND.setdefault((a, b), "component->root complex")
    variants = {"v1": fwd, "v2": merge(fwd, sets), "v3": merge(fwd, sets, up_root), "v4": merge(fwd, sets, up_all),
                "v5": merge(fwd, set_up, up_root)}
    return variants, has_gene, gene_leaf, has_db, len(rows)

def bfs(adj, starts):
    seen = set(starts); fr = list(starts)
    while fr:
        nx = []
        for u in fr:
            for v in adj.get(u, ()):
                if v not in seen: seen.add(v); nx.append(v)
        fr = nx
    return seen

def our_adj(pw):
    adj = defaultdict(set)
    for f in ("logic_network.csv", "boundary_edges.csv"):
        fp = f"{BUILD}/{pw}/{f}"
        if os.path.exists(fp):
            for e in csv.DictReader(open(fp)): adj[e["source_id"]].add(e["target_id"])
    return adj


def path(adj, starts, targets):
    prev = {s: None for s in starts}; fr = list(starts)
    while fr:
        nx = []
        for u in fr:
            if u in targets:
                out = [u]
                while prev[out[-1]] is not None: out.append(prev[out[-1]])
                return out[::-1]
            for v in adj.get(u, ()):
                if v not in prev: prev[v] = u; nx.append(v)
        fr = nx
    return None

cases = defaultdict(list)
for r in csv.DictReader(open(REACH), delimiter="\t"):
    if r["db_v5"] == "1" and r["our_reach"] == "0":
        cases[(name2id[r["pathway"]], r["gene"], r["key_output"])].append(r)
w = csv.writer(open(OUT, "w", newline=""), delimiter="\t")
w.writerow(["pathway_id", "gene", "key_output", "cases", "wrong", "break_kind", "last_reached", "first_missed", "missed_in_network"])
by_pw = defaultdict(list)
for k in cases: by_pw[k[0]].append(k)
for pw, ks in by_pw.items():
    var, hg, gl, hd, n = db_graph(pw)
    d = f"{BUILD}/{pw}"
    sid = {r["uuid"]: r["diagram_entity_id"] for r in csv.DictReader(open(d + "/nodes.csv"))}
    for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv")):
        sid.setdefault(r["uuid"], r["stable_id"].split("::")[0])
    present = set(sid.values()); adj = our_adj(pw)
    for k in ks:
        _, gene, ko = k
        starts = hg.get(gene, set()) | gl.get(gene, set())
        p = path(var["v5"], starts, hd.get(ko, set()))
        reach = bfs(adj, {u for u, s_ in sid.items() if s_ in starts}); rs = {sid.get(u) for u in reach}
        kind, last, miss = "start", "", p[0]
        # a node our network does not have (a set it expanded, a reaction it merged)
        # is transparent: the break is the first node it HAS and does not reach
        prev_ = p[0]
        for a, b in zip(p, p[1:]):
            if b in present and b not in rs:
                kind, last, miss = KIND.get((a, b), "?"), prev_, b
                break
            if b in rs: prev_ = b
        else:
            kind, miss = "readout form", p[-1]
        rr = cases[k]
        w.writerow([pw, gene, ko, len(rr), sum(x["category"] == "no_path" for x in rr), kind, last, miss, int(miss in present)])
    print(pw, len(ks), file=sys.stderr)
