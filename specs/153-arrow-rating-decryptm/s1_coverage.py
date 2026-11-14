#!/usr/bin/env python3
"""specs/153 S1: structural coverage screen (uses no measurements).

For each served arrow (centre, source pathway, entity): which decryptM drugs
have a kinobeads target present in the source pathway's network with a
directed route to the arrow entity, and how many decryptM sites resolve in
the centre. Pure reachability on the sound Release-98 builds.

Usage: s1_coverage.py LEADS_JSON CASES_TSV BUILD92 BUILD33 OUT_TSV
"""
import csv, json, sys
from collections import defaultdict, deque
from pathlib import Path

csv.field_size_limit(10**7)
leads_p, cases_p, b92, b33, out_p = sys.argv[1:6]

leads = json.load(open(leads_p))
arrows = []
for cid, c in leads["centres"].items():
    for a in c["arrows"]:
        arrows.append({"centre": cid, "centre_name": c["name"], "source": a["source_id"],
                       "source_name": a["source_pathway"], "entity": a["entity_id"],
                       "entity_name": a["entity"], "rate_class": a["rate"]["class"]})

# decryptM cases per centre pathway
cases = list(csv.DictReader(open(cases_p), delimiter="\t"))
by_centre = defaultdict(list)
for r in cases:
    by_centre[r["pathway_id"]].append(r)

def pathway_dir(sid):
    for b in (b92, b33):
        d = Path(b) / sid
        if (d / "nodes.csv").exists():
            return d
    return None

# gene -> EWAS stIds via Neo4j 98 (one query per gene set)
sys.path.insert(0, "bench/analysis")
from curator_oracle import _cypher
genes = sorted({g for rows in by_centre.values() for r in rows for g in r["targets"].split(";")})
gene_stids = defaultdict(set)
for i in range(0, len(genes), 40):
    batch = genes[i:i+40]
    q = ("MATCH (e:EntityWithAccessionedSequence)-[:referenceEntity]->(r:ReferenceGeneProduct) "
         "WHERE e.speciesName='Homo sapiens' AND any(g IN r.geneName WHERE g IN [%s]) "
         "RETURN head([g IN r.geneName WHERE g IN [%s]]) + '|' + e.stId"
         % (",".join(f"'{g}'" for g in batch), ",".join(f"'{g}'" for g in batch)))
    for row in _cypher(q):
        g, s = row.split("|")
        gene_stids[g].add(s)

# an arrow entity may be a set/complex: anchor on its leaf descendants too
_entity_leaves = {}
def entity_leaves(stid):
    if stid not in _entity_leaves:
        rows = _cypher(
            "MATCH (e {stId:'%s'})-[:hasComponent|hasMember|hasCandidate*0..10]->(x) "
            "WHERE NOT (x)-[:hasComponent|hasMember|hasCandidate]->() RETURN DISTINCT x.stId" % stid)
        _entity_leaves[stid] = set(rows) | {stid}
    return _entity_leaves[stid]

def load_net(d):
    nodes, by_ent, leaves_of = {}, defaultdict(list), {}
    for r in csv.DictReader(open(d / "nodes.csv")):
        nodes[r["uuid"]] = r
        by_ent[r["diagram_entity_id"]].append(r["uuid"])
        leaves_of[r["uuid"]] = set(r["member_leaves"].split("|")) if r["member_leaves"] else set()
    out = defaultdict(list)
    with open(d / "logic_network.csv") as fh:
        rd = csv.DictReader(fh)
        for r in rd:
            out[r["source_id"]].append(r["target_id"])
    return nodes, by_ent, leaves_of, out

def reaches(starts, targets, out_edges):
    seen, dq = set(starts), deque(starts)
    tset = set(targets)
    while dq:
        u = dq.popleft()
        if u in tset:
            return True
        for v in out_edges.get(u, []):
            if v not in seen:
                seen.add(v); dq.append(v)
    return False

w = csv.writer(open(out_p, "w"), delimiter="\t")
w.writerow(["centre", "centre_name", "source", "source_name", "entity_name", "rate_class",
            "drug", "targets_in_source", "target_reaches_entity", "centre_sites", "centre_sites_regulated"])
nets = {}
summary = []
for a in arrows:
    d = pathway_dir(a["source"])
    rows = by_centre.get(a["centre"], [])
    site_n = len({(r["uniprot"], r["position"]) for r in rows})
    site_reg = len({(r["uniprot"], r["position"]) for r in rows if r["measured"] != "not"})
    if d is None or not rows:
        w.writerow([a["centre"], a["centre_name"], a["source"], a["source_name"], a["entity_name"],
                    a["rate_class"], "-", "-", "no_source_network" if d is None else "no_centre_sites", site_n, site_reg])
        summary.append((a, 0))
        continue
    if d not in nets:
        nets[d] = load_net(d)
    nodes, by_ent, leaves_of, out_edges = nets[d]
    eset = entity_leaves(a["entity"])
    ent_nodes = list(by_ent.get(a["entity"], [])) or \
        [u for u, ls in leaves_of.items() if a["entity"] in ls] or \
        [u for u, ls in leaves_of.items() if ls and ls <= eset | {a["entity"]} and ls & eset]
    drugs = defaultdict(set)
    for r in rows:
        for g in r["targets"].split(";"):
            drugs[r["drug"]].add(g)
    hit = 0
    for drug, tg in sorted(drugs.items()):
        present = sorted({g for g in tg if any(gene_stids.get(g, set()) & leaves_of[u] for u in nodes)})
        if not present or not ent_nodes:
            continue
        starts = [u for u in nodes if any(gene_stids.get(g, set()) & leaves_of[u] for g in present)]
        rch = reaches(starts, ent_nodes, out_edges)
        w.writerow([a["centre"], a["centre_name"], a["source"], a["source_name"], a["entity_name"],
                    a["rate_class"], drug, ";".join(present), rch, site_n, site_reg])
        if rch:
            hit += 1
    summary.append((a, hit))

print(f"{'arrow':60} {'rate':9} drugs_reaching")
for a, hit in summary:
    print(f"{a['source_name'][:30]:30} -> {a['entity_name'][:12]:12} -> {a['centre_name'][:14]:14} {a['rate_class']:9} {hit}")
