"""Curation-review list: held-out no_path perturbations whose route ends in
Reactome itself.

For each (pathway, gene) with no_path cases in a case log, walk forward from
the gene's pinned nodes over every edge except dissociation, and collect the
"sinks": reached entity nodes with no out-edge except dissociation. Each sink's
diagram entity is checked against the pathway in Neo4j (any reaction in the
pathway's event hierarchy using it, or a complex/set containing it, as input,
catalyst or regulator). A perturbation whose sinks are all unused in the
pathway is a curation dead end: the network is faithful, Reactome has no
route on.

    python3 dead_ends.py CASES.tsv BUILD_DIR CATALOG_PATHWAYS.tsv OUT.tsv

Needs cypher-shell in the reactome-neo4j container (password read at run
time inside the container, never printed).
"""
import collections
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from holdout_report import TUNING_PATHWAYS  # noqa: E402


def cypher(q):
    r = subprocess.run(["docker", "exec", "-i", "reactome-neo4j", "sh", "-c",
                        'cypher-shell -u neo4j -p "${NEO4J_AUTH#neo4j/}" --format plain'],
                       input=q, capture_output=True, text=True, check=True)
    return r.stdout.splitlines()[1:]


def main(cases, build, catalog, out):
    pid = {}
    for r in csv.reader(open(catalog), delimiter="\t"):
        if r and r[0].startswith("R-HSA"):
            pid[r[1]] = r[0]
    rows = [r for r in csv.DictReader(open(cases), delimiter="\t")
            if r["valid"] == "1" and not r.get("exclusion") and r["pathway"] not in TUNING_PATHWAYS]
    nopath = collections.defaultdict(list)
    for r in rows:
        if r["category"] == "no_path":
            nopath[(r["pathway"], r["gene"])].append(r)
    results = []
    for (pw, gene), rs in sorted(nopath.items()):
        d = os.path.join(build, pid[pw])
        st = {a: b for a, b in csv.reader(open(os.path.join(d, "stid_to_uuid_mapping.csv")))}
        kind, ident = {}, {}
        for n in csv.DictReader(open(os.path.join(d, "nodes.csv"))):
            kind[n["uuid"]], ident[n["uuid"]] = n["node_kind"], n["diagram_entity_id"]
        outs = collections.defaultdict(list)
        for f in ("logic_network.csv", "boundary_edges.csv"):
            p = os.path.join(d, f)
            if os.path.exists(p):
                for e in csv.DictReader(open(p)):
                    outs[e["source_id"]].append((e["target_id"], e["edge_type"]))
        pins = [u for u in rs[0]["gene_uuids"].split("|") if u]
        seen, stack = set(pins), list(pins)
        while stack:
            u = stack.pop()
            for v, et in outs[u]:
                if et != "dissociation" and v not in seen:
                    seen.add(v)
                    stack.append(v)
        sinks = sorted({ident.get(u) or st.get(u, "").split("::")[0] for u in seen
                        if kind.get(u) not in ("reaction", "dissociation")
                        and not any(et != "dissociation" for _, et in outs[u])})
        sinks = [s for s in sinks if s]
        results.append((pw, gene, len(rs), len(pins), len(seen), sinks))
    allsinks = sorted({(pid[pw], s) for pw, _, _, _, _, ss in results for s in ss})
    used = set()
    for i in range(0, len(allsinks), 200):
        chunk = [list(x) for x in allsinks[i:i + 200]]
        q = ("UNWIND %s AS x MATCH (p:Pathway {stId:x[0]}), (e:PhysicalEntity {stId:x[1]}) "
             "WHERE EXISTS { MATCH (p)-[:hasEvent*]->(r:ReactionLikeEvent), "
             "(c:PhysicalEntity)-[:hasComponent|hasMember|hasCandidate*0..4]->(e) "
             "WHERE (r)-[:input]->(c) OR (r)-[:catalystActivity]->()-[:physicalEntity]->(c) "
             "OR (r)-[:regulatedBy]->()-[:regulator]->(c) } "
             "RETURN x[0], x[1];") % json.dumps(chunk)
        for line in cypher(q):
            a, b = json.loads("[" + line + "]")
            used.add((a, b))
    names = {}
    ids = sorted({s for _, s in allsinks})
    for line in cypher("UNWIND %s AS s MATCH (e {stId:s}) RETURN s, e.displayName;" % json.dumps(ids)):
        a, b = json.loads("[" + line + "]")
        names[a] = b
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pathway", "gene", "no_path_cases", "pins", "reached", "class",
                    "sinks_unused_in_pathway", "sinks_used_in_pathway"])
        tally = collections.Counter()
        for pw, gene, n, npins, nreach, sinks in results:
            unused = [s for s in sinks if (pid[pw], s) not in used]
            usedl = [s for s in sinks if (pid[pw], s) in used]
            cls = ("no_pin" if not npins else "dead_end" if sinks and not usedl
                   else "no_sink" if not sinks else "severed_in_network")
            tally[cls] += n
            fmt = lambda xs: "; ".join(f"{s} {names.get(s, '')}" for s in xs)
            w.writerow([pw, gene, n, npins, nreach, cls, fmt(unused), fmt(usedl)])
    print(dict(tally))


if __name__ == "__main__":
    main(*sys.argv[1:5])
