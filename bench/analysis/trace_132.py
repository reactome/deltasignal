#!/usr/bin/env python3
"""specs/132: for every case an arm changes (under F, so only cases the pathway alone cannot route), which curated
connections the arm's route crosses.

A bridge (a `context_bridge` edge) is on a case's route when the case's sources reach its tail and its head reaches the
case's readout, in the arm's bundle. Sources are the perturbed gene's nodes (experimental, curator: the bench's
`gene_uuids`; korkut: every node of the drug targets' entities in the bundle, a superset of the roots the bench pins,
so a korkut case may list a bridge its pins do not use). Each bridge is named by its census rows (s in Q precedes r in
P, through entity e in role): the rows of P whose entity is the head node's entity, or a set holding it.

Writes one TSV row per changed case, and prints the cases grouped by the connections they cross.
Usage: trace_132.py CONTROL_RESULTS ARM_RESULTS CONTROL_BUNDLE ARM_BUNDLE MEMBERS_TSV OUT.tsv
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import combine_132 as c  # noqa: E402
import korkut_bench as kb  # noqa: E402


def bridges(pdir):
    """[(tail, head)] of the bundle's context_bridge edges."""
    with bm.open_network(pdir) as f:
        return [(r["source_id"], r["target_id"]) for r in csv.DictReader(f) if r["edge_type"] == "context_bridge"]


def reverse(adj):
    rev = defaultdict(list)
    for u, vs in adj.items():
        for v in vs:
            rev[v].append(u)
    return rev


def on_route(br, fwd, rev, sources, readout):
    ahead, behind = bm.reachable_from(fwd, set(sources)), bm.reachable_from(rev, set(readout))
    return [b for b in br if b[0] in ahead and b[1] in behind]


def name_bridge(head, census, stids_of, members):
    """Census rows (s, r, q, entity, role) whose entity is the head's entity or a set holding it."""
    own = stids_of.get(head, set())
    return sorted({(row["s"], row["r"], row["q"], row["entity"], row["role"]) for row in census
                   if row["entity"] in own or members.get(row["entity"], set()) & own})


def names(stids):
    from py2neo import Graph
    graph = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    rows = graph.run("UNWIND $ids AS i MATCH (n:DatabaseObject {stId: i}) RETURN i, n.displayName AS name",
                     ids=sorted(stids)).data()
    return {r["i"]: r["name"] for r in rows}


def korkut_sources():
    """{drug set: target genes} and {gene: stIds}, {readout: stIds}, from the bench's own maps."""
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    genes = sorted({g for gs in drugs.values() for g in gs}
                   | {s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(genes)
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    return drugs, g2s, kb.readout_stids(readouts, g2s, kb.neo4j_site_stids(sites))


def main():
    c_res, a_res, c_bundle, a_bundle, members_tsv, out = sys.argv[1:]
    a_bundle, ids = Path(a_bundle), c.pathway_ids()
    control, arm = c.load_all(c_res, c_bundle), c.load_all(a_res, a_bundle)
    census = defaultdict(list)
    for row in csv.DictReader(open(a_bundle / "curated_census.tsv"), delimiter="\t"):
        if row["status"] == "bridged":
            census[row["pathway"]].append(row)
    members = defaultdict(set)
    for row in csv.DictReader(open(members_tsv), delimiter="\t"):
        members[row["entity"]].add(row["member"])
    cell_rows = {}      # (pathway, gene, direction, key_output) -> (sources, readout), from the arm's case file
    for kind in ("experimental", "curator"):
        for r in csv.DictReader(open(Path(a_res) / f"{kind}_cases.tsv", newline=""), delimiter="\t"):
            cell_rows[(r["pathway"], r["gene"], r["direction"], r["key_output"])] = (
                [u for u in r["gene_uuids"].split("|") if u], [u for u in r["output_uuids"].split("|") if u])
    drugs = g2s = r_stids = None
    cache, rows = {}, []
    for d in ("experimental", "curator_heldout", "korkut"):
        moved = sorted((k for k in control[d].keys() & arm[d].keys()
                        if not control[d][k]["route"] and arm[d][k]["predicted"] != control[d][k]["predicted"]), key=str)
        for k in moved:
            x, y = control[d][k], arm[d][k]
            pid = ids[k[0]] if d != "korkut" else k[0]
            if pid not in cache:
                pdir = a_bundle / pid
                fwd = bm.build_adjacency(pdir)
                s2u = bm.load_stid_to_uuids(pdir)
                stids_of = defaultdict(set)
                for s, us in s2u.items():
                    for u in us:
                        stids_of[u].add(s)
                cache[pid] = (fwd, reverse(fwd), bridges(pdir), s2u, stids_of, bm.context_own_uuids(pdir))
            fwd, rev, br, s2u, stids_of, own = cache[pid]
            if d == "korkut":
                if drugs is None:
                    drugs, g2s, r_stids = korkut_sources()
                genes = {g for code in k[1].split("+") for g in drugs[code]}
                sources = [u for g in genes for s in g2s.get(g, []) for u in s2u.get(s, [])]
                readout = kb.readout_nodes({k[2]: r_stids.get(k[2], [])}, s2u, own).get(k[2], [])
                case = f"{k[1]} -> {k[2]}"
            else:
                sources, readout = cell_rows[k]
                case = f"{k[1]} {'KO' if k[2] == '0' else 'OE' if k[2] == '2' else k[2]} -> {k[3]}"
            route = on_route(br, fwd, rev, sources, readout)
            conns = sorted({n for b in route for n in name_bridge(b[1], census[pid], stids_of, members)})
            verdict = ("fixed" if y["predicted"] == y["measured"] else
                       "broken" if x["predicted"] == x["measured"] else "wrong_to_wrong")
            rows.append({"dataset": d, "pathway": pid, "case": case, "control": x["predicted"],
                         "arm": y["predicted"], "measured": y["measured"], "verdict": verdict,
                         "bridges": len(route), "connections": conns})
    label = names({v for r in rows for t in r["connections"] for v in t[:4]} | {r["pathway"] for r in rows})
    fmt = lambda t: f"{label.get(t[2], t[2])}: {label.get(t[0], t[0])} -> {label.get(t[1], t[1])} [{label.get(t[3], t[3])}, {t[4]}]"
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["dataset", "pathway", "pathway_name", "case", "control", "arm", "measured", "verdict", "bridges",
                    "connections"])
        for r in rows:
            w.writerow([r["dataset"], r["pathway"], label.get(r["pathway"], r["pathway"]), r["case"], r["control"],
                        r["arm"], r["measured"], r["verdict"], r["bridges"], " | ".join(fmt(t) for t in r["connections"])])
    for d in ("experimental", "curator_heldout", "korkut"):
        by = defaultdict(Counter)
        for r in rows:
            if r["dataset"] == d:
                by[(label.get(r["pathway"], r["pathway"]), " | ".join(fmt(t) for t in r["connections"]) or "(none found)")][r["verdict"]] += 1
        print(f"== {d}")
        for (p, conn), v in sorted(by.items(), key=lambda kv: -sum(kv[1].values())):
            print(f"  [{p}] fixed {v['fixed']} broken {v['broken']} wrong->wrong {v['wrong_to_wrong']}\n      {conn}")


if __name__ == "__main__":
    main()
