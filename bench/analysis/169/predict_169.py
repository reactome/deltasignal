#!/usr/bin/env python3
"""specs/169: the served model's predictions for split-A CRISPRi knockdowns at curated transcription targets.

For each pathway of BUILD and each split-A perturbed gene G of the K562 genome-wide screen whose forms the protocol
of record can pin there (/api/pins basis root or root_cycle), the solve pins G's roots at 0 (the knockout protocol,
DS_PERTURB_UI_DOWN=0) and /api/readouts values every measured target gene T != G: the output entity of a curated
gene-expression-regulated event whose product is T, present in the pathway (forms mode, as the curator bench reads a
curated entity). Writes one row per (pathway, G, T): fold, call (served cutoffs), converged.

Usage: predict_169.py BUILD_DIR PORT H5AD OUT_TSV      (needs Neo4j; the server on PORT serves BUILD)
"""
import csv
import json
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent))
from common_169 import split_of, gene_of, h5_col, obs_index  # noqa: E402

TARGETS = """MATCH (r:ReactionLikeEvent)-[:regulatedBy]->(reg)
WHERE r.speciesName = 'Homo sapiens' AND (reg:PositiveGeneExpressionRegulation OR reg:NegativeGeneExpressionRegulation)
MATCH (r)-[:output]->(o)-[:referenceEntity]->(tre) WHERE tre.geneName IS NOT NULL
RETURN DISTINCT o.stId AS s, tre.geneName[0] AS g"""


def post(port, ep, body):
    req = urllib.request.Request(f"http://localhost:{port}/api/{ep}", data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=3600) as r:
        return json.load(r)


GENE_STIDS = """MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
WHERE re.geneName IS NOT NULL RETURN re.geneName AS names, pe.stId AS s"""


def gene_stids(graph, genes):
    """{gene: [stIds]} as bm.neo4j_gene_to_stids (a gene in a reference entity's geneName list), in one pass: the
    UNWIND form scans every reference entity once per gene, about half an hour for 5,000 genes."""
    out = defaultdict(set)
    for r in graph.run(GENE_STIDS):
        for name in r["names"] or ():
            if name in genes:
                out[name].add(r["s"])
    return {g: sorted(v) for g, v in out.items()}


def screen_genes(h5):
    import h5py
    with h5py.File(h5, "r") as f:
        ids, ctrl = obs_index(f), h5_col(f, "obs", "core_control")
        measured = set(h5_col(f, "var", "gene_name"))
    perturbed = {gene_of(i) for i, c in zip(ids, ctrl) if not c and gene_of(i)}
    return perturbed, measured


def main():
    import dataset_roles
    import benchmark_vs_mpbiopath as bm
    from py2neo import Graph
    dataset_roles.require_role("replogle2022_a", "score")
    build, port, h5, out = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
    perturbed, measured = screen_genes(h5)
    genes_a = sorted(g for g in perturbed if split_of(g) == "A")
    g = Graph(bm.NEO4J_URL, auth=(bm.NEO4J_USER, bm.NEO4J_PASSWORD))
    target_stids = defaultdict(set)
    for r in g.run(TARGETS):
        if r["g"] in measured:
            target_stids[r["g"]].add(r["s"])
    g2s = gene_stids(g, set(genes_a))
    rows, n_solves = [], 0
    for pdir in sorted(p for p in build.iterdir() if (p / "logic_network.csv").exists()):
        s2u = bm.load_stid_to_uuids(pdir)
        here_t = {t: sorted(s for s in ss if s2u.get(s)) for t, ss in target_stids.items()}
        here_t = {t: ss for t, ss in here_t.items() if ss}
        leaves = set()
        for r in csv.DictReader(open(pdir / "nodes.csv")):
            leaves.add(r["diagram_entity_id"])
            leaves.update(r["member_leaves"].split("|"))
        cand = [gn for gn in genes_a if set(g2s.get(gn, [])) & leaves]
        if not here_t or not cand:
            continue
        nid = post(port, "parse", {"pathway_id": pdir.name})["network_id"]
        pins = post(port, "pins", {"network_id": nid, "queries": [{"id": gn, "stable_ids": g2s[gn]} for gn in cand]})
        for res in pins["results"]:
            gn = res["id"]
            if res["basis"] not in ("root", "root_cycle") or not res["pins"]:
                continue
            pinned = [p["uuid"] for p in res["pins"]]
            sol = post(port, "solve", {"network_id": nid, "observations": {u: [0.0, 1.0] for u in pinned}})
            n_solves += 1
            qs = [{"id": t, "stable_ids": ss} for t, ss in sorted(here_t.items()) if t != gn]
            if not qs:                       # the pathway's only target is the perturbed gene itself
                continue
            ro = post(port, "readouts", {"network_id": nid, "activities": sol.get("node_activities") or {},
                                          "pinned": pinned, "queries": qs, "converged": bool(sol.get("converged"))})
            for e in ro["results"]:
                if e.get("fold") is None:
                    continue
                rows.append({"pathway": pdir.name, "gene": gn, "target": e["id"], "fold": repr(e["fold"]),
                             "call": {"normal": "not"}.get(e["call"], e["call"]),
                             "converged": int(bool(sol.get("converged")))})
        print(f"  {pdir.name}: {len(cand)} candidate genes", flush=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pathway", "gene", "target", "fold", "call", "converged"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{n_solves} solves, {len(rows)} (pathway, gene, target) predictions -> {out}")


if __name__ == "__main__":
    main()
