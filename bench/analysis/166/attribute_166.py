#!/usr/bin/env python3
"""specs/166 follow-up (Adam, 2026-11-20): attribute change calls in the `with` arm to the literature edges on
their routes, on the FIT axes only (korkut, own markers excluded; MP-BioPath experimental). Curator held-out (SELECT)
is never read.

A literature edge (u -> v) is ON a case's route when u is reachable from a node carrying the perturbed gene(s) and a
readout node is reachable from v, in the `with` network (directed, every edge type). This over-approximates the
bench's pins (any form of the gene, not only pinned roots). Each change call (predicted down/up) in `with` is
classed:
  new    `without` called no change on the same case (the literature made the call)
  kept   `without` already called a change
and right / wrong against the measured label. Per edge_key: counts per class, axis.

Usage: attribute_166.py WITH_BUILD WITH_RUN WITHOUT_RUN OUT_TSV     (needs Neo4j for korkut gene/site resolution)
"""
import csv
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "161"))

CODE = {"0": "down", "1": "not", "2": "up"}


def reach(adj, starts):
    seen, q = set(starts), deque(starts)
    while q:
        for w in adj.get(q.popleft(), ()):
            if w not in seen:
                seen.add(w)
                q.append(w)
    return seen


def on_route(lit, adj, radj, sources, readouts):
    """The literature edge keys (u, v, key) with u reachable from sources and a readout reachable from v."""
    fwd, back = reach(adj, sources), reach(radj, readouts)
    return {k for (u, v, k) in lit if u in fwd and v in back}


def classify(with_call, without_call, truth):
    if with_call == "not":
        return None
    return ("new" if without_call == "not" else "kept") + ("_right" if with_call == truth else "_wrong")


def load_pathway(pdir):
    adj, radj = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(pdir / "logic_network.csv")):
        adj[r["source_id"]].append(r["target_id"])
        radj[r["target_id"]].append(r["source_id"])
    lit = [(r["source_id"], r["target_id"], r["edge_key"])
           for r in csv.DictReader(open(pdir / "literature_origin.csv"))]
    return adj, radj, lit


def main():
    import dataset_roles
    import benchmark_vs_mpbiopath as bm
    import korkut_bench as kb
    from engagement_161 import marker_names
    dataset_roles.require_role("korkut2015", "fit")
    dataset_roles.require_role("mpbiopath_experimental", "fit")
    build, wrun, worun, out = (Path(a) for a in sys.argv[1:5])
    name_to_id = {}
    for r in csv.reader(open(HERE.parents[2] / "catalog_pathways.tsv"), delimiter="\t"):
        if len(r) >= 2 and r[0].startswith("R-HSA"):
            name_to_id[r[1]] = r[0]
    pdirs = {p.name: p for p in build.iterdir() if (p / "literature_origin.csv").exists()}
    nets = {pid: load_pathway(p) for pid, p in pdirs.items()}
    counts = defaultdict(Counter)        # edge_key -> Counter of "<axis>:<class>"
    cases = Counter()

    # experimental: the bench's own gene and readout node ids
    wo = {}
    for r in csv.DictReader(open(worun / "experimental_cases.tsv"), delimiter="\t"):
        if r.get("valid") == "1":
            wo[(r["pathway"], r["gene"], r["direction"], r["key_output"])] = CODE[r["predicted"]]
    for r in csv.DictReader(open(wrun / "experimental_cases.tsv"), delimiter="\t"):
        k = (r["pathway"], r["gene"], r["direction"], r["key_output"])
        pid = name_to_id.get(r["pathway"])
        if r.get("valid") != "1" or k not in wo or pid not in nets:
            continue
        cls = classify(CODE[r["predicted"]], wo[k], CODE[r["expected"]])
        if cls is None:
            continue
        cases["experimental:" + cls] += 1
        adj, radj, lit = nets[pid]
        srcs = [u for u in r["gene_uuids"].split("|") if u]
        outs = [u for u in r["output_uuids"].split("|") if u]
        for key in on_route(lit, adj, radj, srcs, outs):
            counts[key]["experimental:" + cls] += 1

    # korkut: targets and readouts resolved as the bench does, own markers excluded
    marks = marker_names()
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    genes = sorted({g for gs in drugs.values() for g in gs} |
                   {s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(genes)
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = kb.readout_stids(readouts, kb.total_gene_stids(readouts, g2s), kb.neo4j_site_stids(sites))
    totals = {r["readout"] for r in readouts if r["kind"] == "total"}
    wo = {(r["pathway_id"], r["drug_set"], r["readout"]): r["predicted"]
          for r in csv.DictReader(open(worun / "korkut_cases.tsv"), delimiter="\t") if r["status"] == "scored"}
    by_pw = defaultdict(list)
    for r in csv.DictReader(open(wrun / "korkut_cases.tsv"), delimiter="\t"):
        k = (r["pathway_id"], r["drug_set"], r["readout"])
        if r["status"] != "scored" or k not in wo or r["pathway_id"] not in nets:
            continue
        if any(r["readout"] in marks.get(d, set()) for d in r["drug_set"].split("+")):
            continue
        cls = classify(r["predicted"], wo[k], r["measured"])
        if cls is not None:
            by_pw[r["pathway_id"]].append((r, cls))
    for pid, rows in by_pw.items():
        s2u = bm.load_stid_to_uuids(pdirs[pid])
        ro_nodes = kb.total_readout_nodes(kb.readout_nodes(r_stids, s2u, None), totals, pdirs[pid])
        adj, radj, lit = nets[pid]
        for r, cls in rows:
            cases["korkut:" + cls] += 1
            gs = {g for c in r["drug_set"].split("+") for g in drugs[c]}
            srcs = [u for g in gs for s in g2s.get(g, []) for u in s2u.get(s, [])]
            for key in on_route(lit, adj, radj, srcs, list(ro_nodes.get(r["readout"], []))):
                counts[key]["korkut:" + cls] += 1

    cols = [f"{ax}:{c}" for ax in ("korkut", "experimental")
            for c in ("new_wrong", "new_right", "kept_wrong", "kept_right")]
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["edge_key"] + cols)
        for key in sorted(counts, key=lambda k: (-counts[k]["korkut:new_wrong"] - counts[k]["experimental:new_wrong"], k)):
            w.writerow([key] + [counts[key][c] for c in cols])
    print("change calls in `with` (FIT axes):", dict(cases))
    print(f"{len(counts)} literature edge keys on at least one route -> {out}")


if __name__ == "__main__":
    main()
