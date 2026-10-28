#!/usr/bin/env python3
"""Hop distance from a korkut case's pinned targets to its readout nodes (graph only; specs/115, specs/116).

Resolves pins and readouts as `korkut_bench.py` does; shortest path in `build_adjacency` from any pinned node to any
unpinned readout node; blank if no route.

Usage: korkut_hops.py BUILD CASES_TSV OUT_TSV
"""
import collections
import csv
import sys
from pathlib import Path


def hops_from(adj, sources):
    """{node: shortest hop count from any source}, sources at 0."""
    dist = {u: 0 for u in sources}
    frontier, d = list(sources), 0
    while frontier:
        d += 1
        nxt = []
        for u in frontier:
            for v in adj.get(u, ()):
                if v not in dist:
                    dist[v] = d
                    nxt.append(v)
        frontier = nxt
    return dist


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import korkut_bench as kb
    bm = kb.bm
    build, cases, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
    drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
    readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
    drug_sets = sorted({ds for _, ds in kb.conditions()})
    tg = sorted({g for gs in drugs.values() for g in gs})
    rg = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
    g2s = bm.neo4j_gene_to_stids(sorted(set(tg) | set(rg)))
    sites = sorted({(s.split(":")[0], int(s.split(":")[1]))
                    for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
    r_stids = kb.readout_stids(readouts, g2s, kb.neo4j_site_stids(sites))
    need = collections.defaultdict(set)
    for r in csv.DictReader(open(cases), delimiter="\t"):
        if r["status"] == "scored":
            need[r["pathway_id"]].add(r["drug_set"])
    res = {}
    for pid in sorted(need):
        pdir = build / pid
        s2u = bm.load_stid_to_uuids(pdir)
        ro = {n: list(dict.fromkeys(u for s in st for u in s2u.get(s, []))) for n, st in r_stids.items()}
        adj = bm.build_adjacency(pdir)
        for ds in drug_sets:
            key = "+".join(ds)
            if key not in need[pid]:
                continue
            genes = sorted({g for c in ds for g in drugs[c]})
            g2u = {g: [u for s in g2s.get(g, []) for u in s2u.get(s, [])] for g in genes}
            pins, _ = bm.pin_occurrences(pdir, g2u, g2s)
            pinned = set().union(*pins.values()) if pins else set()
            dist = hops_from(adj, pinned)
            for n, nodes in ro.items():
                free = [u for u in nodes if u not in pinned and u in dist]
                res[(pid, key, n)] = min(dist[u] for u in free) if free else ""
        print(pid, flush=True)
    with open(out, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pathway_id", "drug_set", "readout", "hops"])
        for k, v in sorted(res.items()):
            w.writerow([*k, v])


if __name__ == "__main__":
    main()
