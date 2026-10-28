"""Hop distance from pinned targets to each korkut readout, per (network, drug set); joined to the case log."""
import collections, csv, sys
from pathlib import Path
sys.path.insert(0, "/home/awright/gitroot/deltasignal/bench")
import korkut_bench as kb
bm = kb.bm
build, cases, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
drugs = {d["code"]: d["genes"].split(",") for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}
readouts = [r for r in kb.read_tsv(kb.MAPS / "readouts.tsv") if r["kind"] != "exclude"]
drug_sets = sorted({ds for _, ds in kb.conditions()})
tg = sorted({g for gs in drugs.values() for g in gs})
rg = sorted({s.split(":")[0] for r in readouts for s in r["sites"].split(",")})
g2s = bm.neo4j_gene_to_stids(sorted(set(tg) | set(rg)))
sites = sorted({(s.split(":")[0], int(s.split(":")[1])) for r in readouts if r["kind"] == "phospho" for s in r["sites"].split(",")})
r_stids = kb.readout_stids(readouts, g2s, kb.neo4j_site_stids(sites))
need = collections.defaultdict(set)
for r in csv.DictReader(open(cases), delimiter="\t"):
    if r["status"] == "scored": need[r["pathway_id"]].add(r["drug_set"])
res = {}
for pid in sorted(need):
    pdir = build / pid
    s2u = bm.load_stid_to_uuids(pdir)
    ro = {n: list(dict.fromkeys(u for s in st for u in s2u.get(s, []))) for n, st in r_stids.items()}
    adj = bm.build_adjacency(pdir)
    for ds in drug_sets:
        key = "+".join(ds)
        if key not in need[pid]: continue
        genes = sorted({g for c in ds for g in drugs[c]})
        g2u = {g: [u for s in g2s.get(g, []) for u in s2u.get(s, [])] for g in genes}
        pins, _ = bm.pin_occurrences(pdir, g2u, g2s)
        pinned = set().union(*pins.values()) if pins else set()
        dist = {u: 0 for u in pinned}; fr = list(pinned); d = 0
        while fr:
            d += 1; nx = []
            for u in fr:
                for v in adj.get(u, ()):
                    if v not in dist: dist[v] = d; nx.append(v)
            fr = nx
        for n, nodes in ro.items():
            free = [u for u in nodes if u not in pinned and u in dist]
            res[(pid, key, n)] = min(dist[u] for u in free) if free else ""
    print(pid, flush=True)
with open(out, "w") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["pathway_id", "drug_set", "readout", "hops"])
    for k, v in sorted(res.items()): w.writerow([*k, v])
