# specs/080: for the readout-form breaks, the gene-containing nodes that DO reach the
# benchmark's readout copies, classified: a root (no in-edge: the pin choice missed it),
# regenerated (every in-edge comes from its own downstream: a closed cycle, effectively
# a starting form), or produced (made from upstream the pinned copies do not reach).
# Usage: python3 classify_readout_form.py breaks_<sha>.tsv <build dir> <results dir>
import csv, os, sys
from collections import Counter, defaultdict
brk, B, R = sys.argv[1:4]
name = {}
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../bench/catalog_pathways.tsv")):
    if line.startswith(("#", "id\t")) or not line.strip(): continue
    i, n = line.rstrip("\n").split("\t")[:2]; name[i] = n
cases = defaultdict(list)
for r in csv.DictReader(open(f"{R}/curator_cases.tsv"), delimiter="\t"):
    if r["valid"] == "1": cases[(r["pathway"], r["gene"], r["key_output"])].append(r)
want = [r for r in csv.DictReader(open(brk), delimiter="\t") if r["break_kind"] == "readout form" and int(r["wrong"]) > 0]
tally = Counter(); wrong = Counter(); ex = defaultdict(list)
cache = {}
for w in want:
    pw = w["pathway_id"]; d = f"{B}/{pw}"
    if pw not in cache:
        sid = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
        cont = defaultdict(set)
        for x in csv.DictReader(open(d + "/containment.csv")): cont[x["stable_id"]].add(x["contains_stable_id"])
        adj = defaultdict(set); inc = defaultdict(set)
        for f in ("logic_network.csv", "boundary_edges.csv"):
            if os.path.exists(f"{d}/{f}"):
                for e in csv.DictReader(open(f"{d}/{f}")):
                    adj[e["source_id"]].add(e["target_id"]); inc[e["target_id"]].add(e["source_id"])
        cache[pw] = (sid, cont, adj, inc)
    sid, cont, adj, inc = cache[pw]
    def reach(starts):
        seen = set(starts); fr = list(starts)
        while fr:
            nx = []
            for u in fr:
                for v in adj[u]:
                    if v not in seen: seen.add(v); nx.append(v)
            fr = nx
        return seen
    rs = [r for r in cases[(name[pw], w["gene"], w["key_output"])]]
    if not rs: continue
    r0 = rs[0]
    g = r0["gene_uuids"].split("|"); o = set(r0["output_uuids"].split("|"))
    gst = {sid[u].split("::")[0] for u in g if u in sid}
    holders = [u for u, s in sid.items() if any(x in cont.get(s, ()) or x == s.split("::")[0] for x in gst)]
    hit = [u for u in holders if o & reach([u])]
    kinds = set()
    for u in hit:
        if not inc[u]: kinds.add("root")
        elif all(u in reach([a]) for a in inc[u]): kinds.add("regenerated")
        else: kinds.add("produced")
    k = "root" if "root" in kinds else ("regenerated" if "regenerated" in kinds else ("produced" if kinds else "none reaches"))
    n_wrong = sum(x["predicted"] != x["expected"] for x in rs)
    tally[k] += len(rs); wrong[k] += n_wrong
    ex[k].append((name[pw][:28], w["gene"], w["key_output"], n_wrong))
for k in tally:
    print(f"{k:14s} cases {tally[k]:4d} wrong {wrong[k]:4d}  e.g.", ex[k][:5])
