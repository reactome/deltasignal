# specs/088: open carrier-loop releases that DS_CARRIER_LOOP_OPEN=fed_single keeps as
# "displacements" (the release reaction has >= 2 distinct non-cofactor input parents), split
# by whether the second input is only ubiquitin (a donor co-substrate, not a displacer).
# Usage: python3 count_ub_displacers.py <build dir>
import csv, os, sys
from collections import Counter, defaultdict
UB = {"R-HSA-68524", "R-HSA-113595", "R-HSA-9660007"}
b = sys.argv[1]
tot = Counter(); per = Counter(); ex = []
for pw in sorted(p for p in os.listdir(b) if p.startswith("R-HSA-")):
    d = f"{b}/{pw}"
    if not os.path.exists(d + "/carrier_loops.csv"): continue
    nodes = {r["uuid"]: r for r in csv.DictReader(open(d + "/nodes.csv"))}
    ent = lambda u: nodes.get(u, {}).get("diagram_entity_id", "")
    cof = set()
    if os.path.exists(d + "/cofactors.csv"):
        cof = {list(r.values())[0] for r in csv.DictReader(open(d + "/cofactors.csv"))}
    rows = list(csv.DictReader(open(d + "/carrier_loops.csv")))
    listed = defaultdict(set)
    for r in rows: listed[r["free_uuid"]].add(r["loop_source_uuid"])
    pos = defaultdict(set); ins = defaultdict(set)
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if not os.path.exists(f"{d}/{f}"): continue
        for e in csv.DictReader(open(f"{d}/{f}")):
            if e["pos_neg"] != "pos": continue
            pos[e["target_id"]].add(e["source_id"])
            if e["edge_type"] == "input" and ent(e["source_id"]) not in cof:
                ins[e["target_id"]].add(e["source_id"])
    for r in rows:
        f, s = r["free_uuid"], r["loop_source_uuid"]
        if r["closed"] != "False" or not (pos[f] - listed[f]): continue
        if len(ins[s]) < 2: continue                       # dropped by fed_single
        others = {ent(u) for u in ins[s]}
        kind = "ub only" if others & UB and len(others - UB) <= 1 else "other displacer"
        tot[kind] += 1; per[(pw, kind)] += 1
        if kind == "ub only": ex.append((pw, ent(f), ent(s)))
print(dict(tot))
for (pw, k), n in sorted(per.items(), key=lambda x: -x[1])[:12]: print("  ", pw, k, n)
for x in sorted(set(ex))[:15]: print("  ", x)
