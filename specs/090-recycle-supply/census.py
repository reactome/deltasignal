"""specs/090 census: pool states whose supply mixes several producers.

The solver's pool supply (cycles.jl update_pool!) is the plain mean, over every
state, of each producer's fold, where a state's producers are its positive
in-edges other than the pool's own step reaction copies. A state assembled from
a protein and its partners therefore averages them. Per pool kind (pool,
recycle, release): pools, states, states with more than one producer, and the
edge types and AND/OR of those producers.
Usage: python3 census.py <build dir>
"""
import csv, os, sys
from collections import Counter, defaultdict

build = sys.argv[1]
tot = Counter(); types = Counter(); ex = []
for pw in sorted(p for p in os.listdir(build) if p.startswith("R-HSA-")):
    d = os.path.join(build, pw)
    if not os.path.exists(os.path.join(d, "pools.csv")):
        continue
    forms = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(d, "pools.csv"))):
        forms[r["pool_id"]].append(r)
    steps = {r["reaction_uuid"] for r in csv.DictReader(open(os.path.join(d, "pool_transitions.csv")))}
    into = defaultdict(list)
    for f in ("logic_network.csv", "boundary_edges.csv"):
        p = os.path.join(d, f)
        if not os.path.exists(p):
            continue
        for e in csv.DictReader(open(p)):
            if e["pos_neg"] == "pos" and e["edge_type"] not in ("recycle",) and e["source_id"] not in steps:
                into[e["target_id"]].append((e["edge_type"], e["and_or"]))
    for pid, rows in forms.items():
        kind = pid.rstrip("0123456789")
        tot[(kind, "pools")] += 1
        multi = 0
        for r in rows:
            if r["role"] != "state":
                continue
            tot[(kind, "states")] += 1
            prods = into.get(r["node_uuid"], [])
            if len(prods) > 1:
                multi += 1
                tot[(kind, "states >1 producer")] += 1
                for t in prods:
                    types[(kind,) + t] += 1
        if multi:
            tot[(kind, "pools with such a state")] += 1
            ex.append((pw, pid, multi))
for k in sorted(tot):
    print(*k, tot[k], sep="\t")
print("producer edges of those states (kind, edge_type, and_or):")
for k, n in sorted(types.items()):
    print("  ", *k, n, sep="\t")
print("pools:", ex)
