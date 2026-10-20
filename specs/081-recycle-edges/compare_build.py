# specs/081: does a variant build equal canonical apart from the recycle pools?
# Compares, per pathway and by entity (uuids are re-minted per build), every edge of
# logic_network.csv + boundary_edges.csv and every row of pools.csv / pool_transitions.csv.
# Usage: python3 compare_build.py <canonical build> <variant build>
import csv, os, sys
from collections import Counter
a, b = sys.argv[1:3]
def load(d):
    sid = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
    E = Counter()
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if os.path.exists(f"{d}/{f}"):
            for e in csv.DictReader(open(f"{d}/{f}")):
                E[(sid.get(e["source_id"], "?"), sid.get(e["target_id"], "?"), e["edge_type"], e["pos_neg"])] += 1
    P = Counter((r["pool_id"].rstrip("0123456789"), sid.get(r["node_uuid"], "?"), r["role"], r["is_base"])
                for r in csv.DictReader(open(d + "/pools.csv")))
    T = Counter((r["pool_id"].rstrip("0123456789"), sid.get(r["source_uuid"], "?"), sid.get(r["target_uuid"], "?"),
                 r["reaction_stid"]) for r in csv.DictReader(open(d + "/pool_transitions.csv")))
    return E, P, T
diff = 0
for pw in sorted(p for p in os.listdir(a) if p.startswith("R-HSA-")):
    (ea, pa, ta), (eb, pb, tb) = load(f"{a}/{pw}"), load(f"{b}/{pw}")
    d = [(name, x - y, y - x) for name, x, y in (("edges", ea, eb), ("pools", pa, pb), ("transitions", ta, tb)) if x != y]
    for name, minus, plus in d:
        diff += 1
        kinds = Counter(k[2] if name == "edges" else k[0] for k in (plus + minus).elements())
        print(pw, name, "removed", sum(minus.values()), "added", sum(plus.values()), dict(kinds))
print("pathway tables that differ:", diff)
