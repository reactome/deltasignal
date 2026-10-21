# specs/084: closed carrier-loop free forms with NO root copy of their own entity (so
# closed_twin gives them nothing) whose contained entities DO have root copies in the
# pathway. The rule would feed each such form from those root copies, keeping only the
# maximal ones (a root entity contained in another contained root entity is skipped).
# Usage: python3 count_component_twins.py <build dir> [<results dir>]
import csv, os, sys
from collections import Counter, defaultdict
b = sys.argv[1]; res = sys.argv[2] if len(sys.argv) > 2 else None
tot = Counter(); ex = []; touched = set()
for pw in sorted(p for p in os.listdir(b) if p.startswith("R-HSA-")):
    d = f"{b}/{pw}"
    if not os.path.exists(d + "/carrier_loops.csv"): continue
    full = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
    cont = defaultdict(set)
    for r in csv.DictReader(open(d + "/containment.csv")): cont[r["stable_id"]].add(r["contains_stable_id"])
    cof = set()
    if os.path.exists(d + "/cofactors.csv"):
        cof = {list(r.values())[0] for r in csv.DictReader(open(d + "/cofactors.csv"))}
    hasin = set()
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if os.path.exists(f"{d}/{f}"):
            for e in csv.DictReader(open(f"{d}/{f}")):
                if e["edge_type"] != "recycle": hasin.add(e["target_id"])
    rows = list(csv.DictReader(open(d + "/carrier_loops.csv")))
    listed = {r["free_uuid"] for r in rows}
    roots = defaultdict(list)
    for u, s in full.items():
        if u not in hasin and u not in listed: roots[s].append(u)
    closed = sorted({r["free_uuid"] for r in rows if r["closed"] == "True"})
    for f in closed:
        s = full.get(f, "")
        if roots.get(s): continue                     # closed_twin already feeds it
        inside = {c for c in cont.get(s, ()) if c != s and c not in cof and roots.get(c)}
        maximal = {c for c in inside if not any(c != o and c in cont.get(o, ()) for o in inside)}
        if not maximal: tot["closed forms with no root of any kind"] += 1; continue
        tot["closed forms fed by component roots"] += 1
        tot["component-root edges"] += sum(len(roots[c]) for c in maximal)
        touched.add((pw, s.split("::")[0]))
        ex.append((pw, s.split("::")[0], sorted(c.split("::")[0] for c in maximal)[:4]))
print(dict(tot), "distinct (pathway, entity)", len(touched), "pathways", len({p for p, _ in touched}))
for x in ex[:20]: print("  ", x)
