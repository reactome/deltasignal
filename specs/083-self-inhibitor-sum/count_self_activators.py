# specs/083: reactions with SELF-CONTAINED ACTIVATORS: a positive regulator or catalyst
# whose entity contains one of the reaction's `input` entities (containment.csv) - a
# TF:gene complex regulating that gene's transcription. Each multiplies the input's fold
# again under the AND product. The mirror of specs/022's self-contained inhibitors.
# Usage: python3 count_self_activators.py <build dir> [<results dir>]
import csv, os, sys
from collections import Counter, defaultdict
b = sys.argv[1]
tot = Counter(); per = Counter(); ex = []
by_pw_rx = defaultdict(set)
for pw in sorted(p for p in os.listdir(b) if p.startswith("R-HSA-")):
    d = f"{b}/{pw}"
    full = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
    cont = defaultdict(set)
    for r in csv.DictReader(open(d + "/containment.csv")):
        cont[r["stable_id"]].add(r["contains_stable_id"])
    rx = defaultdict(lambda: {"in": set(), "act": []})
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if not os.path.exists(f"{d}/{f}"): continue
        for e in csv.DictReader(open(f"{d}/{f}")):
            if e["pos_neg"] != "pos": continue
            s = full.get(e["source_id"], "")
            if e["edge_type"] == "input": rx[e["target_id"]]["in"].add(s)
            elif e["edge_type"] in ("regulator", "catalyst"): rx[e["target_id"]]["act"].append(s)
    for t, x in rx.items():
        selfs = [a for a in x["act"] if any(i and i != a and i in cont.get(a, ()) for i in x["in"])]
        if selfs:
            tot["reaction nodes"] += 1; tot["self-contained activator slots"] += len(selfs)
            k = full.get(t, t).split("::")[0]
            if k not in by_pw_rx[pw]:
                by_pw_rx[pw].add(k); per[pw] += 1
                ex.append((pw, k, len({a.split("::")[0] for a in selfs})))
print(dict(tot), "distinct reactions", sum(per.values()), "pathways", len(per))
for pw, n in per.most_common(10): print("  ", pw, n)
for x in sorted(ex, key=lambda x: -x[2])[:12]: print("  ", x)
