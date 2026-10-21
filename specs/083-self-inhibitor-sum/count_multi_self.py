# specs/083: reactions with two or more SELF-CONTAINED inhibitors: negative (inhibitor,
# not depletion) in-edges whose source entity contains one of the reaction's activating
# inputs' entity (containment.csv), the specs/022 pair. Counted per reaction node, then by
# entity in each pathway.
# Usage: python3 count_multi_self.py <build dir>
import csv, os, sys
from collections import Counter, defaultdict
b = sys.argv[1]
tot = Counter(); per = Counter(); examples = []
for pw in sorted(p for p in os.listdir(b) if p.startswith("R-HSA-")):
    d = f"{b}/{pw}"
    full = {r["uuid"]: r["stable_id"] for r in csv.DictReader(open(d + "/stid_to_uuid_mapping.csv"))}
    sid = {u: s.split("::")[0] for u, s in full.items()}
    cont = defaultdict(set)
    for r in csv.DictReader(open(d + "/containment.csv")):
        cont[r["stable_id"]].add(r["contains_stable_id"])
    ins = defaultdict(lambda: {"pos": set(), "neg": []})
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if not os.path.exists(f"{d}/{f}"): continue
        for e in csv.DictReader(open(f"{d}/{f}")):
            t = e["target_id"]
            if e["pos_neg"] == "pos": ins[t]["pos"].add(full.get(e["source_id"], ""))
            elif e["edge_type"] != "depletion": ins[t]["neg"].append(full.get(e["source_id"], ""))
    seen = set()
    for t, x in ins.items():
        selfs = [n for n in x["neg"] if any(a and a in cont.get(n, ()) for a in x["pos"])]
        if len(selfs) >= 2:
            tot["reaction nodes"] += 1; tot["self inhibitor slots"] += len(selfs)
            base = tuple(sorted({n.split("::")[0] for n in selfs}))
            key = (sid.get(t, t), base)
            if key not in seen:
                seen.add(key); per[pw] += 1
                examples.append((pw, sid.get(t, t), len(base)))
print(dict(tot), "pathways", len(per))
for pw, n in per.most_common(12): print("  ", pw, n)
for x in sorted(examples, key=lambda x: -x[2])[:12]: print("  ", x)
