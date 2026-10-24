"""Census: reaction copies whose protein content does not balance when a sibling copy does."""
import csv, sys, collections, subprocess, json, os
B = sys.argv[1]
pws = sorted(p for p in os.listdir(B) if p.startswith("R-HSA-") and os.path.isfile(f"{B}/{p}/logic_network.csv"))
leaves_all = set(); data = {}
for pw in pws:
    nodes = {r["uuid"]: r for r in csv.DictReader(open(f"{B}/{pw}/nodes.csv"))}
    ins = collections.defaultdict(list); outs = collections.defaultdict(list); rx = {}
    for e in csv.DictReader(open(f"{B}/{pw}/logic_network.csv")):
        if e["edge_type"] == "input": ins[e["target_id"]].append(e["source_id"]); rx[e["target_id"]] = e["edge_reaction_id"]
        elif e["edge_type"] == "output": outs[e["source_id"]].append(e["target_id"]); rx[e["source_id"]] = e["edge_reaction_id"]
    lv = lambda u, nodes=nodes: set((nodes.get(u, {}).get("member_leaves") or nodes.get(u, {}).get("diagram_entity_id") or "").split("|")) - {""}
    data[pw] = (ins, outs, rx, lv)
    for u in nodes.values(): leaves_all |= set(u["member_leaves"].split("|")) | {u["diagram_entity_id"]}
leaves_all.discard("")
ref = {}
for line in open(sys.argv[3]).read().splitlines()[1:]:
    a, b = [x.strip().strip('"') for x in line.split(",", 1)]
    ref[a] = b
ref = {k: v for k, v in ref.items() if k in leaves_all}
tot_e = tot_bad = 0; per_pw = {}; examples = []
for pw, (ins, outs, rx, lv) in data.items():
    groups = collections.defaultdict(list)
    for node, r in rx.items():
        if ins[node] and outs[node]: groups[r].append(node)
    bad = 0; rset = set()
    for r, copies in groups.items():
        inref = {c: {ref[l] for u in ins[c] for l in lv(u) if l in ref} for c in copies}
        for c in copies:
            for t in outs[c]:
                tot_e += 1
                o = {ref[l] for l in lv(t) if l in ref}
                created = o - inref[c]
                if not created: continue
                # another copy of the reaction that outputs t supplies it
                if any(created <= inref[c2] for c2 in copies if c2 != c and t in outs[c2]):
                    bad += 1; rset.add(r)
                    if len(examples) < 4000: examples.append((pw, r, c, t, sorted(created)))
    per_pw[pw] = (len(rset), bad); tot_bad += bad
print(f"leaves mapped {len(ref)}/{len(leaves_all)}; output edges {tot_e}; creating a protein a sibling copy supplies: {tot_bad} in {sum(v[0] for v in per_pw.values())} reactions; pathways {sum(1 for v in per_pw.values() if v[0])}/{len(per_pw)}")
for pw, v in sorted(per_pw.items(), key=lambda kv: -kv[1][1])[:20]: print(pw, v)
json.dump(examples, open(sys.argv[2], "w"))
