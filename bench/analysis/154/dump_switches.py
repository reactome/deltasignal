#!/usr/bin/env python3
"""specs/154 FR-006, dump side: drive the SERVED stack end to end for every
benchmark (pathway, gene, direction) and record each recorded readout's switch
fields. parse -> solve (pins from the specs/140 dump, ko=0 / oe=80) ->
readouts with the solve's switch_components passed back. The engine must be
the feature branch and the catalog the one the dumps used.

Usage: dump_switches.py PINS_JSON RAILS_JSON READOUTS_JSON OUT_TSV [BASE_URL]
"""
import csv, json, sys, urllib.request
from collections import defaultdict

pins_p, rails_p, reads_p, out_p = sys.argv[1:5]
base = sys.argv[5] if len(sys.argv) > 5 else "http://localhost:8080"

def post(ep, body):
    rq = urllib.request.Request(f"{base}/api/{ep}", json.dumps(body).encode(),
                                {"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(rq, timeout=600).read())

pins = {p["pathway_name"]: p for p in json.load(open(pins_p))}
rails = json.load(open(rails_p))
reads = json.load(open(reads_p))
stids_of = defaultdict(dict)   # pathway_name -> ko -> (stable_ids, members)
for p in reads["curator"]:
    for q in p["readouts"]:
        stids_of[p["pathway_name"]][q["ko"]] = (list(q["stable_ids"]), list(q["members"]))

cases = defaultdict(set)       # (pathway, gene, direction) -> set of (axis, ko)
for c in rails["cases"]:
    cases[(c["pathway"], c["gene"], c["direction"])].add((c["axis"], c["key_output"]))
fold_rec = {(c["pathway"], c["gene"], c["direction"], c["key_output"], c["axis"]):
            (c["fold"], c["railed"]) for c in rails["cases"]}

w = csv.writer(open(out_p, "w"), delimiter="\t")
w.writerow(["axis", "pathway", "gene", "direction", "key_output", "fold_recorded",
            "railed_recorded", "fold_served", "switch", "state", "other_state_fold",
            "components_found", "converged"])
nid_cache = {}
solves = 0
by_pw = defaultdict(list)
for (pw, gene, d), kos in cases.items():
    by_pw[pw].append((gene, d, kos))

for pw in sorted(by_pw):
    pq = pins.get(pw)
    if pq is None:
        continue
    if pw not in nid_cache:
        nid_cache[pw] = post("parse", {"pathway_id": pq["dir"]})["network_id"]
    nid = nid_cache[pw]
    pin_of = {q["gene"]: q["pins"] for q in pq["queries"]}
    for gene, d, kos in sorted(by_pw[pw]):
        pu = pin_of.get(gene, [])
        if not pu:
            continue
        act = 0.0 if d == "0" else 80.0
        sol = post("solve", {"network_id": nid, "observations": {u: [act, 0.9] for u in pu}})
        solves += 1
        acts = sol.get("node_activities") or {}
        comps = sol.get("switch_components", [])
        queries = []
        for axis, ko in sorted(kos):
            st, mem = stids_of[pw].get(ko, ([f"R-HSA-{ko}"], []))
            queries.append({"id": f"{axis}|{ko}", "stable_ids": st, "members": mem})
        ro = post("readouts", {"network_id": nid, "activities": acts,
                               "pinned": pu, "queries": queries,
                               "converged": bool(sol.get("converged", False)),
                               "switch_components": comps})
        for e in ro["results"]:
            axis, ko = e["id"].split("|", 1)
            fr = fold_rec.get((pw, gene, d, ko, axis), (None, None))
            w.writerow([axis, pw, gene, d, ko, fr[0], fr[1],
                        e.get("fold"), e.get("switch"), e.get("state"),
                        e.get("other_state_fold"), len(comps),
                        sol.get("converged")])
    print(f"{pw}: done ({solves} solves so far)", flush=True)
print(f"total solves {solves} -> {out_p}")
