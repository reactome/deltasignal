#!/usr/bin/env python3
"""specs/147 P4: census the diagram-only Interaction links for the benchmark
pathways. Usage: census.py DIAGRAM_DIR OUT_TSV"""
import collections
import csv
import json
import sys
from pathlib import Path

DIA, OUT = Path(sys.argv[1]), sys.argv[2]
HERE = Path(__file__).resolve()
pathways = [r for r in csv.reader((l for l in open(HERE.parents[3] / 'bench/catalog_pathways.tsv')
                                   if not l.startswith('#')), delimiter='\t') if r and r[0] != 'id']
rows, total = [], collections.Counter()
covered = 0
for stid, name in pathways:
    p = DIA / f"{stid}.json"
    if not p.exists():
        continue
    covered += 1
    d = json.load(open(p))
    byid = {n["id"]: n for n in d.get("nodes", [])}
    procs = {n["id"] for n in d.get("nodes", []) if n.get("renderableClass") in ("ProcessNode", "EncapsulatedNode")}
    for l in d.get("links") or []:
        rc = l.get("renderableClass")
        # release 98 re-renders some drawn connection arrows as FlowLine
        if rc not in ("Interaction", "FlowLine"):
            continue
        for i in [q.get("id") for q in (l.get("inputs") or [])]:
            for o in [q.get("id") for q in (l.get("outputs") or [])]:
                ni, no = byid.get(i) or {}, byid.get(o) or {}
                kind = ("pathway->entity" if i in procs and o not in procs else
                        "entity->pathway" if o in procs and i not in procs else
                        "pathway->pathway" if i in procs and o in procs else "entity->entity")
                total[kind] += 1
                rows.append({"pathway": stid, "pathway_name": name, "kind": kind, "link_class": rc,
                             "from": ni.get("displayName", "?"), "from_id": ni.get("reactomeId", ""),
                             "to": no.get("displayName", "?"), "to_id": no.get("reactomeId", "")})
with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, ["pathway", "pathway_name", "kind", "link_class", "from", "from_id", "to", "to_id"], delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(f"{covered} diagrams; {len(rows)} interaction arrows; {dict(total)}")
