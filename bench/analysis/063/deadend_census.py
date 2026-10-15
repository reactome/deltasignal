#!/usr/bin/env python3
"""specs/063 census: binding reactions whose complex is a dead end.

A candidate is a reaction R with
- at least two non-cofactor entity inputs (a binding),
- exactly one non-cofactor output C,
- C has no out-edge except into dissociation sinks (nothing in the pathway
  uses it).
Its partners are R's non-cofactor inputs with at least one other use (an
out-edge into anything but R or a dissociation sink). The rule under study
adds C -| P (depletion) for each such partner.

Usage: deadend_census.py BUILD_DIR [OUT_TSV]
"""
import collections
import csv
import sys
from pathlib import Path

ENTITY_SKIP = {"reaction", "dissociation_sink"}


def load(pdir):
    kind = {n["uuid"]: n["node_kind"] for n in csv.DictReader(open(pdir / "nodes.csv"))}
    dia = {n["uuid"]: n["diagram_entity_id"] for n in csv.DictReader(open(pdir / "nodes.csv"))}
    cof = set()
    if (pdir / "cofactors.csv").exists():
        stid = {a: b for a, b in csv.reader(open(pdir / "stid_to_uuid_mapping.csv"))}
        cs = {r[0] for r in csv.reader(open(pdir / "cofactors.csv"))}
        cof = {u for u, s in stid.items() if s in cs}
    out, inc = collections.defaultdict(list), collections.defaultdict(list)
    for f in ("logic_network.csv", "boundary_edges.csv"):
        if (pdir / f).exists():
            for e in csv.DictReader(open(pdir / f)):
                out[e["source_id"]].append(e)
                inc[e["target_id"]].append(e)
    return kind, dia, cof, out, inc


def candidates(pdir):
    kind, dia, cof, out, inc = load(pdir)
    rows = []
    for r, k in sorted(kind.items()):
        if k != "reaction":
            continue
        ins = sorted({e["source_id"] for e in inc[r] if e["edge_type"] == "input"} - cof)
        outs = sorted({e["target_id"] for e in out[r] if e["edge_type"] == "output"} - cof)
        if len(ins) < 2 or len(outs) != 1:
            continue
        c = outs[0]
        if kind.get(c) in ENTITY_SKIP or any(kind.get(e["target_id"]) != "dissociation_sink" for e in out[c]):
            continue
        partners = [p for p in ins
                    if any(e["target_id"] != r and kind.get(e["target_id"]) != "dissociation_sink" for e in out[p])]
        for p in partners:
            rows.append({"pathway_id": pdir.name, "reaction": r, "complex": c, "complex_entity": dia.get(c, ""),
                         "partner": p, "partner_entity": dia.get(p, ""), "partner_kind": kind.get(p, "")})
    return rows


def main():
    build = Path(sys.argv[1])
    rows = [x for p in sorted(build.iterdir()) if (p / "nodes.csv").exists() for x in candidates(p)]
    print(f"{len(rows)} sequestration edges, {len({(r['pathway_id'], r['complex']) for r in rows})} dead-end "
          f"complexes, {len({r['pathway_id'] for r in rows})} pathways")
    for pw, n in collections.Counter(r["pathway_id"] for r in rows).most_common(12):
        print(f"  {pw}\t{n}")
    if len(sys.argv) > 2:
        with open(sys.argv[2], "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
