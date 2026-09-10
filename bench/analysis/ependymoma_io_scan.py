#!/usr/bin/env python3
"""Feasibility scan for Marija's approach #4 on the Taylor PFA hypoxia DEGs.

For each Reactome pathway in the catalog, classify every DEG that appears as:
  - input-only  : feeds reactions but is not produced in-pathway (a root/source)
  - output-only : produced by a reaction but feeds nothing (a terminal/sink)
  - both        : intermediate
Then find pathways that contain input-only DEGs AND output-only DEGs with a
directed path between them -- these are the testable "pin inputs, predict
outputs, compare to measured direction" cases.

Resolves gene->stids from Reactome Neo4j (same source as build_gene_cache.py).
"""
import csv, json, os, sys
from collections import defaultdict, deque
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common as C

# Persistent home (survives reboot / /tmp cleanup); was the session scratchpad.
SCRATCH = Path("/home/awright/proteomics_analysis/taylor-lab/deltasignal_scan")
DEGS = SCRATCH / "degs.csv"
CACHE = SCRATCH / "deg_gene_to_stids.json"


def load_one(dirname):
    """Self-contained loader tolerant of the current 7-col logic_network.csv
    (source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id)."""
    incoming = defaultdict(list)
    outgoing = defaultdict(list)
    with open(C.CATALOG / dirname / "logic_network.csv") as fh:
        rd = csv.reader(fh)
        next(rd, None)
        for r in rd:
            if len(r) < 5:
                continue
            src, tgt, pn, ao, et = r[0], r[1], r[2], r[3], r[4]
            incoming[tgt].append((src, pn, ao, et))
            outgoing[src].append((tgt, pn, ao, et))
    rev = defaultdict(list)
    for row in csv.DictReader(open(C.CATALOG / dirname / "stid_to_uuid_mapping.csv")):
        rev[row["stable_id"]].append(row["uuid"])
    return incoming, outgoing, dict(rev)


def load_degs():
    rows = list(csv.DictReader(open(DEGS)))
    return {r["gene"]: {"log2FC": float(r["log2FC"]), "dir": r["dir"]} for r in rows}


def build_cache(genes):
    if CACHE.exists():
        return {k: set(v) for k, v in json.load(open(CACHE)).items()}
    from py2neo import Graph
    g = Graph("bolt://localhost:7687", auth=("neo4j", "reactome"))
    rows = g.run(
        "UNWIND $names AS gn "
        "MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity) "
        "WHERE gn IN re.geneName "
        "RETURN gn AS gene, COLLECT(DISTINCT pe.stId) AS stids",
        names=list(genes)).data()
    m = {r["gene"]: list(r["stids"]) for r in rows}
    json.dump(m, open(CACHE, "w"))
    return {k: set(v) for k, v in m.items()}


def forward_reachable(outgoing, sources, targets, max_depth=40):
    """Which target uuids are forward-reachable from any source uuid."""
    tset = set(targets)
    seen = set(sources)
    hit = set()
    q = deque((s, 0) for s in sources)
    while q:
        u, d = q.popleft()
        if u in tset:
            hit.add(u)
        if d >= max_depth:
            continue
        for tgt, _, _, _ in outgoing.get(u, []):
            if tgt not in seen:
                seen.add(tgt)
                q.append((tgt, d + 1))
    return hit


def main():
    degs = load_degs()
    gene_stids = build_cache(set(degs))
    mapped = {g for g in degs if gene_stids.get(g)}
    print(f"DEGs: {len(degs)} total; {len(mapped)} map to >=1 Reactome PhysicalEntity stId "
          f"({len(degs)-len(mapped)} unmapped)\n")

    catalog_dirs = sorted(p.name for p in C.CATALOG.iterdir()
                          if p.is_dir() and (p / "logic_network.csv").exists())
    print(f"Scanning {len(catalog_dirs)} catalog pathways...\n")

    summary = []  # per pathway with both roles
    for d in catalog_dirs:
        incoming, outgoing, rev = load_one(d)
        # classify each DEG present in this pathway
        input_only = {}   # gene -> [uuids]
        output_only = {}
        both = set()
        for g in mapped:
            uuids = C.resolve_gene_uuids(g, gene_stids, rev)
            if not uuids:
                continue
            io = oo = False
            in_uuids, out_uuids = [], []
            for u in uuids:
                has_in = u in incoming
                has_out = u in outgoing
                if has_out and not has_in:
                    io = True; in_uuids.append(u)
                elif has_in and not has_out:
                    oo = True; out_uuids.append(u)
                elif has_in and has_out:
                    pass
            if io and not oo:
                input_only[g] = in_uuids
            elif oo and not io:
                output_only[g] = out_uuids
            elif io or oo:
                both.add(g)
        if not (input_only and output_only):
            continue
        # directed input-DEG -> output-DEG reachability
        src_uuids = [u for us in input_only.values() for u in us]
        tgt_uuids = {u: g for g, us in output_only.items() for u in us}
        hit = forward_reachable(outgoing, src_uuids, set(tgt_uuids))
        reached_genes = sorted({tgt_uuids[u] for u in hit})
        summary.append({
            "pathway": d,
            "n_input_only": len(input_only),
            "n_output_only": len(output_only),
            "input_genes": sorted(input_only),
            "output_genes": sorted(output_only),
            "reachable_output_genes": reached_genes,
            "n_testable_pairs_genes": len(reached_genes),
        })

    summary.sort(key=lambda s: (-s["n_testable_pairs_genes"], -s["n_output_only"]))
    print(f"Pathways with BOTH input-only and output-only DEGs: {len(summary)}")
    tot_testable = sum(s["n_testable_pairs_genes"] for s in summary)
    print(f"Pathways with >=1 DIRECTED input-DEG -> output-DEG path: "
          f"{sum(1 for s in summary if s['n_testable_pairs_genes'])}")
    print(f"Total testable output-DEG predictions (gene-level, path-backed): {tot_testable}\n")

    print("=== Top pathways (by testable input->output DEG predictions) ===")
    for s in summary[:20]:
        print(f"\n* {s['pathway']}")
        print(f"    input-only DEGs ({s['n_input_only']}): {', '.join(s['input_genes'][:12])}"
              + (" ..." if s['n_input_only'] > 12 else ""))
        print(f"    output-only DEGs ({s['n_output_only']}): {', '.join(s['output_genes'][:12])}"
              + (" ..." if s['n_output_only'] > 12 else ""))
        if s["reachable_output_genes"]:
            print(f"    >>> REACHABLE output DEGs from input DEGs: "
                  f"{', '.join(s['reachable_output_genes'])}")

    json.dump(summary, open(SCRATCH / "io_scan_result.json", "w"), indent=2)
    print(f"\nFull result -> {SCRATCH/'io_scan_result.json'}")


if __name__ == "__main__":
    main()
