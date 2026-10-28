#!/usr/bin/env python3
"""Hop distance from a decryptM case's pinned targets to its readout nodes (graph only; specs/116).

Resolves pins and readouts exactly as `decryptm_bench.py` does (`pin_occurrences`, `readout_uuids`), then takes
the shortest path in `build_adjacency` from any pinned node to any unpinned readout node; blank if no route.

Usage: decryptm_hops.py BUILD CASES_TSV PHOSPHO_CSV OUT_TSV
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
from decryptm_bench import readout_uuids  # noqa: E402
from decryptm_coverage import RESIDUE, neo4j_rows  # noqa: E402
from korkut_hops import hops_from  # noqa: E402


def main():
    build, cases_tsv, phospho, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
    site_stids = collections.defaultdict(set)
    for st, acc, _gene, coord, mod in neo4j_rows(phospho):
        if mod in RESIDUE:
            site_stids[(acc, int(coord))].add(st)
    cases = [c for c in csv.DictReader(open(cases_tsv), delimiter="\t") if c["status"] == "scored"]
    gene_to_stids = bm.neo4j_gene_to_stids(sorted({g for c in cases for g in c["targets"].split(";")}))
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for c in cases:
        by[c["pathway_id"]][(c["drug"], c["targets"])].append(c)
    rows = []
    for pid in sorted(by):
        pdir = build / pid
        s2u = bm.load_stid_to_uuids(pdir)
        adj = bm.build_adjacency(pdir)
        for (drug, targets), group in sorted(by[pid].items()):
            g2u = {g: [u for s in gene_to_stids.get(g, []) for u in s2u.get(s, [])] for g in targets.split(";")}
            pins, _ = bm.pin_occurrences(pdir, g2u, gene_to_stids)
            pinned = set().union(*pins.values()) if pins else set()
            dist = hops_from(adj, pinned)
            for c in group:
                ro = readout_uuids(site_stids.get((c["uniprot"], int(c["position"])), ()), s2u)
                d = [dist[u] for u in ro if u not in pinned and u in dist]
                rows.append([pid, drug, c["uniprot"], c["position"], min(d) if d else ""])
        print(pid, flush=True)
    with open(out, "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pathway_id", "drug", "uniprot", "position", "hops"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
