#!/usr/bin/env python3
"""specs/140 FR-003, bench side: dump the pins the BENCH resolves for every
benchmark (pathway, gene), via the very functions the benchmarks run
(`pin_occurrences`), so the engine's port can be held equal to them.

Writes JSON:
  [{"pathway_id", "pathway_name", "dir",
    "queries": [{"gene", "stable_ids", "pins", "basis", "occurrences"}]}]

`pins` are sorted (the gate compares sets); `basis` is derived the way
resolve_pins reports it: "root" when the gene has a root occurrence,
"root_cycle" when the specs/024 fallback produced the pins, "none" otherwise.

Run with DS_CATALOG_ROOT pointing at the catalog under test (the frozen build
for the registered parity run):
  DS_CATALOG_ROOT=~/deltasignal-catalogs/current \
    python bench/analysis/140/dump_pins.py bench/analysis/140/bench_pins.json
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from benchmark_vs_mpbiopath import (  # noqa: E402
    CATALOG_ROOT, GENE_NAME_CORRECTIONS, MPBIO_ROOT, apply_catalog_ids,
    find_pathway_dir, load_curator, load_stid_to_uuids, neo4j_gene_to_stids,
    parse_perturbation_columns, pin_occurrences, root_occurrences,
)


def pathway_list():
    with open(MPBIO_ROOT / "pathway_list.tsv") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx_id, idx_name = header.index("pathway_id"), header.index("pathway_name")
        pathways = [(c[idx_id], c[idx_name]) for c in
                    (line.rstrip("\n").split("\t") for line in f)
                    if len(c) > max(idx_id, idx_name)]
    return apply_catalog_ids(pathways,
                             Path(__file__).resolve().parents[2] / "catalog_pathways.tsv")


def main():
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "bench_pins.json"
    print(f"catalog: {CATALOG_ROOT}", flush=True)
    dumped, missing = [], 0
    for pid, pname in pathway_list():
        pathway_dir = find_pathway_dir(pid)
        if pathway_dir is None:
            missing += 1
            continue
        curator = load_curator(pname)
        if curator is None:
            continue
        header, _rows = curator
        gene_names = sorted({g for g, _, _ in parse_perturbation_columns(header)})
        lookup = {g: GENE_NAME_CORRECTIONS.get((pname, g), g) for g in gene_names}
        found = neo4j_gene_to_stids(sorted(set(lookup.values())))
        gene_to_stids = {g: found.get(lookup[g], []) for g in gene_names}

        stid_to_uuids = load_stid_to_uuids(pathway_dir)
        gene_to_uuids = {g: [u for sid in gene_to_stids[g] for u in stid_to_uuids.get(sid, [])]
                         for g in gene_names}
        occurrences = {g: list(dict.fromkeys(us)) for g, us in gene_to_uuids.items()}
        resolved, indeg = pin_occurrences(pathway_dir, gene_to_uuids, gene_to_stids)

        queries = []
        for g in gene_names:
            occ = occurrences[g]
            pins = resolved[g]
            roots = root_occurrences(occ, indeg)
            basis = ("none" if not occ or not pins else
                     "root" if roots else "root_cycle")
            queries.append({"gene": g, "stable_ids": gene_to_stids[g],
                            "pins": sorted(set(pins)), "basis": basis,
                            "occurrences": len(occ)})
        dumped.append({"pathway_id": pid, "pathway_name": pname,
                       "dir": pathway_dir.name, "queries": queries})
        print(f"  {pname}: {len(queries)} genes", flush=True)

    out_path.write_text(json.dumps(dumped, indent=1))
    n_q = sum(len(p["queries"]) for p in dumped)
    print(f"wrote {out_path}: {len(dumped)} pathways, {n_q} gene queries "
          f"({missing} pathways without a network)", flush=True)


if __name__ == "__main__":
    main()
