#!/usr/bin/env python3
"""Trace the shortest path from a perturbed gene to a readout node.

For each (pathway, gene, key_output) tuple, walks the BFS shortest path
backward from the readout to any gene UUID and prints, for every edge:
  - source UUID (truncated) + Reactome stid (if known)
  - polarity (pos/neg), cluster type (and/or/single), edge category (input/
    catalyst/regulator/output/depletion/etc.)
  - branching factor of the target's incoming siblings (how many alternative
    AND inputs or OR alternatives exist at that step)

The shortest path is a sample — the model evaluates every parallel path in
parallel, so signal strength at the readout is determined by ALL upstream
paths combined, not just this one. The trace just makes it easy to eyeball
where signal might be diluted or lost.

Usage
-----
  python trace_paths.py PATHWAY GENE KEY_OUTPUT
  python trace_paths.py --from-dump /tmp/ds_classified.tsv PATH_EXISTS_SIGNAL_LOST
"""
import csv, sys
from collections import defaultdict
from _common import (Networks, load_gene_stids, shortest_path_back,
                     resolve_gene_uuids, resolve_readout_uuids, read_case_dump)


def trace(pname, gene, key_output, nets, gene_stids, max_depth=30):
    inc, _, stids, rev, proxies = nets.get(pname)
    gene_uuids = resolve_gene_uuids(gene, gene_stids, rev)
    ko_uuids = resolve_readout_uuids(key_output, rev, proxies)
    if not gene_uuids:
        print(f"  ! gene {gene} has no UUIDs in {pname}")
        return
    if not ko_uuids:
        print(f"  ! readout R-HSA-{key_output} has no UUIDs in {pname}")
        return
    path = shortest_path_back(inc, set(gene_uuids), set(ko_uuids), max_depth)
    if path is None:
        print(f"  ! no path {gene} → R-HSA-{key_output} in {pname}")
        return
    print(f"\n>> {pname}  {gene} → R-HSA-{key_output}  ({len(path)} edges)")
    src0 = path[0][0]
    print(f"   [0] {src0[:8]} {stids.get(src0, '?')}")
    for i, (src, tgt, pn, ao, et) in enumerate(path):
        siblings = inc.get(tgt, [])
        and_n = sum(1 for s in siblings if s[2] == "and")
        or_n  = sum(1 for s in siblings if s[2] == "or")
        sg_n  = sum(1 for s in siblings if s[2] == "")
        cluster = ao or "single"
        print(f"       │  {pn:3s} cluster={cluster:6s} type={et}  "
              f"siblings_at_tgt: and={and_n} or={or_n} single={sg_n}")
        print(f"   [{i+1}] {tgt[:8]} {stids.get(tgt, '?')}")


def from_dump(dump_path, category, limit=None):
    rows = [r for r in read_case_dump(dump_path) if r["category"] == category]
    if limit:
        rows = rows[:int(limit)]
    nets = Networks()
    gene_stids = load_gene_stids()
    for r in rows:
        trace(r["pathway"], r["gene"], r["key_output"], nets, gene_stids)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr); sys.exit(1)
    if sys.argv[1] == "--from-dump":
        # --from-dump DUMP CATEGORY [LIMIT]
        from_dump(*sys.argv[2:])
    else:
        nets = Networks()
        gene_stids = load_gene_stids()
        trace(sys.argv[1], sys.argv[2], sys.argv[3], nets, gene_stids)
