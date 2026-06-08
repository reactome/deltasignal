#!/usr/bin/env python3
"""Summarize topology of the shortest path for each wrong case.

Per case, computes the shortest gene→readout path and reports:
  - hops (number of edges)
  - cluster mix: how many AND, OR, and single-input clusters the path crosses
  - OR-branching: at each OR step, how many alternative inputs exist
  - max single-step OR branching factor (used to identify "swallowed by paralog
    cluster" cases — KOs lost in 6+ way OR clusters can't propagate via OR-max)

Useful for distinguishing PATH_EXISTS_SIGNAL_LOST cases by mechanism:
  - long-path attenuation (many AND hops with multiple inputs each)
  - paralog masking (single big OR cluster with many alternatives)
  - cleanly broken (path exists but logically can't propagate, e.g., regulator-
    only path bypassing the main reaction flow)

Usage
-----
  python path_summary.py /tmp/ds_classified.tsv [CATEGORY]
"""
import csv, sys
from collections import defaultdict
from _common import (Networks, load_gene_stids, shortest_path_back,
                     resolve_gene_uuids, resolve_readout_uuids, read_case_dump,
                     CLS_NAME)


def analyze(row, nets, gene_stids):
    pname = row["pathway"]
    inc, _, stids, rev, proxies = nets.get(pname)
    gene_uuids = resolve_gene_uuids(row["gene"], gene_stids, rev)
    ko_uuids = resolve_readout_uuids(row["key_output"], rev, proxies)
    if not gene_uuids or not ko_uuids:
        return None
    path = shortest_path_back(inc, set(gene_uuids), set(ko_uuids), max_depth=30)
    if path is None:
        return None
    or_branches = []
    cluster_counts = defaultdict(int)
    for src, tgt, pn, ao, et in path:
        siblings = inc.get(tgt, [])
        and_n = sum(1 for s in siblings if s[2] == "and")
        or_n = sum(1 for s in siblings if s[2] == "or")
        sg_n = sum(1 for s in siblings if s[2] == "")
        cluster_counts[ao or "single"] += 1
        if ao == "or" and or_n > 1:
            or_branches.append(or_n)
    return {
        "hops": len(path),
        "and": cluster_counts["and"],
        "or":  cluster_counts["or"],
        "single": cluster_counts["single"],
        "or_branches": or_branches,
        "max_or_branch": max(or_branches) if or_branches else 0,
    }


def main(dump_path, category=None):
    rows = read_case_dump(dump_path)
    if category:
        rows = [r for r in rows if r.get("category") == category]
        print(f"Filtered to {len(rows)} rows of category={category}", file=sys.stderr)
    nets = Networks()
    gene_stids = load_gene_stids()

    summary = []
    for r in rows:
        a = analyze(r, nets, gene_stids)
        if a is None:
            continue
        summary.append({**r, **a})

    print(f"\n{'pathway':30s} {'gene':9s} {'dir':3s} {'→':>8s}  exp pred "
          f"hops and  or single max_or_branch")
    for r in summary:
        d = "KO" if r["direction"] == "0" else "OE"
        e = CLS_NAME.get(r["expected"], "?")
        p = CLS_NAME.get(r["predicted"], "?")
        print(f"{r['pathway'][:28]:30s} {r['gene']:9s} {d:3s} "
              f"{r['key_output']:>8s}  {e:4s} {p:4s} "
              f"{r['hops']:3d}  {r['and']:3d}  {r['or']:3d}  {r['single']:3d}  "
              f"{r['max_or_branch']:3d}")

    # Aggregate stats
    if summary:
        print("\n=== Aggregates ===")
        avg = lambda k: sum(r[k] for r in summary) / len(summary)
        print(f"  Median hops: {sorted(r['hops'] for r in summary)[len(summary)//2]}")
        print(f"  Mean hops:   {avg('hops'):.1f}")
        print(f"  Mean AND clusters per path:    {avg('and'):.2f}")
        print(f"  Mean OR  clusters per path:    {avg('or'):.2f}")
        with_or = [r for r in summary if r["max_or_branch"]]
        if with_or:
            print(f"  Cases crossing ≥1 OR cluster: {len(with_or)}")
            print(f"    Max OR branching (largest): "
                  f"{max(r['max_or_branch'] for r in with_or)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr); sys.exit(1)
    main(*sys.argv[1:])
