#!/usr/bin/env python3
"""Classify each wrong case in a per-case dump by root-cause category.

Categories
----------
  GENE_NOT_IN_PATHWAY         — Perturbed gene has no UUID in the pathway network.
  READOUT_NOT_IN_PATHWAY      — Readout entity has no UUID and no proxy.
  NO_PATH_GENE_TO_READOUT     — Both endpoints exist, but no directed path.
  PATH_EXISTS_DIRECTION_FLIPPED  — Path exists, model prediction has the OPPOSITE
                                   sign of the expected change (DOWN↔UP).
  PATH_EXISTS_FALSE_POSITIVE  — Expected NORM, path exists, model predicted change.
  PATH_EXISTS_SIGNAL_LOST     — Expected non-NORM, path exists, model predicted NORM.
  PATH_EXISTS_OTHER_…         — Any remaining residual exp/pred combination.

For Mit_G1 and S_Phase cases the JOINT network is used (matching the
benchmark's selective-joint configuration). All other pathways use the
single-pathway network.

Usage
-----
  python classify_failures.py /tmp/ds_best_cases.tsv [out.tsv]
"""
import csv, sys
from collections import defaultdict
from _common import (Networks, load_gene_stids, bfs_upstream,
                     resolve_gene_uuids, resolve_readout_uuids,
                     read_case_dump, CLS_NAME)


def classify(row, nets, gene_stids):
    pname = row["pathway"]
    inc, _, _, rev, proxies = nets.get(pname)
    gene_uuids = resolve_gene_uuids(row["gene"], gene_stids, rev)
    ko_uuids = resolve_readout_uuids(row["key_output"], rev, proxies)
    if not gene_uuids:
        return "GENE_NOT_IN_PATHWAY"
    if not ko_uuids:
        return "READOUT_NOT_IN_PATHWAY"
    upstream = bfs_upstream(inc, set(ko_uuids), max_depth=30)
    if not (set(gene_uuids) & upstream):
        return "NO_PATH_GENE_TO_READOUT"
    e, p = row["expected"], row["predicted"]
    if (e, p) in (("0", "2"), ("2", "0")):
        return "PATH_EXISTS_DIRECTION_FLIPPED"
    if e == "1":
        return "PATH_EXISTS_FALSE_POSITIVE"
    if p == "1":
        return "PATH_EXISTS_SIGNAL_LOST"
    return f"PATH_EXISTS_OTHER_E{e}_P{p}"


def main(dump_path, out_path=None):
    rows = read_case_dump(dump_path)
    gene_stids = load_gene_stids()
    nets = Networks()
    wrong = [r for r in rows if r["predicted"] != r["expected"]]
    print(f"Total cases: {len(rows)}  Wrong: {len(wrong)}  "
          f"Acc: {(len(rows) - len(wrong)) * 100 / len(rows):.2f}%",
          file=sys.stderr)

    classified = []
    cnt = defaultdict(int)
    for r in wrong:
        cat = classify(r, nets, gene_stids)
        cnt[cat] += 1
        classified.append({**r, "category": cat})

    print("\n=== Failure root-cause categories ===")
    for k, v in sorted(cnt.items(), key=lambda x: -x[1]):
        print(f"  {k:35s} {v}")

    # Per-pathway breakdown
    by_pw = defaultdict(lambda: defaultdict(int))
    for r in classified:
        by_pw[r["pathway"]][r["category"]] += 1
    print("\n=== Per-pathway × category ===")
    for pw in sorted(by_pw, key=lambda p: -sum(by_pw[p].values())):
        cs = by_pw[pw]
        total = sum(cs.values())
        bits = " ".join(f"{k.replace('PATH_EXISTS_', 'P_'):12s}={v}"
                        for k, v in sorted(cs.items(), key=lambda x: -x[1]))
        print(f"  {pw[:38]:40s} ({total:3d}) {bits}")

    if out_path:
        with open(out_path, "w") as f:
            cols = ["pathway", "gene", "direction", "key_output",
                    "predicted", "expected", "method", "pred_ui", "category"]
            f.write("\t".join(cols) + "\n")
            for r in classified:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        print(f"\nWrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: classify_failures.py <case_dump.tsv> [out.tsv]",
              file=sys.stderr)
        sys.exit(1)
    main(*sys.argv[1:])
