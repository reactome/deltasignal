#!/usr/bin/env python3
"""Did the PATHWAY change under the ground truth?

The curator predictions are from 2019. The networks we solve are Release97.
Where Reactome recurated a pathway in between, a disagreement with the curator
is not necessarily our error: the curator answered about the structure as it
then was, and we answer about the structure as it is now. Both can be right.

MP-BioPath shipped the networks it actually used, as dbId edge lists — and
dbIds are stable across releases where stIds are not. So the same reachability
question can be put to the 2019 structure and to ours, and the disagreements
sorted:

  lost_since_2019   the 2019 network connects gene to readout, ours does not.
                    A curator "this changes" reflects a route Reactome has
                    since removed. Our NO_CHANGE may be right about today.

  gained_since_2019 ours connects them, the 2019 network did not. A curator
                    "no change" reflects the absence of a route we now have.
                    Our change-call may be right about today.

  same_structure    both agree on reachability; the disagreement is ours to
                    explain.

This does not prove the curator wrong — recuration can also correct an error,
and a route existing is not a route mattering. It separates "we disagree with
a 2019 fact" from "we disagree about the same structure", which are different
problems.
"""

from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

MPBIO = Path.home() / "gitroot" / "mp-biopath-pathways"


def load_mpbio_network(pathway_name: str):
    """dbId adjacency for the network MP-BioPath actually ran."""
    f = MPBIO / "pathways" / f"{pathway_name}.tsv"
    if not f.exists():
        return None
    adj = collections.defaultdict(set)
    with f.open() as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            adj[cols[0].strip()].add(cols[1].strip())
    return adj


def reachable(adj, sources, cap=200000):
    seen = set(sources)
    stack = list(sources)
    while stack:
        u = stack.pop()
        for v in adj.get(u, ()):
            if v not in seen:
                seen.add(v)
                stack.append(v)
                if len(seen) > cap:
                    return seen
    return seen


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import benchmark_vs_mpbiopath as B
    from py2neo import Graph
    graph = Graph(B.NEO4J_URL, auth=(B.NEO4J_USER, B.NEO4J_PASSWORD))

    with args.cases.open() as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("valid") == "1"]
    wrong = [r for r in rows if r["predicted"] != r["expected"]]

    genes = sorted({r["gene"] for r in wrong})
    gene_dbids = collections.defaultdict(set)
    for rec in graph.run(
        """UNWIND $names AS gene
           MATCH (re:ReferenceEntity)<-[:referenceEntity]-(pe:PhysicalEntity)
           WHERE gene IN re.geneName
           RETURN gene AS gene, COLLECT(DISTINCT toString(pe.dbId)) AS dbids""",
        names=genes,
    ).data():
        gene_dbids[rec["gene"]] = set(rec["dbids"])

    verdicts = []
    nets: dict = {}
    for r in wrong:
        pw = r["pathway"]
        if pw not in nets:
            nets[pw] = load_mpbio_network(pw)
        adj = nets[pw]
        if adj is None:
            verdicts.append((r, "no_2019_network"))
            continue
        srcs = {d for d in gene_dbids.get(r["gene"], set()) if d in adj}
        target = r["key_output"].strip()
        if not srcs:
            verdicts.append((r, "gene_absent_from_2019_network"))
            continue
        then = target in reachable(adj, srcs)
        now = r["category"] != "no_path"
        if then and not now:
            verdicts.append((r, "lost_since_2019"))
        elif now and not then:
            verdicts.append((r, "gained_since_2019"))
        else:
            verdicts.append((r, "same_structure"))

    counts = collections.Counter(v for _, v in verdicts)
    total = len(verdicts)
    print(f"failing cases examined: {total}\n")
    for k, n in counts.most_common():
        print(f"   {k:<32} {n:>6}  ({100 * n / total:.1f}%)")

    # the two buckets where the ground truth may simply be out of date
    stale = counts["lost_since_2019"] + counts["gained_since_2019"]
    print(f"\n   structure CHANGED under the ground truth: {stale} "
          f"({100 * stale / total:.1f}% of failures)")

    if args.out:
        with args.out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["pathway", "gene", "direction", "key_output",
                        "predicted", "expected", "category", "skew"])
            for r, v in verdicts:
                w.writerow([r["pathway"], r["gene"], r["direction"], r["key_output"],
                            r["predicted"], r["expected"], r["category"], v])
        print(f"\nper-case: {args.out}")


if __name__ == "__main__":
    main()
