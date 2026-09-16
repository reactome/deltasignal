#!/usr/bin/env python3
"""Split the failing cases into "ours to fix" and "not fixable by us".

An accuracy number says how many cases are wrong. It does not say how many
COULD be right, and those are different questions with different work behind
them. This asks, per failing case, whether the failure is:

  network_missing_path   the logic network has no route but REACTOME DOES.
                         The generator lost a connection curators made. Fixable
                         upstream, and a bug rather than a modelling choice.

  truly_unreachable      neither the logic network nor Reactome has a directed
                         route. The curator's call rests on something the
                         pathway graph does not encode — co-regulation, a
                         cross-pathway effect, or knowledge outside the
                         diagram. Not fixable by connectivity in either repo.

  propagator_*           a route exists and the propagator got it wrong. Ours.

The Reactome side walks the same bipartite structure the logic network encodes
(entity feeds a reaction; the reaction produces entities), restricted to events
inside the pathway, so the two are asked an equivalent question.
"""

from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from curated_reachability import build_pathway_graph, reachable  # noqa: E402


def load_cases(path: Path):
    with path.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r.get("valid") == "1"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--limit-pathways", type=int, default=0)
    ap.add_argument("--containment", choices=("none", "down", "both"),
                    default="both",
                    help="how permissive the Reactome-side walk is; "
                         "'both' also admits co-membership routes")
    args = ap.parse_args()

    import benchmark_vs_mpbiopath as B

    cases = load_cases(args.cases)
    wrong = [r for r in cases if r["predicted"] != r["expected"]]
    by_pathway = collections.defaultdict(list)
    for r in wrong:
        by_pathway[r["pathway"]].append(r)

    # Resolve every gene and every readout once.
    genes = sorted({r["gene"] for r in wrong})
    gene_stids = B.neo4j_gene_to_stids(genes)
    dbids = sorted({r["key_output"] for r in wrong})
    dbid_stid = B.neo4j_dbid_to_stid(dbids)

    dirs = {d.name.rsplit("_R-HSA-", 1)[0]: d
            for d in args.catalog.iterdir() if d.is_dir()}
    stid_of_dir = {d.name.rsplit("_R-HSA-", 1)[0]: "R-HSA-" + d.name.rsplit("_R-HSA-", 1)[1]
                   for d in args.catalog.iterdir() if d.is_dir()}

    verdicts = []
    names = sorted(by_pathway, key=lambda n: -len(by_pathway[n]))
    if args.limit_pathways:
        names = names[:args.limit_pathways]

    for i, name in enumerate(names, 1):
        rows = by_pathway[name]
        if name not in dirs:
            for r in rows:
                verdicts.append((r, "pathway_dir_missing"))
            continue
        pid = stid_of_dir[name]
        try:
            feeds, produces, contains, entities = build_pathway_graph(pid)
        except Exception as exc:                       # pragma: no cover
            print(f"  [skip] {name}: {exc}", flush=True)
            for r in rows:
                verdicts.append((r, "curated_graph_error"))
            continue

        cache: dict[tuple, set] = {}
        for r in rows:
            cat = r.get("category")
            if cat != "no_path":
                # A route exists in the logic network; the propagator owns it.
                verdicts.append((r, "propagator_" + (cat or "unknown")))
                continue
            srcs = {s for s in gene_stids.get(r["gene"], []) if s in entities}
            tgt = dbid_stid.get(r["key_output"])
            if not srcs:
                verdicts.append((r, "gene_absent_from_curated_pathway"))
                continue
            if not tgt or tgt not in entities:
                verdicts.append((r, "readout_absent_from_curated_pathway"))
                continue
            key = tuple(sorted(srcs))
            if key not in cache:
                cache[key] = reachable(feeds, produces, contains, set(srcs),
                                       containment=args.containment)
            verdicts.append((r, "network_missing_path" if tgt in cache[key]
                             else "truly_unreachable"))
        print(f"  [{i}/{len(names)}] {name[:46]:<48} {len(rows):>5} wrong", flush=True)

    counts = collections.Counter(v for _, v in verdicts)
    total = len(verdicts)
    print(f"\ntriaged {total} failing cases\n")
    ours = sum(n for k, n in counts.items() if k.startswith("propagator_"))
    upstream = counts["network_missing_path"]
    for k, n in counts.most_common():
        print(f"   {k:<42} {n:>6}  ({100 * n / total:.1f}%)")
    print(f"\n   FIXABLE IN DELTASIGNAL (a route exists):      {ours:>6}  "
          f"({100 * ours / total:.1f}% of failures)")
    print(f"   FIXABLE UPSTREAM (Reactome has the route):    {upstream:>6}  "
          f"({100 * upstream / total:.1f}% of failures)")

    if args.out:
        with args.out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["pathway", "gene", "direction", "key_output",
                        "predicted", "expected", "verdict"])
            for r, v in verdicts:
                w.writerow([r["pathway"], r["gene"], r["direction"],
                            r["key_output"], r["predicted"], r["expected"], v])
        print(f"\nper-case verdicts: {args.out}")


if __name__ == "__main__":
    main()
