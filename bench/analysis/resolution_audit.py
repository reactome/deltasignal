#!/usr/bin/env python3
"""Audit the entity-to-node mapping in both directions.

Feature 005 requires that every Reactome entity resolve to nodes and every
node resolve back, with anything unresolved declared rather than absent.
This reports where a catalog stands against that, and is the measurement
behind LNG #67's before/after.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def orphan_report(catalog: Path) -> dict:
    """Context rows naming a node that is not in the network. #67's metric."""
    total: Counter = Counter()
    orphaned: Counter = Counter()
    per_pathway: dict[str, tuple[int, int]] = {}
    for directory in sorted(p for p in catalog.iterdir() if p.is_dir()):
        net = directory / "logic_network.csv"
        ctx = directory / "node_reaction_context.csv"
        if not (net.exists() and ctx.exists()):
            continue
        live: set[str] = set()
        with net.open(newline="") as handle:
            for row in csv.DictReader(handle):
                live.add((row["source_id"] or "").strip())
                live.add((row["target_id"] or "").strip())
        n = o = 0
        with ctx.open(newline="") as handle:
            for row in csv.DictReader(handle):
                role = (row.get("role") or "").strip()
                total[role] += 1
                n += 1
                if (row.get("context_node") or "").strip() not in live:
                    orphaned[role] += 1
                    o += 1
        per_pathway[directory.name] = (n, o)
    return {"total": total, "orphaned": orphaned, "per_pathway": per_pathway}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--compare", type=Path,
                        help="A second catalog to diff against (the before)")
    args = parser.parse_args()

    def show(label: str, rep: dict) -> tuple[int, int]:
        print(f"\n# {label}")
        print(f"{'role':12s} {'rows':>8} {'orphaned':>9} {'%':>7}")
        tot = orph = 0
        for role in sorted(rep["total"]):
            t, o = rep["total"][role], rep["orphaned"][role]
            tot += t
            orph += o
            print(f"{role:12s} {t:8d} {o:9d} {100*o/t if t else 0:6.1f}%")
        print(f"{'TOTAL':12s} {tot:8d} {orph:9d} {100*orph/tot if tot else 0:6.1f}%")
        bad = [k for k, (n, o) in rep["per_pathway"].items() if o]
        if bad:
            print(f"  pathways with orphans: {len(bad)} — {', '.join(sorted(bad)[:3])}...")
        return tot, orph

    after = orphan_report(args.catalog)
    if args.compare:
        show("BEFORE", orphan_report(args.compare))
    tot, orph = show("AFTER" if args.compare else str(args.catalog), after)

    if orph:
        print(f"\nFAIL: {orph} orphaned context rows remain")
        return 1
    print("\nOK: every context row names a node present in its network")
    return 0


if __name__ == "__main__":
    sys.exit(main())
