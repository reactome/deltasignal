#!/usr/bin/env python3
"""Report the difference between the two repositories' cofactor lists.

Two repositories hold a list of the same chemistry for two different questions.
The generator's `_COFACTOR_STIDS` decides whether a diagram bridge may be drawn
across a node; DeltaSignal's `COFACTOR_STIDS` decides whether a perturbation may
travel through one during a solve. Neither has to contain the other.

An earlier version of this script asserted that the solver's list must be a
superset, and that was wrong: it would import four defective entries. Six of the
generator's thirteen ids are stale or mislabelled against Release97, so a plain
set difference reads as "the solver is missing things" when the opposite is
true. Every known difference is therefore annotated below with what the id
actually is, and only an UNEXPLAINED difference is an error — that is the case
worth waking someone for, because it means the lists have drifted again.

The solver's list is derived from the release by ChEBI identity, so it is the
one that self-corrects. The generator's is still hand-written; correcting it
changes which bridges are suppressed, which changes the networks, which needs
its own A/B, so it is deliberately left alone and explained here instead.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DEFAULT_GENERATOR = Path.home() / "gitroot" / "logic-network-generator"

# Generator ids the solver deliberately does NOT carry, with what Release97
# says each one actually is. Verified against Neo4j on 2026-09-14.
KNOWN_DIFFERENCES = {
    "R-ALL-29390": "commented 'Pi variant'; is PXLP (pyridoxal 5'-phosphate)",
    "R-ALL-217093": "commented 'NADP+'; absent from Release97",
    "R-ALL-110114": "commented 'NADPH'; absent from Release97",
    "R-ALL-29986": "commented 'NAD+'; absent from Release97",
}


def extract(text: str, pattern: str, what: str) -> set[str]:
    match = re.search(pattern, text, re.S)
    if not match:
        raise SystemExit(f"could not locate the {what} list")
    return set(re.findall(r'"(R-[A-Z]+-\d+)"', match.group(1)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--generator-root", type=Path, default=DEFAULT_GENERATOR)
    ap.add_argument("--deltasignal-root", type=Path,
                    default=Path(__file__).resolve().parents[2])
    args = ap.parse_args()

    generator = args.generator_root / "src" / "logic_network_generator.py"
    solver = args.deltasignal_root / "src" / "core" / "cofactors.jl"
    if not generator.exists():
        print(f"generator checkout not found at {args.generator_root}; skipping")
        return 0

    upstream = extract(generator.read_text(),
                       r"_COFACTOR_STIDS: frozenset = frozenset\(\{(.*?)\}\)",
                       "generator")
    downstream = extract(solver.read_text(),
                         r"const COFACTOR_STIDS = Set\(\[(.*?)\]\)",
                         "solver")

    print(f"generator list: {len(upstream)}   solver list: {len(downstream)}")

    only_upstream = sorted(upstream - downstream)
    explained = [s for s in only_upstream if s in KNOWN_DIFFERENCES]
    unexplained = [s for s in only_upstream if s not in KNOWN_DIFFERENCES]

    if explained:
        print(f"\nin the generator's list only, and deliberately so ({len(explained)}):")
        for stid in explained:
            print(f"  {stid}  — {KNOWN_DIFFERENCES[stid]}")

    if unexplained:
        print(f"\nUNEXPLAINED — in the generator's list, not in the solver's "
              f"({len(unexplained)}):")
        for stid in unexplained:
            print(f"  {stid}")
        print("\nThe lists have drifted. Either the solver's derivation no longer "
              "covers a molecule the generator treats as a cofactor, or the "
              "generator has gained a new entry. Check which before changing "
              "either list.")
        return 1

    print("\nno unexplained differences")
    return 0


if __name__ == "__main__":
    sys.exit(main())
