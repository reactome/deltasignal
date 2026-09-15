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


def compare_molecules(generator_src: str, solver_src: str) -> int:
    """Compare the molecule sets behind the two derived lists.

    The generator names molecules as `_COFACTOR_CHEBI` keys; the solver groups
    its baked stable ids under `# Molecule` comment headers. If the generator
    gains or loses a molecule, every freshly generated bundle starts declaring
    a different set, and because the bundle overrides the solver's built-in
    list, the solver's own copy silently becomes the stale fallback for older
    bundles only. That is the drift this feature introduced.
    """
    gen = re.search(r"_COFACTOR_CHEBI: Dict\[str, List\[str\]\] = \{(.*?)\n\}",
                    generator_src, re.S)
    if not gen:
        print("\ngenerator has no _COFACTOR_CHEBI (pre-dates the derived list); "
              "skipping the molecule audit")
        return 0
    gen_molecules = set(re.findall(r'"([^"]+)":\s*\[', gen.group(1)))

    block = re.search(r"const COFACTOR_STIDS = Set\(\[(.*?)\]\)", solver_src, re.S)
    sol_molecules = set(re.findall(r"^\s*#\s*(.+?)\s*$", block.group(1), re.M))

    print(f"\ngenerator ChEBI molecules: {len(gen_molecules)}   "
          f"solver molecule groups: {len(sol_molecules)}")
    only_gen = sorted(gen_molecules - sol_molecules)
    only_sol = sorted(sol_molecules - gen_molecules)
    if not only_gen and not only_sol:
        print("  the two derived lists cover the same molecules")
        return 0
    for name in only_gen:
        print(f"  ONLY the generator derives: {name}  "
              f"(new bundles will declare it; the solver's fallback will not)")
    for name in only_sol:
        print(f"  ONLY the solver carries: {name}  "
              f"(new bundles will NOT declare it, so it stops being pinned)")
    return 1


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

    print(f"generator bridge list: {len(upstream)}   solver list: {len(downstream)}")

    # The list that can actually change a solve is the generator's ChEBI set,
    # because it is what `get_cofactor_species` derives and `cofactors.csv`
    # ships, and a bundled list OVERRIDES the solver's own. Comparing only the
    # 13-id bridge list audits the pair that no longer matters most.
    chebi_src = args.generator_root / "src" / "neo4j_connector.py"
    rc = compare_molecules(
        chebi_src.read_text() if chebi_src.exists() else "",
        solver.read_text())

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
    if rc:
        return rc

    print("\nno unexplained differences in the bridge list")
    return rc


if __name__ == "__main__":
    sys.exit(main())
