#!/usr/bin/env python3
"""Check that DeltaSignal's cofactor list still covers the generator's.

Two repositories hold a list of the same chemistry for two different purposes.
The generator's `_COFACTOR_STIDS` stops the diagram-bridge pass asserting that
a producer of ATP feeds a consumer of ATP; DeltaSignal's `COFACTOR_STIDS`
stops a perturbation travelling through one during a solve. They are not the
same decision, but the generator's set must be a SUBSET of the solver's: a
molecule shared enough that bridging on it is meaningless is certainly shared
enough that conducting through it is meaningless.

They diverged once already — the solver's list was derived from a Neo4j query
by molecule name and silently omitted NAD+ entirely, plus three compartment
variants the generator's hand-written list had. Nothing caught it, because
each list is only read by its own repo.

This is a cross-repo check, so it lives here rather than in either test suite:
it needs both checkouts present and must not fail CI when only one is.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DEFAULT_GENERATOR = Path.home() / "gitroot" / "logic-network-generator"


def extract(text: str, pattern: str) -> set[str]:
    match = re.search(pattern, text, re.S)
    if not match:
        raise SystemExit(f"could not locate the list matching {pattern!r}")
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
                       r"_COFACTOR_STIDS: frozenset = frozenset\(\{(.*?)\}\)")
    downstream = extract(solver.read_text(),
                         r"const COFACTOR_STIDS = Set\(\[(.*?)\]\)")

    missing = sorted(upstream - downstream)
    print(f"generator list: {len(upstream)}   solver list: {len(downstream)}")
    if missing:
        print(f"MISSING from the solver's list ({len(missing)}):")
        for stid in missing:
            print(f"  {stid}")
        return 1
    print("solver list covers the generator's list")
    return 0


if __name__ == "__main__":
    sys.exit(main())
