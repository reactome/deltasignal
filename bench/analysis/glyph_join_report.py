#!/usr/bin/env python3
"""How completely do generated nodes join to diagram glyphs, and both ways?

Adam's requirement: clicking one of the two ATP glyphs in the pathway browser
must resolve to what THAT glyph means, and a solved node must resolve back to
the glyph to highlight.

Three populations, and all three are reported because only the difference
between them distinguishes an intended absence from a defect:

  joinable    a diagram triple with a generated node
  diagram-only  drawn but not generated
  LNG-only      generated but not drawn — EXPECTED, because generation
                descends below the diagram's own reactions
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def load_resolution(pathway_dir: Path) -> list[dict]:
    path = pathway_dir / "node_resolution.csv"
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--diagrams", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(Path.home() / "gitroot" / "logic-network-generator"))
    import os
    os.environ.setdefault("LNG_DIAGRAM_DIR", str(args.diagrams))
    from src.diagram_connectivity import diagram_glyph_positions

    totals = Counter()
    print(f"{'pathway':34s} {'joinable':>9} {'diagram-only':>13} {'LNG-only':>9} {'cover':>7}")
    for directory in sorted(p for p in args.catalog.iterdir() if p.is_dir()):
        pathway_id = "R-HSA-" + directory.name.rsplit("_R-HSA-", 1)[-1]
        rows = load_resolution(directory)
        if not rows:
            continue
        try:
            diagram = diagram_glyph_positions(pathway_id)
        except Exception:
            continue
        # LNG's own (reaction, entity, role) triples, from the rows that name
        # a reaction position at all.
        lng = {(r["reaction_stid"], r["stable_id"], r["role"])
               for r in rows if r["reaction_stid"] and r["role"]}
        joinable = lng & set(diagram)
        diagram_only = set(diagram) - lng
        lng_only = lng - set(diagram)
        totals["joinable"] += len(joinable)
        totals["diagram_only"] += len(diagram_only)
        totals["lng_only"] += len(lng_only)
        cover = len(joinable) / len(diagram) if diagram else 0.0
        print(f"{directory.name.split('_R-HSA-')[0][:34]:34s} {len(joinable):9d} "
              f"{len(diagram_only):13d} {len(lng_only):9d} {cover:6.1%}")

    diagram_total = totals["joinable"] + totals["diagram_only"]
    print(f"\n  joinable {totals['joinable']} of {diagram_total} diagram triples "
          f"({totals['joinable']/max(1,diagram_total):.1%})")
    print(f"  LNG-only {totals['lng_only']} — EXPECTED: generation descends below the "
          f"diagram's own reactions,\n    so these are not failures and are reported "
          f"so that 'expected' is checkable rather than assumed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
