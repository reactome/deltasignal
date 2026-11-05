#!/usr/bin/env python3
"""specs/134 amendment 3: which total readouts' abundance copies differ between two versions of the rule, per network.

Structural only (no measurement read, no solve): for each pathway of a build, each total readout's nodes as the korkut
or CPPA bench resolves them, cut by `abundance_copies` under versions A and B. Prints one row per (pathway, readout)
whose copies differ, then counts. Run with the generator venv python (py2neo).

Usage: copies_diff.py BUILD_DIR READOUTS_TSV [VERSION_A VERSION_B]   (default abundance_a2 abundance_a3)
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import benchmark_vs_mpbiopath as bm  # noqa: E402
import korkut_bench as kb  # noqa: E402


def main():
    build, readouts_tsv = Path(sys.argv[1]), sys.argv[2]
    va, vb = (sys.argv[3:5] if len(sys.argv) > 4 else ("abundance_a2", "abundance_a3"))
    readouts = [r for r in kb.read_tsv(readouts_tsv) if r["kind"] == "total"]
    genes = sorted({g for r in readouts for g in r["sites"].split(",")})
    gene_to_stids = kb.neo4j_protein_stids(bm.neo4j_gene_to_stids(genes))
    r_stids = kb.readout_stids(readouts, gene_to_stids, {})
    counts = collections.Counter()
    print("pathway\treadout\tcopies_" + va + "\tcopies_" + vb)
    for pdir in sorted(p for p in build.iterdir() if (p / "stid_to_uuid_mapping.csv").exists()):
        ro = kb.readout_nodes(r_stids, bm.load_stid_to_uuids(pdir), bm.context_own_uuids(pdir))
        if not ro:
            continue
        with bm.open_network(pdir) as f:
            rows = list(csv.DictReader(f))
        for name, nodes in sorted(ro.items()):
            a, b = kb.abundance_copies(nodes, rows, va), kb.abundance_copies(nodes, rows, vb)
            counts["pairs"] += 1
            if a != b:
                counts["differ"] += 1
                counts["empty_to_some" if not a else "some_to_empty" if not b else "both_some"] += 1
                print(f"{pdir.name}\t{name}\t{len(a)}\t{len(b)}")
    print("#", dict(counts), file=sys.stderr)


if __name__ == "__main__":
    main()
