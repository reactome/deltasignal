#!/usr/bin/env python3
"""Does the no-path wall hold beyond the ten benchmark pathways?

The decomposition on the ten-pathway set says 18.7% of cases have no directed
path from perturbation to readout, and that on those we score ~0.22 because
the model correctly answers "no change" while the ground truth says the
readout moved. Two pathways supply 60% of that case set, so the finding needs
checking at scale before anything is built on it.

No solver is needed. With no path the prediction is NO_CHANGE by
construction, so accuracy on those cases is exactly the fraction whose ground
truth is also "no change" — computable from the graph and the truth table
alone. That makes an 85-pathway check cheap.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

DOWN, NORMAL, UP = "0", "1", "2"


def load_network(pathway_dir: Path):
    adjacency: dict[str, list[str]] = defaultdict(list)
    sources: set[str] = set()
    targets: set[str] = set()
    with (pathway_dir / "logic_network.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            s, t = row["source_id"].strip(), row["target_id"].strip()
            adjacency[s].append(t)
            sources.add(s)
            targets.add(t)
    return adjacency, sources - targets


def load_dbid_map(pathway_dir: Path) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = defaultdict(list)
    nodes = pathway_dir / "nodes.csv"
    if nodes.exists():
        with nodes.open(newline="") as handle:
            for row in csv.DictReader(handle):
                ids = {(row.get("diagram_entity_id") or "").strip()}
                ids |= {x.strip() for x in (row.get("member_leaves") or "").split("|")}
                for stable_id in ids - {""}:
                    mapping[stable_id.rsplit("-", 1)[-1]].append(row["uuid"])
    resolution = pathway_dir / "node_resolution.csv"
    if resolution.exists():
        with resolution.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["relation"] == "set_member":
                    dbid = row["stable_id"].rsplit("-", 1)[-1]
                    if row["uuid"] not in mapping[dbid]:
                        mapping[dbid].append(row["uuid"])
    return mapping


def reachable(adjacency, sources, targets) -> bool:
    target_set = set(targets)
    seen = set(sources)
    queue = deque(sources)
    while queue:
        node = queue.popleft()
        if node in target_set:
            return True
        for nxt in adjacency.get(node, ()):
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return bool(seen & target_set)


def load_curator(path: Path):
    """key_output rows against GENE_0 / GENE_2 columns, values 0/1/2."""
    aliases = ("key_output", "key output", "key_outout", "key outout", "key_ouput")
    lines = [ln.rstrip("\n").split("\t") for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0]
    key_col = next((header.index(a) for a in aliases if a in header), None)
    if key_col is None:
        return []
    cases = []
    for row in lines[1:]:
        if len(row) <= key_col or not row[key_col].strip():
            continue
        key_output = row[key_col].strip()
        for i, column in enumerate(header):
            if i == key_col or column.strip().lower() == "control":
                continue
            if "_" not in column:
                continue
            gene, _, direction = column.rpartition("_")
            if direction not in ("0", "2") or i >= len(row):
                continue
            value = row[i].strip()
            if value in (DOWN, NORMAL, UP):
                cases.append((key_output, gene, direction, value))
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--curator-dir", type=Path, required=True)
    parser.add_argument("--id-map", type=Path, required=True)
    parser.add_argument("--per-pathway", action="store_true")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmark_mpbiopath_cases import load_gene_dbids

    dirs = {p.name.rsplit("_R-HSA-", 1)[0]: p for p in args.catalog.iterdir() if p.is_dir()}
    totals = Counter()
    truth_when_nopath = Counter()
    truth_when_path = Counter()
    per_pathway = []

    for name, directory in sorted(dirs.items()):
        curator = args.curator_dir / f"{name}_reactome_curator_results.tsv"
        if not curator.exists() or not (directory / "logic_network.csv").exists():
            continue
        cases = load_curator(curator)
        if not cases:
            continue
        adjacency, roots = load_network(directory)
        dbid_map = load_dbid_map(directory)
        gene_dbids = load_gene_dbids(args.id_map, {c[1] for c in cases})
        n = np_ = 0
        for key_output, gene, _direction, truth in cases:
            out_uuids = dbid_map.get(key_output, [])
            gene_uuids = sorted({u for db in gene_dbids.get(gene, ())
                                 for u in dbid_map.get(db, ())} & roots)
            if not (out_uuids and gene_uuids):
                totals["unmappable"] += 1
                continue
            n += 1
            totals["scoreable"] += 1
            if reachable(adjacency, gene_uuids, out_uuids):
                totals["path"] += 1
                truth_when_path[truth] += 1
            else:
                np_ += 1
                totals["nopath"] += 1
                truth_when_nopath[truth] += 1
        if n:
            per_pathway.append((name, n, np_ / n))

    scoreable = totals["scoreable"]
    if not scoreable:
        print("no scoreable cases — is the catalog generated?", file=sys.stderr)
        return 1
    nopath = totals["nopath"]
    print(f"pathways with both a network and curator truth: {len(per_pathway)}")
    print(f"cases: {scoreable} scoreable, {totals['unmappable']} unmappable\n")
    print(f"  no directed path : {nopath:6d}  ({100*nopath/scoreable:5.1f}%)")
    print(f"  path exists      : {totals['path']:6d}  ({100*totals['path']/scoreable:5.1f}%)")

    # With no path the model answers NO_CHANGE, so this IS its accuracy there.
    nc = truth_when_nopath[NORMAL]
    print(f"\n  ground truth on the no-path cases: "
          f"DOWN {truth_when_nopath[DOWN]}, NO_CHANGE {nc}, UP {truth_when_nopath[UP]}")
    print(f"  => deterministic accuracy ceiling on no-path cases: "
          f"{nc}/{nopath} = {nc/nopath:.3f}" if nopath else "")
    pnc = truth_when_path[NORMAL]
    print(f"  ground truth where a path exists:  "
          f"DOWN {truth_when_path[DOWN]}, NO_CHANGE {pnc}, UP {truth_when_path[UP]}")

    if args.per_pathway:
        per_pathway.sort(key=lambda t: -t[2])
        print(f"\n  {'pathway':52s} {'n':>6} {'no-path%':>9}")
        for name, n, rate in per_pathway:
            print(f"  {name[:52]:52s} {n:6d} {100*rate:8.1f}%")
        rates = sorted(r for _, _, r in per_pathway)
        mid = rates[len(rates)//2]
        print(f"\n  median pathway no-path rate: {100*mid:.1f}%  "
              f"(range {100*rates[0]:.1f}% - {100*rates[-1]:.1f}%)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
