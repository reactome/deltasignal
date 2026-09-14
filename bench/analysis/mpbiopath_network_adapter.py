#!/usr/bin/env python3
"""Convert MP-BioPath's published networks into DeltaSignal catalog directories.

The point is a control that has never been run: DeltaSignal's propagator on
MP-BioPath's OWN networks, scoring the same cases. DeltaSignal scores 365 of
564 on its networks and MP-BioPath scores 407 on its own, but MP-BioPath's
published pathways had loops removed by hand — four of the nine benchmark
networks are fully acyclic and the largest component anywhere is eleven
nodes, against 836 in ours. Those two numbers have never described the same
task. Running our propagator on their networks splits the 42-case deficit
into a propagator share and a network share.

Format (no header, four columns), confirmed against MP-BioPath's own reader
in mp-biopath/src/pi.jl:

    parent dbid | child dbid | polarity: 1 = POS, -1 = NEG | 0 = AND, 1 = OR

The conjunction convention is the one thing here that could silently invert
the model, so it was verified rather than assumed: pi.jl documents "if there
is only one parent it will be an AND relation", and of children with exactly
one parent 5,148 carry flag 0 against 47 carrying flag 1 (99.1%).

Node identity needs no translation. MP-BioPath's node ids are Reactome
database identifiers, which is the same key the benchmark cases already use
(211 of 212 ERBB2 node ids resolve in db_id_to_name_mapping.txt). Emitting
stable_id as "R-HSA-<dbid>" makes load_dbid_to_uuids' stable_id.rsplit("-",
1)[-1] recover exactly that dbid.

The uuid column carries an "mpb-" prefix rather than the bare dbid. Bare
dbids are all-numeric, CSV.jl types that column as Int64, and the parser's
String(row.source_id) then throws a MethodError that surfaces as an opaque
HTTP 400 (src/io/tsv_parser.jl:98). LNG uuids are never all-numeric so this
had never been exercised. The prefix is on the uuid only; stable_id keeps
the bare dbid, so the benchmark's mapping is unaffected.

Only logic_network.csv and stid_to_uuid_mapping.csv are written: those are
the server's stated minimum (src/api/server.jl:423-425). nodes.csv and
entity_reaction_proxy_mapping.csv are optional fallbacks and are deliberately
omitted rather than fabricated, since a fabricated proxy map would let a case
score through a mapping MP-BioPath never had.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

UUID_PREFIX = "mpb-"
POLARITY = {"1": "pos", "-1": "neg"}
CONJUNCTION = {"0": "and", "1": "or"}


def parse_network(path: Path) -> list[tuple[str, str, str, str]]:
    edges: list[tuple[str, str, str, str]] = []
    with path.open() as handle:
        for lineno, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"{path.name}:{lineno}: expected 4 columns, got {len(parts)}"
                )
            parent, child, polarity, conjunction = (p.strip() for p in parts[:4])
            if polarity not in POLARITY:
                raise ValueError(f"{path.name}:{lineno}: unknown polarity {polarity!r}")
            if conjunction not in CONJUNCTION:
                raise ValueError(
                    f"{path.name}:{lineno}: unknown conjunction flag {conjunction!r}"
                )
            edges.append((parent, child, POLARITY[polarity], CONJUNCTION[conjunction]))
    return edges


def catalog_ids(catalog: Path) -> dict[str, str]:
    """Map each catalog directory's leading name to its R-HSA suffix."""
    mapping: dict[str, str] = {}
    for directory in sorted(p for p in catalog.iterdir() if p.is_dir()):
        name = directory.name
        if "_R-HSA-" not in name:
            continue
        prefix, _, numeric = name.rpartition("_R-HSA-")
        mapping[prefix] = f"R-HSA-{numeric}"
    return mapping


def case_pathway_names(cases_tsv: Path) -> dict[str, str]:
    """Map R-HSA id -> the benchmark's own pathway_name.

    The LNG directory prefix is NOT a reliable key for finding the
    MP-BioPath file: R-HSA-453279 is "Mitotic_G1_phase_and_G1_S_transition"
    in the catalog and "Mitotic_G1-G1_S_phases" in MP-BioPath, and
    R-HSA-5693567's file carries a trailing underscore the directory drops.
    The benchmark's pathway_name column matches MP-BioPath's filenames for
    all ten, so use it as the bridge rather than guessing from directory
    names, which silently lost two pathways on the first attempt.
    """
    mapping: dict[str, str] = {}
    with cases_tsv.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pid = (row.get("pathway_id") or "").strip()
            name = (row.get("pathway_name") or "").strip()
            if pid and name:
                mapping[pid] = name
    return mapping


def write_catalog_dir(out: Path, edges: list[tuple[str, str, str, str]]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "logic_network.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_id", "target_id", "pos_neg", "and_or",
                         "edge_type", "stoichiometry", "edge_reaction_id"])
        for parent, child, sign, conj in edges:
            parent, child = f"{UUID_PREFIX}{parent}", f"{UUID_PREFIX}{child}"
            # edge_type is left empty rather than guessed. The solver routes
            # "assembly" to a min and treats "depletion" specially, and
            # inventing either would be modelling MP-BioPath's networks as
            # something they do not say they are.
            writer.writerow([parent, child, sign, conj, "", 1, ""])

    nodes = sorted({e[0] for e in edges} | {e[1] for e in edges})
    with (out / "stid_to_uuid_mapping.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["uuid", "stable_id"])
        for node in nodes:
            writer.writerow([f"{UUID_PREFIX}{node}", f"R-HSA-{node}"])


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pathways", type=Path, required=True,
                        help="Directory of MP-BioPath four-column TSVs")
    parser.add_argument("--catalog", type=Path, required=True,
                        help="LNG catalog, read only, to source the R-HSA ids")
    parser.add_argument("--out", type=Path, required=True,
                        help="Directory to write the adapted catalog into")
    parser.add_argument("--cases", type=Path, required=True,
                        help="A benchmark_cases.tsv, used to map R-HSA ids to "
                             "the pathway names MP-BioPath's files carry")
    args = parser.parse_args()

    ids = catalog_ids(args.catalog)
    if not ids:
        print(f"No R-HSA-suffixed directories under {args.catalog}", file=sys.stderr)
        return 1

    names = case_pathway_names(args.cases)

    args.out.mkdir(parents=True, exist_ok=True)
    written, missing_mpb = [], []
    for prefix, stable in sorted(ids.items()):
        # Try the benchmark's own name first, then the directory prefix.
        candidates = [names.get(stable), prefix]
        source = next(
            (args.pathways / f"{c}.tsv" for c in candidates
             if c and (args.pathways / f"{c}.tsv").exists()),
            None,
        )
        if source is None:
            missing_mpb.append(f"{prefix} ({stable})")
            continue
        edges = parse_network(source)
        target = args.out / f"{prefix}_{stable}"
        write_catalog_dir(target, edges)
        written.append((target.name, len(edges)))

    mpb_files = {p.stem for p in args.pathways.glob("*.tsv")}
    unmatched = sorted(mpb_files - set(names.values()) - set(ids))

    for name, n in written:
        print(f"wrote {name:56s} {n:6d} edges")
    print(f"\n{len(written)} of {len(ids)} catalog pathways adapted")
    if missing_mpb:
        # Reported, not silently dropped: a pathway absent here is absent
        # from the control and must be excluded from the LNG arm too.
        print(f"NO MP-BioPath network for: {', '.join(missing_mpb)}")
    print(f"MP-BioPath pathways with no catalog counterpart: {len(unmatched)}"
          " (not part of this benchmark)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
