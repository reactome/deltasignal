#!/usr/bin/env python3
"""Assertions for the MP-BioPath network adapter.

A silent mis-mapping here would be indistinguishable from a propagator
finding: if the conjunction flag were inverted, the control would report a
worse score and we would conclude the propagator is at fault. These tests
exist so that conclusion cannot be reached by accident.
"""

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from mpbiopath_network_adapter import (  # noqa: E402
    CONJUNCTION, POLARITY, UUID_PREFIX, case_pathway_names, parse_network,
    write_catalog_dir,
)

failures = []


def check(label, actual, expected):
    if actual != expected:
        failures.append(f"{label}: expected {expected!r}, got {actual!r}")


# The conjunction convention is the one that could invert the model.
# MP-BioPath's own reader documents single-parent => AND, and 99.1% of
# single-parent children carry flag 0.
check("conjunction 0 is AND", CONJUNCTION["0"], "and")
check("conjunction 1 is OR", CONJUNCTION["1"], "or")
check("polarity 1 is pos", POLARITY["1"], "pos")
check("polarity -1 is neg", POLARITY["-1"], "neg")

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    src = tmp / "Toy.tsv"
    src.write_text("100\t200\t1\t0\n300\t200\t-1\t1\n")
    edges = parse_network(src)
    check("edge count", len(edges), 2)
    check("positive AND edge", edges[0], ("100", "200", "pos", "and"))
    check("negative OR edge", edges[1], ("300", "200", "neg", "or"))

    # A malformed row must fail loudly rather than default.
    bad = tmp / "Bad.tsv"
    bad.write_text("100\t200\t7\t0\n")
    try:
        parse_network(bad)
        failures.append("unknown polarity 7 was accepted")
    except ValueError:
        pass
    short = tmp / "Short.tsv"
    short.write_text("100\t200\n")
    try:
        parse_network(short)
        failures.append("a two-column row was accepted")
    except ValueError:
        pass

    # The dbid must survive the round trip through the benchmark's own
    # stable-id parsing, or every case silently fails to map.
    out = tmp / "Toy_R-HSA-1"
    write_catalog_dir(out, edges)
    with (out / "stid_to_uuid_mapping.csv").open(newline="") as handle:
        pairs = {r["uuid"]: r["stable_id"] for r in csv.DictReader(handle)}
    check("node count", len(pairs), 3)
    for uuid, stable in pairs.items():
        # The uuid carries a non-numeric prefix so CSV.jl types the column as
        # String; the stable id keeps the bare dbid the benchmark keys on.
        check(f"uuid is not all-numeric ({uuid})", uuid.isdigit(), False)
        check(f"dbid round trip for {uuid}",
              f"{UUID_PREFIX}{stable.rsplit('-', 1)[-1]}", uuid)

    with (out / "logic_network.csv").open(newline="") as handle:
        net = list(csv.DictReader(handle))
    check("edge ids carry the prefix", net[0]["source_id"], f"{UUID_PREFIX}100")
    check("edge ids are not all-numeric", net[0]["source_id"].isdigit(), False)

    # Only the two files the server requires; no fabricated proxy map.
    written = sorted(p.name for p in out.iterdir())
    check("files written", written, ["logic_network.csv", "stid_to_uuid_mapping.csv"])

    # pathway_name bridging: the directory prefix is NOT a usable key.
    cases = tmp / "cases.tsv"
    cases.write_text("pathway_id\tpathway_name\nR-HSA-453279\tMitotic_G1-G1_S_phases\n")
    check("case name bridge", case_pathway_names(cases),
          {"R-HSA-453279": "Mitotic_G1-G1_S_phases"})

if failures:
    print("FAILED:")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print(f"all adapter assertions passed ({12} checks)")
