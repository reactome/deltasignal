#!/usr/bin/env python3
"""specs/142 parity gate, the comparison: any flag difference fails (exit 1).

Compares bench_rails.json (the bench's railed(), imported from score_124.py)
against engine_rails.json (the engine's served readout_railed) case by case on
identical folds. Also checks the engine's stops are the bench's 0.01 / 99.0.
"""

import json
import sys
from pathlib import Path


def key(c):
    return (c["axis"], c["pathway"], c["gene"], c["direction"], c["key_output"])


def main():
    here = Path(__file__).parent
    bench = json.loads((here / "bench_rails.json").read_text())
    engine = json.loads((here / "engine_rails.json").read_text())

    failures = []
    if engine["rail_low"] != 0.01 or engine["rail_high"] != 0.99 * bench["ceiling"]:
        failures.append(f"stops differ: engine ({engine['rail_low']}, {engine['rail_high']}) "
                        f"vs bench ceiling {bench['ceiling']}")

    e_by_key = {key(c): c for c in engine["cases"]}
    if len(e_by_key) != len(engine["cases"]):
        failures.append("duplicate keys in engine dump")
    agree = 0
    railed = 0
    for c in bench["cases"]:
        e = e_by_key.get(key(c))
        if e is None:
            failures.append(f"missing from engine: {key(c)}")
            continue
        if e["fold"] != c["fold"]:
            failures.append(f"fold differs at {key(c)}: {e['fold']} vs {c['fold']}")
        if bool(e["railed"]) != bool(c["railed"]):
            failures.append(f"flag differs at {key(c)}: fold {c['fold']}, "
                            f"engine {e['railed']}, bench {c['railed']}")
        else:
            agree += 1
            railed += bool(c["railed"])

    print(f"{len(bench['cases'])} case folds compared; agreeing: {agree} "
          f"({railed} railed) / differences: {len(failures)}")
    for f in failures[:20]:
        print("  DIFF:", f)
    if failures:
        print("THE GATE FAILS")
        sys.exit(1)
    print("no differences — the gate passes")


if __name__ == "__main__":
    main()
