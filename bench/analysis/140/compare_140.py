#!/usr/bin/env python3
"""specs/140 FR-003, the gate: the engine's pins EQUAL the bench's pins (as
sets), and the basis agrees, for every benchmark (pathway, gene). Any
difference is listed in full and the exit code is non-zero: a port that is
"close" is a second protocol.

  python bench/analysis/140/compare_140.py \
    bench/analysis/140/bench_pins.json bench/analysis/140/engine_pins.json
"""
import json
import sys
from pathlib import Path


def by_key(dump):
    out = {}
    for p in dump:
        for q in p["queries"]:
            out[(p["dir"], q["gene"])] = q
    return out


def main():
    bench_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "bench_pins.json"
    engine_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "engine_pins.json"
    bench = by_key(json.loads(bench_path.read_text()))
    engine = by_key(json.loads(engine_path.read_text()))

    diffs = []
    if set(bench) != set(engine):
        diffs.append(f"case sets differ: bench-only {sorted(set(bench) - set(engine))[:5]}, "
                     f"engine-only {sorted(set(engine) - set(bench))[:5]}")
    agree_basis = {"root": 0, "root_cycle": 0, "none": 0}
    for key in sorted(set(bench) & set(engine)):
        b, e = bench[key], engine[key]
        if sorted(set(b["pins"])) != sorted(set(e["pins"])):
            diffs.append(f"{key}: pins differ — bench {sorted(set(b['pins']))} vs engine {sorted(set(e['pins']))}")
        elif b["basis"] != e["basis"]:
            diffs.append(f"{key}: basis differs — bench {b['basis']} vs engine {e['basis']}")
        elif b["occurrences"] != e["occurrences"]:
            diffs.append(f"{key}: occurrences differ — bench {b['occurrences']} vs engine {e['occurrences']}")
        else:
            agree_basis[b["basis"]] += 1

    n = len(set(bench) & set(engine))
    print(f"{n} (pathway, gene) cases compared; agreeing: {sum(agree_basis.values())} "
          f"(root {agree_basis['root']}, root_cycle {agree_basis['root_cycle']}, none {agree_basis['none']})")
    if diffs:
        print(f"\n{len(diffs)} DIFFERENCES — the gate FAILS:")
        for d in diffs:
            print(f"  {d}")
        sys.exit(1)
    print("no differences — the gate passes")


if __name__ == "__main__":
    main()
