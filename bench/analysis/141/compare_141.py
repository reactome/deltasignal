#!/usr/bin/env python3
"""specs/141 FR-003/FR-004: the parity gate. Compares the bench's readout
resolutions (dump_readouts.py, the bench's own functions) with the engine's
(readouts_141.jl) field by field: curator basis + curated-entity groups,
korkut occurrence and kept-copy lists. ANY difference is listed and fails the
gate (non-zero exit).

  python bench/analysis/141/compare_141.py \
    bench/analysis/141/bench_readouts.json bench/analysis/141/engine_readouts.json
"""
import json
import sys
from pathlib import Path


def main():
    bench_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "bench_readouts.json"
    engine_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "engine_readouts.json"
    bench = json.loads(bench_path.read_text())
    engine = json.loads(engine_path.read_text())

    diffs = []
    total = agree = 0
    basis_counts = {}

    eng_cur = {p["dir"]: {q["ko"]: q for q in p["readouts"]} for p in engine["curator"]}
    for p in bench["curator"]:
        ep = eng_cur.get(p["dir"])
        if ep is None:
            diffs.append(f"curator {p['dir']}: missing from the engine dump")
            continue
        for q in p["readouts"]:
            total += 1
            eq = ep.get(q["ko"])
            if eq is None:
                diffs.append(f"curator {p['dir']} ko {q['ko']}: missing from the engine dump")
                continue
            problems = []
            if q["basis"] != eq["basis"]:
                problems.append(f"basis {q['basis']} != {eq['basis']}")
            if q["groups"] != eq["groups"]:
                problems.append(f"groups differ: bench {q['groups']} engine {eq['groups']}")
            if problems:
                diffs.append(f"curator {p['dir']} ko {q['ko']}: " + "; ".join(problems))
            else:
                agree += 1
                basis_counts[q["basis"]] = basis_counts.get(q["basis"], 0) + 1

    k_total = k_agree = 0
    eng_kor = {p["dir"]: {q["readout"]: q for q in p["readouts"]} for p in engine["korkut"]}
    for p in bench["korkut"]:
        ep = eng_kor.get(p["dir"])
        if ep is None:
            diffs.append(f"korkut {p['dir']}: missing from the engine dump")
            continue
        for q in p["readouts"]:
            k_total += 1
            eq = ep.get(q["readout"])
            if eq is None:
                diffs.append(f"korkut {p['dir']} {q['readout']}: missing from the engine dump")
                continue
            problems = []
            if q["occurrences"] != eq["occurrences"]:
                problems.append("occurrences differ")
            if q["kept"] != eq["kept"]:
                problems.append(f"kept differ: bench {q['kept']} engine {eq['kept']}")
            if problems:
                diffs.append(f"korkut {p['dir']} {q['readout']}: " + "; ".join(problems))
            else:
                k_agree += 1

    basis_str = ", ".join(f"{b} {n}" for b, n in sorted(basis_counts.items()))
    print(f"curator: {total} (pathway, key_output) readouts compared; agreeing: {agree} ({basis_str})")
    print(f"korkut:  {k_total} (pathway, total readout) cases compared; agreeing: {k_agree}")
    if diffs:
        print(f"\n{len(diffs)} DIFFERENCES — the gate FAILS:")
        for d in diffs:
            print(f"  {d}")
        sys.exit(1)
    print("no differences — the gate passes")


if __name__ == "__main__":
    main()
