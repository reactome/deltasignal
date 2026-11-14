#!/usr/bin/env python3
"""specs/154 FR-005: compare the detection pass on the frozen 97 catalog with
the identical-content rebuild (different uuids). The check: the same
(axis, pathway, gene, direction, readout) keys flag switch with the same
state, and folds agree to served precision.

Usage: compare_inv.py FROZEN_TSV REBUILD_TSV
"""
import csv, sys

def load(p):
    rows = {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        key = (r["axis"], r["pathway"], r["gene"], r["direction"], r["key_output"])
        rows[key] = r
    return rows

a = load(sys.argv[1])
b = load(sys.argv[2])

only_a = set(a) - set(b)
only_b = set(b) - set(a)
print(f"keys: frozen {len(a)}, rebuild {len(b)}, only-frozen {len(only_a)}, only-rebuild {len(only_b)}")

def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

flag_diff, state_diff = [], []
max_fold, max_other = 0.0, 0.0
switches = 0
for k in set(a) & set(b):
    ra, rb = a[k], b[k]
    if ra["switch"] != rb["switch"]:
        flag_diff.append((k, ra["switch"], rb["switch"]))
    elif ra["switch"] == "True":
        switches += 1
        if ra["state"] != rb["state"]:
            state_diff.append((k, ra["state"], rb["state"]))
        fa, fb = num(ra["other_state_fold"]), num(rb["other_state_fold"])
        if fa is not None and fb is not None:
            max_other = max(max_other, abs(fa - fb))
    fa, fb = num(ra["fold_served"]), num(rb["fold_served"])
    if fa is not None and fb is not None:
        max_fold = max(max_fold, abs(fa - fb))

print(f"switch flags differing: {len(flag_diff)}")
for k, x, y in flag_diff[:20]:
    print("  FLAG", k, x, "->", y)
print(f"states differing on shared switches ({switches} flagged in both): {len(state_diff)}")
for k, x, y in state_diff[:20]:
    print("  STATE", k, x, "->", y)
print(f"max |fold| difference: {max_fold:.3e}")
print(f"max |other_state_fold| difference on switches: {max_other:.3e}")

# Bonus parity: each arm against the recorded folds (the frozen record).
for name, rows in (("frozen97", a), ("rebuild97", b)):
    worst = 0.0
    n = 0
    for r in rows.values():
        fr, fs = num(r["fold_recorded"]), num(r["fold_served"])
        if fr is None or fs is None:
            continue
        n += 1
        worst = max(worst, abs(fr - fs))
    print(f"{name}: {n} folds vs record, max |diff| {worst:.3e}")

# The registered check (research.md FR-005): the same (pathway, readout) set
# flags, with the same STATES, to served precision; the canonical answer is
# covered by the specs/013 record (pred_ui precision). The alternate state's
# exact coordinates are reported above, not gated: a structurally symmetric
# component (SMAD2/SMAD3 in TGF-beta) holds near-degenerate alternates whose
# selection follows the specs/013 tie order — see the research.md note.
ok = not flag_diff and not state_diff and max_fold < 1e-6
print("FR-005 (flag set + states + canonical folds):", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
