#!/usr/bin/env python3
"""specs/148 scoring: the F discipline against the frozen control.

The pathway's own (control) answer stands wherever the pathway alone routes
the readout (control fold != 1.0); the arrow-bridged arm answers only the
control-untouched readouts. H1: korkut balanced accuracy on that class vs
always-no-change, by-drug-set bootstrap. H2: per curator/experimental axis,
fixed vs broken among control-untouched cases.

Usage: score_148.py CONTROL_DIR KORKUT_CONTROL_TSV ARM_SUBSET_DIR ARM_KORKUT_DIR
"""
import csv
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

control_dir, kctl, arm_s, arm_k = map(Path, sys.argv[1:5])
CENTRES = {"Signaling_by_ERBB2", "Cell_Cycle_Checkpoints", "Mitotic_G1-G1_S_phases",
           "Intrinsic_Pathway_for_Apoptosis", "Transcriptional_regulation_by_RUNX1",
           "Transcriptional_regulation_by_RUNX2"}
CENTRE_DIRS = {"R-HSA-1227986", "R-HSA-69620", "R-HSA-453279", "R-HSA-109606",
               "R-HSA-8878166", "R-HSA-8878171"}


def load(path, key_cols, val_cols, filt=None):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if filt and not filt(r):
            continue
        out[tuple(r[c] for c in key_cols)] = {c: r[c] for c in val_cols}
    return out


def untouched(fold):
    return float(fold) == 1.0


print("== curator / experimental (H2) ==")
for axis in ("curator", "experimental"):
    kc = ("pathway", "gene", "direction", "key_output")
    vc = ("predicted", "expected", "pred_ui", "valid")
    ctl = load(control_dir / f"{axis}_cases.tsv", kc, vc, lambda r: r["pathway"] in CENTRES)
    arm = load(arm_s / f"{axis}_cases.tsv", kc, vc)
    both = [k for k in ctl if k in arm and ctl[k]["valid"] == "1" and arm[k]["valid"] == "1"]
    un = [k for k in both if untouched(ctl[k]["pred_ui"])]
    fixed = sum(1 for k in un if arm[k]["predicted"] == ctl[k]["expected"] != ctl[k]["predicted"])
    broke = sum(1 for k in un if ctl[k]["predicted"] == ctl[k]["expected"] != arm[k]["predicted"])
    changed = sum(1 for k in un if arm[k]["predicted"] != ctl[k]["predicted"])
    acc_c = sum(ctl[k]["predicted"] == ctl[k]["expected"] for k in both) / len(both)
    f_pred = {k: (ctl[k]["predicted"] if not untouched(ctl[k]["pred_ui"]) else arm[k]["predicted"]) for k in both}
    acc_f = sum(f_pred[k] == ctl[k]["expected"] for k in both) / len(both)
    print(f"{axis}: cases {len(both)} (untouched {len(un)}); F changed {changed}, "
          f"fixed {fixed} / broke {broke}; accuracy control {acc_c:.4f} -> F {acc_f:.4f}")

print("\n== korkut (H1) ==")
kc = ("pathway_id", "drug_set", "readout")
vc = ("predicted", "measured", "pred_ui", "status")
ctl = load(kctl, kc, vc, lambda r: r["pathway_id"] in CENTRE_DIRS and r["status"] == "scored" and r["measured"])
arm = load(arm_k / "korkut_cases.tsv", kc, vc, lambda r: r["status"] == "scored")
both = [k for k in ctl if k in arm]
un = [k for k in both if untouched(ctl[k]["pred_ui"])]
print(f"cases {len(both)}; control-untouched {len(un)}; arm answers differing there "
      f"{sum(1 for k in un if arm[k]['predicted'] != ctl[k]['predicted'])}")


def balanced(keys, pred):
    per = defaultdict(lambda: [0, 0])
    for k in keys:
        t = ctl[k]["measured"]
        per[t][0] += pred(k) == t
        per[t][1] += 1
    rec = [a / b for a, b in per.values() if b]
    return sum(rec) / len(rec)


if un:
    b_arm = balanced(un, lambda k: arm[k]["predicted"])
    b_not = balanced(un, lambda k: "not")
    sets = sorted({k[1] for k in un})
    by = defaultdict(list)
    for k in un:
        by[k[1]].append(k)
    rng = random.Random(148)
    diffs = []
    for _ in range(10000):
        ks = [k for s in (rng.choice(sets) for _ in sets) for k in by[s]]
        diffs.append(balanced(ks, lambda k: arm[k]["predicted"]) - balanced(ks, lambda k: "not"))
    diffs.sort()
    print(f"H1 untouched class: arm balanced {b_arm:.4f} vs always-NOT {b_not:.4f} "
          f"(diff {b_arm-b_not:+.4f}, 95% [{diffs[249]:+.4f}, {diffs[9749]:+.4f}])"
          f" -> {'PASSES' if diffs[249] > 0 else 'FAILS'}")
    mc = Counter((ctl[k]["measured"], arm[k]["predicted"]) for k in un)
    print("untouched (measured, arm_pred):", dict(mc))
