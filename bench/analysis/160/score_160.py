#!/usr/bin/env python3
"""specs/160: score the per-protein turnover bound on the recorded korkut TSV.

For each scored row on a TOTAL (abundance) readout, the predicted measured fold
at T = 24 h is L = f + (1 - f) * 2^(-T / t_eff), with f the recorded pred_ui
and t_eff the median effective half-life of the readout's gene(s) from
specs/160-turnover-bound/halflives.tsv; a readout with no covered gene is
unchanged. Calls at the bench's cutoffs (0.85 / 1.15). Phospho rows keep
their recorded call. Control: the same rule with every protein at the 46 h
anchor (tier 1). Gate (spec.md US1): balanced accuracy beats the recorded
calls, 95% bootstrap over drug sets (10k, seed 160), interval excluding zero.

Usage: score_160.py KORKUT_CASES_TSV HALFLIVES_TSV
"""
import csv
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import dataset_roles  # noqa: E402
import korkut_bench as kb  # noqa: E402

T_H = 24.0
ANCHOR_H = 46.0
DOWN_CUT, UP_CUT = 0.85, 1.15
N_BOOT = 10_000
SEED = 160
RAF = "R-HSA-5673001"


def call(fold):
    return "down" if fold < DOWN_CUT else ("up" if fold >= UP_CUT else "not")


def bound(f, t_eff):
    return f + (1.0 - f) * 2.0 ** (-T_H / t_eff)


def ba(rows, col):
    per = []
    for m in ("down", "not", "up"):
        sub = [r for r in rows if r["measured"] == m]
        if sub:
            per.append(sum(r[col] == m for r in sub) / len(sub))
    return sum(per) / len(per)


def acc(rows, col):
    return sum(r[col] == r["measured"] for r in rows) / len(rows)


def load_halflives(path):
    hl = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("gene\t"):
            continue
        g, _, _, te = line.rstrip("\n").split("\t")
        hl[g] = float(te)
    return hl


def gate(rows, col, rng):
    by_set = defaultdict(list)
    for r in rows:
        by_set[r["drug_set"]].append(r)
    sets = sorted(by_set)
    deltas = []
    for _ in range(N_BOOT):
        picked = [by_set[rng.choice(sets)] for _ in sets]
        flat = [r for g in picked for r in g]
        deltas.append(ba(flat, col) - ba(flat, "predicted"))
    deltas.sort()
    return (ba(rows, col) - ba(rows, "predicted"),
            deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1], len(sets))


def main():
    dataset_roles.require_role("korkut2015", "fit")
    cases_tsv, hl_path = sys.argv[1], sys.argv[2]
    hl = load_halflives(hl_path)
    readouts = {r["readout"]: r for r in kb.read_tsv(kb.MAPS / "readouts.tsv")}
    t_of = {}
    for name, r in readouts.items():
        if r["kind"] != "total":
            continue
        ts = [hl[g] for g in r["sites"].split(",") if g in hl]
        t_of[name] = statistics.median(ts) if ts else None

    rows = []
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored":
            continue
        f = float(r["pred_ui"])
        name = r["readout"]
        rec = r["predicted"]
        if name in t_of:
            t = t_of[name]
            rule = call(bound(f, t)) if t is not None else rec
            ctrl = call(bound(f, ANCHOR_H))
        else:
            rule = ctrl = rec
        rows.append({"drug_set": r["drug_set"], "pathway": r["pathway_id"],
                     "readout": name, "measured": r["measured"], "predicted": rec,
                     "rule": rule, "ctrl": ctrl, "pred_ui": f,
                     "total": name in t_of})

    rng = random.Random(SEED)
    print(f"{len(rows)} scored rows; total-readout rows {sum(r['total'] for r in rows)}; "
          f"covered total readouts {sum(1 for t in t_of.values() if t is not None)} of {len(t_of)}")
    print(f"recorded:              acc {acc(rows, 'predicted'):.4f}  balanced {ba(rows, 'predicted'):.4f}")
    for col, name in (("rule", "per-protein bound"), ("ctrl", "tier-1 control (46 h)")):
        d, lo, hi, n = gate(rows, col, rng)
        print(f"{name:>22}: acc {acc(rows, col):.4f}  balanced {ba(rows, col):.4f}  "
              f"delta {d:+.4f} [{lo:+.4f}, {hi:+.4f}] over {n} drug sets"
              + (f"  -> {'PASS' if lo > 0 else 'FAIL'}" if col == "rule" else ""))
        moved = [r for r in rows if r[col] != r["predicted"]]
        fixed = sum(r[col] == r["measured"] for r in moved)
        broke = sum(r["predicted"] == r["measured"] for r in moved)
        print(f"{'':>24}calls moved {len(moved)}: fixed {fixed}, broken {broke}")
    nr = [r for r in rows if r["pathway"] != RAF]
    d, lo, hi, _ = gate(nr, "rule", random.Random(SEED))
    print(f"  (RAF set aside: per-protein delta {d:+.4f} [{lo:+.4f}, {hi:+.4f}])")

    railed = [r for r in rows if r["total"] and (r["pred_ui"] <= 0.01 or r["pred_ui"] >= 99)]
    print(f"\nrailed total-readout rows {len(railed)}: recorded right "
          f"{sum(r['predicted'] == r['measured'] for r in railed)}, per-protein right "
          f"{sum(r['rule'] == r['measured'] for r in railed)}")


if __name__ == "__main__":
    main()
