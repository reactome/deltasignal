#!/usr/bin/env python3
"""Dose-response: does moving an input further from baseline move the outputs
further from baseline? Pre-registered in specs/021.

The benchmark perturbs every gene at one fixed strength, so this reruns the
identical pipeline at graded strengths (DS_PERTURB_UI_DOWN / _UP) and joins the
results per (pathway, gene, direction, readout).

  M1 monotonicity -- the pass/fail. Across increasing input strength a readout
     should move consistently (up, or down if it is inhibited). A reversal is a
     solver defect.
  M2 graded vs switched -- is a moving readout already at a rail (0 or the 100
     cap) at the mildest input? If most are, magnitude carries no dose
     information for them.
  M3 transfer -- for readouts that never rail, the slope of |log output fold|
     against |log input fold|.

Usage:
  python bench/analysis/dose_response.py --dir <results>/dose
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from holdout_report import TUNING_PATHWAYS  # noqa: E402

# (knockdown value, overexpression value) per run, mildest first.
LADDER = [(0.5, 2.0), (0.2, 5.0), (0.05, 20.0), (0.0, 80.0)]
BASELINE = 1.0
MOVE_TOL = 1e-6       # an output within this of baseline has not moved
MONO_TOL = 1e-9       # the pre-registered tolerance between steps
PRINT_UNIT = 1e-6     # pred_ui is written with %.6f, so steps below this are
                      # print noise; reported as a sensitivity, not the decision
IDENTITY_REL = 1e-4   # an output within this RELATIVE distance of the input at
                      # every step copies it (absolute would be 2e-3 relative at
                      # KD 0.05 but 1e-6 at OE 80)
LOG_PRINT_TOL = 1e-5  # |log10 fold| steps smaller than this are print noise
ZERO_RAIL = 1e-6      # output at or below this is railed at zero
CAP_RAIL = 99.99      # output at or above this is railed at the 100 cap
ZERO_FLOOR = 1e-6     # for logs of a zero output or input


def fname(down: float, up: float) -> str:
    fmt = lambda v: f"{v:g}"
    return f"cases_down{fmt(down)}_up{fmt(up)}.tsv"


def load_ladder(directory: str) -> dict[tuple, list[float]]:
    """case key -> [output at each strength, mildest first], valid cases only."""
    per_run = []
    for down, up in LADDER:
        path = os.path.join(directory, fname(down, up))
        with open(path, newline="") as fh:
            rows = {}
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("valid") != "1" or r["direction"] not in ("0", "2"):
                    continue
                rows[(r["pathway"], r["gene"], r["direction"], r["key_output"])] = float(r["pred_ui"])
        per_run.append(rows)
    common = set(per_run[0]).intersection(*per_run[1:])
    union = set().union(*per_run)
    if len(union) != len(common):
        # A case valid at one strength but not another would bias M1 silently.
        print(f"WARNING: {len(union) - len(common)} cases are not valid in every run "
              f"and are dropped ({len(common)} kept)", file=sys.stderr)
    return {k: [run[k] for run in per_run] for k in common}


def input_values(direction: str) -> list[float]:
    idx = 0 if direction == "0" else 1
    return [pair[idx] for pair in LADDER]


def moves(outs: list[float]) -> bool:
    return any(abs(o - BASELINE) > MOVE_TOL for o in outs)


def monotone(outs: list[float], tol: float = MONO_TOL) -> bool:
    """Consistent in one direction as the input strengthens (either way).
    NOTE: this is the pre-registered M1 and checks the raw value only; it does
    not require the output to get further from baseline (see `away`)."""
    up = all(b >= a - tol for a, b in zip(outs, outs[1:]))
    down = all(b <= a + tol for a, b in zip(outs, outs[1:]))
    return up or down


def away(outs: list[float], tol: float = MONO_TOL, log_tol: float = 1e-9) -> bool:
    """What the question actually asks: every output on one side of baseline,
    and |log fold| never shrinking as the input strengthens."""
    if not (all(o >= BASELINE - tol for o in outs) or all(o <= BASELINE + tol for o in outs)):
        return False
    d = [abs(math.log10(max(o, ZERO_FLOOR) / BASELINE)) for o in outs]
    return all(b >= a - log_tol for a, b in zip(d, d[1:]))


def identity(ins: list[float], outs: list[float]) -> bool:
    """The readout equals the pinned input at every step: it is the perturbed
    node itself, or a single-input pass-through. Monotone by construction, so it
    cannot test anything and is reported separately."""
    return all(abs(o - i) <= max(IDENTITY_REL * i, PRINT_UNIT) for i, o in zip(ins, outs))


def railed(v: float) -> bool:
    return v <= ZERO_RAIL or v >= CAP_RAIL


def transfer_slope(ins: list[float], outs: list[float]) -> float | None:
    """Least-squares slope of |log10 output fold| on |log10 input fold|, over
    the steps with a finite input. A knockdown to 0 has no finite log; flooring
    it at 1e-6 put one point at x = 6 against x <= 1.3 and dominated the fit."""
    pts = [(i, o) for i, o in zip(ins, outs) if i > 0]
    x = np.array([abs(math.log10(i / BASELINE)) for i, _ in pts])
    y = np.array([abs(math.log10(max(o, ZERO_FLOOR) / BASELINE)) for _, o in pts])
    if len(x) < 2:
        return None
    if np.ptp(x) == 0:
        return None
    return float(np.polyfit(x, y, 1)[0])


def summarise(cases: dict[tuple, list[float]]) -> dict:
    moving = {k: v for k, v in cases.items() if moves(v)}
    ident = {k for k, v in moving.items() if identity(input_values(k[2]), v)}
    nontrivial = {k: v for k, v in moving.items() if k not in ident}
    mono = {k: v for k, v in moving.items() if monotone(v)}
    # M2 on the non-identity readouts: an input copy is never railed at 2x, so
    # counting it would halve the railed fraction for no reason.
    railed_mild = [k for k, v in nontrivial.items() if railed(v[0])]
    changes = [k for k, v in nontrivial.items() if abs(v[-1] - v[0]) > MOVE_TOL]
    never_railed = {k: v for k, v in nontrivial.items() if not any(railed(o) for o in v)}
    slopes = [s for k, v in never_railed.items()
              if (s := transfer_slope(input_values(k[2]), v)) is not None]
    reversals = sorted(
        ((k, v) for k, v in nontrivial.items() if not monotone(v, PRINT_UNIT)),
        key=lambda kv: -(max(kv[1]) - min(kv[1])))
    return {
        "cases": len(cases), "moving": len(moving), "monotone": len(mono),
        "identity": len(ident), "nontrivial": len(nontrivial),
        "nt_monotone": sum(monotone(v) for v in nontrivial.values()),
        "nt_monotone_print": sum(monotone(v, PRINT_UNIT) for v in nontrivial.values()),
        "nt_away": sum(away(v) for v in nontrivial.values()),
        "nt_away_print": sum(away(v, PRINT_UNIT, LOG_PRINT_TOL) for v in nontrivial.values()),
        "railed_at_mildest": len(railed_mild), "change_mild_to_strong": len(changes),
        "never_railed": len(never_railed), "slopes": slopes, "reversals": reversals,
    }


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:5.1f}%" if b else "  n/a"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    a = ap.parse_args()
    cases = load_ladder(a.dir)
    split = defaultdict(dict)
    for k, v in cases.items():
        split["HELD-OUT" if k[0] not in TUNING_PATHWAYS else "TUNING"][k] = v
        split["ALL"][k] = v

    for label in ("ALL", "HELD-OUT", "TUNING"):
        s = summarise(split[label])
        print(f"\n{label}: {s['cases']:,} resolved cases, {s['moving']:,} move at some strength")
        print(f"  M1 monotone in input strength     {s['monotone']:>6,}  {pct(s['monotone'], s['moving'])}"
              f"   (P1 needs >= 95%)")
        n = s["nontrivial"]
        print(f"     of which copy the input (identity) {s['identity']:>6,}  -- monotone by construction")
        print(f"  M1 on the {n:,} non-identity readouts {s['nt_monotone']:>6,}  {pct(s['nt_monotone'], n)}"
              f"   (tolerance one print unit: {pct(s['nt_monotone_print'], n)})")
        print(f"     further from baseline, same side  {s['nt_away']:>6,}  {pct(s['nt_away'], n)}"
              f"   (print tolerance: {pct(s['nt_away_print'], n)})")
        print(f"  M2 non-identity railed at mildest {s['railed_at_mildest']:>6,}  {pct(s['railed_at_mildest'], n)}"
              f"   (> 50% = effectively a switch)")
        print(f"     output changes mild -> strong   {s['change_mild_to_strong']:>6,}  {pct(s['change_mild_to_strong'], n)}")
        if s["slopes"]:
            q = np.percentile(s["slopes"], [10, 25, 50, 75, 90])
            neg = sum(x < 0 for x in s["slopes"])
            print(f"     never railed: {len(s['slopes']):,} of {n:,} (slopes below exclude the strongest"
                  f" transmitters); negative slopes {neg:,} ({pct(neg, len(s['slopes'])).strip()})")
            print(f"  M3 transfer slope, {len(s['slopes']):,} never-railed non-identity, finite steps:  "
                  f"p10 {q[0]:.2f}  p25 {q[1]:.2f}  median {q[2]:.2f}  p75 {q[3]:.2f}  p90 {q[4]:.2f}")
        if label == "ALL" and s["reversals"]:
            print(f"  largest reversals ({len(s['reversals'])} non-identity, beyond one print unit):")
            for (pw, gene, d, ko), outs in s["reversals"][:8]:
                ins = input_values(d)
                path = "  ".join(f"{i:g}->{o:.3g}" for i, o in zip(ins, outs))
                print(f"    {pw[:34]:<34} {gene:<8} {'KD' if d == '0' else 'OE'}  ko {ko:<8} {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
