#!/usr/bin/env python3
"""specs/155: fit measured-fold interval tables from the recorded fitting
datasets, and measure the registered gates out-of-sample by whole drug.

Pairs (research.md R1):
  korkut_rppa      — recorded canonical korkut scoring (pred_ui per case) x
                     the dataset's own per-condition log2 values (expr.csv);
                     one pair per (drug set, readout, condition); LODO unit =
                     drug set.
  decryptm_phospho — recorded canonical decryptM scoring x the paper's own
                     fitted curve effect sizes (all_sites.tsv, effect_size =
                     plateau ratio - 1, so measured fold = max(0, 1 + e));
                     LODO unit = drug. The cases file records no convergence
                     flag, so strata resolve at the converged-free levels.

Strata: the specs/146 machinery, imported from the record's own fitter
(stratum_key / magnitude_bin / LEVELS), MIN_CELL 30 (registered for this
spec). Arm (a): per-cell quantile intervals of the measured fold at the
nominal shares (0.5, 0.9). Arm (b): equal-count predicted-fold bands inside a
cell (each band >= MIN_CELL pairs, at most 4 bands), falling back to (a)
outside every band. Widths are measured in log2(hi/lo), lo floored at 0.01.

Gates (registered in spec.md FR-002 before this script ran):
  gate 1 (per axis, per share): the 95% drug-bootstrap interval of held-out
         coverage contains the nominal share;
  gate 2 (per axis): (b) passes gate 1 on both shares AND reduces the median
         held-out width at the 50% share by >= 10% against (a).

Usage:
  fit_155.py KORKUT_CASES_TSV DECRYPTM_CASES_TSV OUT_DIR
             [--korkut-data DIR] [--sites TSV] [--out src/api/fold_ranges.json]
Writes OUT_DIR/report_155.json always; writes the serving table only for axes
that pass gate 1, with the arm gate 2 chose.
"""

import argparse
import csv
import json
import math
import random
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2]))          # bench/
sys.path.insert(0, str(HERE.parents[1]))          # bench/analysis (score_124, holdout_report)
sys.path.insert(0, str(HERE.parents[1] / "146"))  # bench/analysis/146
import dataset_roles  # noqa: E402
from korkut_bench import conditions  # noqa: E402
from fit_confidence import LEVELS, magnitude_bin, stratum_key  # noqa: E402
from score_124 import railed  # noqa: E402

CEILING = 100.0
MIN_CELL = 30
SHARES = (0.5, 0.9)
BANDS_MAX = 4
BOOT = 10_000
WIDTH_FLOOR = 0.01
EXTERNAL = Path.home() / "deltasignal-catalogs" / "external"


def quantile(sorted_vals, q):
    """Linear-interpolated quantile of a pre-sorted list."""
    if not sorted_vals:
        raise ValueError("quantile of nothing")
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def korkut_pairs(cases_tsv, data_dir):
    dataset_roles.require_role("korkut2015", "fit")
    names = [l.strip() for l in open(data_dir / "node_Index.csv")]
    design = conditions(data_dir / "expr_index.txt")
    rows = list(csv.reader(open(data_dir / "expr.csv")))
    if len(rows) != len(design):
        raise SystemExit("expr.csv rows != expr_index.txt lines (specs/062 guard)")
    values = defaultdict(list)  # (drug set string, readout) -> [log2 ...]
    for (_, ds), row in zip(design, rows):
        if len(row) != len(names):
            raise SystemExit("expr.csv columns != node_Index.csv names (specs/062 guard)")
        for name, v in zip(names, row):
            values[("+".join(ds), name)].append(float(v))
    pairs, unmatched = [], set()
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored" or not r["pred_ui"]:
            continue
        key = (r["drug_set"], r["readout"])
        if key not in values:
            unmatched.add(key)
            continue
        fold = float(r["pred_ui"])
        conv = r.get("converged") == "1"
        for log2v in values[key]:
            pairs.append({
                "axis": "korkut_rppa", "unit": r["drug_set"],
                "pred": fold, "meas": 2.0 ** log2v,
                "call": r["predicted"], "railed": bool(railed(fold, CEILING)),
                "basis": None, "mag": magnitude_bin(fold), "conv": conv,
                "conv_known": True,
            })
    if unmatched:
        raise SystemExit(f"{len(unmatched)} korkut (drug set, readout) keys missing from "
                         f"expr.csv — the join assumption is wrong; stop and amend. "
                         f"e.g. {sorted(unmatched)[:3]}")
    return pairs


def decryptm_pairs(cases_tsv, sites_tsv):
    dataset_roles.require_role("decryptm1", "fit")
    # effect sizes per (uniprot, position, drug, cell): the paper's own fit
    effects = defaultdict(list)
    for r in csv.DictReader(open(sites_tsv), delimiter="\t"):
        if not r["effect_size"]:
            continue
        exp = r["experiment"]  # e.g. ddPTM_A549_PD325901_60min_R1
        parts = exp.split("_")
        if len(parts) < 3:
            continue
        cell, drug = parts[1], parts[2]
        effects[(r["uniprot"], r["position"], drug, cell)].append(float(r["effect_size"]))
    pairs = []
    for r in csv.DictReader(open(cases_tsv), delimiter="\t"):
        if r["status"] != "scored" or not r["pred_ui"]:
            continue
        fold = float(r["pred_ui"])
        for cell in r["cell_lines"].split(";"):
            for e in effects.get((r["uniprot"], r["position"], r["drug"], cell), ()):
                pairs.append({
                    "axis": "decryptm_phospho", "unit": r["drug"],
                    "pred": fold, "meas": max(0.0, 1.0 + e),
                    "call": r["predicted"], "railed": bool(railed(fold, CEILING)),
                    "basis": None, "mag": magnitude_bin(fold),
                    "conv": False, "conv_known": False,
                })
    return pairs


def levels_for(p):
    """The lookup walk order for this pair, the readout_confidence rule."""
    out = []
    for lv in LEVELS:
        if lv == "full" and p["basis"] is None:
            continue
        if not p["conv_known"] and lv in ("full", "no_basis", "no_magnitude"):
            continue
        out.append(lv)
    return out


def key_of(lv, p):
    return stratum_key(lv, p["call"], p["railed"], p["basis"], p["mag"], p["conv"])


def fit_cells(train):
    """Per level/key: sorted measured folds and sorted (pred, meas)."""
    cells = {lv: defaultdict(list) for lv in LEVELS}
    for p in train:
        for lv in levels_for(p):
            cells[lv][key_of(lv, p)].append((p["pred"], p["meas"]))
    return cells


def interval(meas_sorted, share):
    lo = quantile(meas_sorted, (1 - share) / 2)
    hi = quantile(meas_sorted, 1 - (1 - share) / 2)
    return lo, hi


def width(lo, hi):
    return math.log2(max(hi, WIDTH_FLOOR) / max(lo, WIDTH_FLOOR))


def cell_a(pm):
    meas = sorted(m for _, m in pm)
    return {s: interval(meas, s) for s in SHARES}


def cell_b(pm):
    """Predicted-fold bands of >= MIN_CELL pairs each, never splitting a tied
    predicted fold across bands (predicted folds are heavily discrete — railed
    0.01s above all — and splitting a tie would stratify the band by the
    MEASURED value, which is circular; the first run did exactly that and
    covered 17% at nominal 50% in-sample). Falls back to (a) when fewer than
    two bands can be formed."""
    groups = defaultdict(list)
    for pred, meas in pm:
        groups[pred].append(meas)
    distinct = sorted(groups)
    if len(distinct) < 2 or len(pm) < 2 * MIN_CELL:
        return None
    target = max(MIN_CELL, len(pm) // BANDS_MAX)
    chunks, cur, cur_n = [], [], 0
    for pred in distinct:
        cur.append(pred)
        cur_n += len(groups[pred])
        if cur_n >= target:
            chunks.append(cur)
            cur, cur_n = [], 0
    if cur:
        if chunks and cur_n < MIN_CELL:
            chunks[-1].extend(cur)
        else:
            chunks.append(cur)
    if len(chunks) < 2:
        return None
    bands = []
    for preds in chunks:
        meas = sorted(m for p in preds for m in groups[p])
        if len(meas) < MIN_CELL:
            return None
        bands.append({
            "fold_lo": preds[0], "fold_hi": preds[-1],
            "n": len(meas),
            "shares": {s: interval(meas, s) for s in SHARES},
        })
    return bands


def serve(cells, p, arm):
    """What the server would serve this pair: {share: (lo, hi)} or None."""
    for lv in levels_for(p):
        pm = cells[lv].get(key_of(lv, p))
        if pm is None or len(pm) < MIN_CELL:
            continue
        if arm == "b":
            bands = cell_b(pm)
            if bands:
                for b in bands:
                    if b["fold_lo"] <= p["pred"] <= b["fold_hi"]:
                        return b["shares"]
        return cell_a(pm)
    return None


def lodo(pairs):
    """Held-out serving per pair per arm: [(unit, share, covered, width)]."""
    units = sorted({p["unit"] for p in pairs})
    out = {arm: [] for arm in ("a", "b")}
    by_unit = defaultdict(list)
    for p in pairs:
        by_unit[p["unit"]].append(p)
    for u in units:
        train = [p for p in pairs if p["unit"] != u]
        cells = fit_cells(train)
        for p in by_unit[u]:
            for arm in ("a", "b"):
                served = serve(cells, p, arm)
                if served is None:
                    continue
                for s, (lo, hi) in served.items():
                    out[arm].append((u, s, lo <= p["meas"] <= hi, width(lo, hi)))
    return out


def bootstrap_coverage(rows, share, seed):
    """95% percentile interval of coverage, resampling units (the only real
    independence structure: conditions of one drug set are not independent)."""
    per_unit = defaultdict(lambda: [0, 0])
    for u, s, cov, _ in rows:
        if s != share:
            continue
        per_unit[u][0] += 1
        per_unit[u][1] += cov
    units = sorted(per_unit)
    if not units:
        return None
    rng = random.Random(seed)
    stats = []
    for _ in range(BOOT):
        n = c = 0
        for u in rng.choices(units, k=len(units)):
            n += per_unit[u][0]
            c += per_unit[u][1]
        stats.append(c / n if n else 0.0)
    stats.sort()
    point = sum(v[1] for v in per_unit.values()) / sum(v[0] for v in per_unit.values())
    return {"point": point, "lo": quantile(stats, 0.025), "hi": quantile(stats, 0.975),
            "units": len(units), "pairs": sum(v[0] for v in per_unit.values())}


def median(vals):
    vals = sorted(vals)
    return vals[len(vals) // 2] if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("korkut_cases")
    ap.add_argument("decryptm_cases")
    ap.add_argument("out_dir")
    ap.add_argument("--korkut-data", default=str(EXTERNAL / "korkut2015"))
    ap.add_argument("--sites", default=str(EXTERNAL / "decryptm" / "all_sites.tsv"))
    ap.add_argument("--out", default=None, help="serving table path (omit to only report)")
    ap.add_argument("--catalog-build", default=None,
                    help="the serving record the table is stamped for (required with --out); "
                         "the pairs' predicted folds are the frozen-97 record, which specs/152 "
                         "proved the sound Release-98 record reproduces exactly")
    ap.add_argument("--solver-commit", default=None)
    ap.add_argument("--seed", type=int, default=155)
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    axes = {
        "korkut_rppa": korkut_pairs(args.korkut_cases, Path(args.korkut_data)),
        "decryptm_phospho": decryptm_pairs(args.decryptm_cases, args.sites),
    }
    report = {"spec": "specs/155-fold-intervals", "date": date.today().isoformat(),
              "min_cell": MIN_CELL, "shares": list(SHARES), "boot": BOOT,
              "width": "log2(hi/lo), floors at 0.01", "axes": {}}
    passing = {}
    for axis, pairs in axes.items():
        print(f"{axis}: {len(pairs)} pairs, {len({p['unit'] for p in pairs})} units", flush=True)
        if not pairs:
            report["axes"][axis] = {"pairs": 0, "verdict": "no pairs"}
            continue
        held = lodo(pairs)
        ax = {"pairs": len(pairs), "units": len({p["unit"] for p in pairs}), "arms": {}}
        for arm in ("a", "b"):
            arm_r = {}
            for s in SHARES:
                arm_r[str(s)] = bootstrap_coverage(held[arm], s, args.seed)
            arm_r["median_width_50"] = median([w for _, s, _, w in held[arm] if s == 0.5])
            arm_r["served_pairs"] = len({(u, s) for u, s, _, _ in held[arm]})
            ax["arms"][arm] = arm_r
        g1 = {arm: all(ax["arms"][arm][str(s)] and
                       ax["arms"][arm][str(s)]["lo"] <= s <= ax["arms"][arm][str(s)]["hi"]
                       for s in SHARES)
              for arm in ("a", "b")}
        wa, wb = ax["arms"]["a"]["median_width_50"], ax["arms"]["b"]["median_width_50"]
        reduction = (wa - wb) / wa if wa and wb is not None else None
        g2 = bool(g1["b"] and reduction is not None and reduction >= 0.10)
        ax["gate1"] = g1
        ax["gate2_width_reduction"] = reduction
        ax["gate2"] = g2
        ax["verdict"] = ("serve (b)" if g2 and g1["b"] else
                         "serve (a)" if g1["a"] else "serve nothing")
        report["axes"][axis] = ax
        if ax["verdict"] != "serve nothing":
            passing[axis] = "b" if ax["verdict"] == "serve (b)" else "a"
        print(f"  gate1 a={g1['a']} b={g1['b']} gate2={g2} "
              f"(width reduction {reduction if reduction is None else round(reduction, 4)}) "
              f"-> {ax['verdict']}", flush=True)
        for arm in ("a", "b"):
            for s in SHARES:
                print(f"    arm {arm} share {s}: {ax['arms'][arm][str(s)]}", flush=True)

    (out_dir / "report_155.json").write_text(json.dumps(report, indent=1))
    print(f"report -> {out_dir / 'report_155.json'}")

    if args.out and passing:
        if not args.catalog_build or not args.solver_commit:
            raise SystemExit("--out needs --catalog-build and --solver-commit (FR-003)")
        table = {"version": date.today().isoformat(), "spec": "specs/155-fold-intervals",
                 "catalog_build": args.catalog_build, "solver_commit": args.solver_commit,
                 "fit_record": {"korkut_cases": str(args.korkut_cases),
                                "decryptm_cases": str(args.decryptm_cases),
                                "datasets": "korkut2015 (FIT), decryptm1 (FIT)",
                                "prediction_identity": "pairs carry the frozen-97 record's "
                                "predicted folds; specs/152 proved the served Release-98 "
                                "record reproduces them exactly (zero differences)"},
                 "min_cell": MIN_CELL, "nominal_shares": list(SHARES),
                 "gates": report["axes"], "arms": passing, "table": {}}
        for axis, arm in passing.items():
            cells = fit_cells(axes[axis])
            ax_t = {}
            for lv in LEVELS:
                lv_t = {}
                for key, pm in cells[lv].items():
                    if len(pm) < MIN_CELL:
                        continue
                    a = cell_a(pm)
                    entry = {"n": len(pm),
                             "shares": {str(s): {"lo": round(a[s][0], 6), "hi": round(a[s][1], 6)}
                                        for s in SHARES}}
                    if arm == "b":
                        bands = cell_b(pm)
                        if bands:
                            entry["bands"] = [
                                {"fold_lo": round(b["fold_lo"], 6), "fold_hi": round(b["fold_hi"], 6),
                                 "n": b["n"],
                                 "shares": {str(s): {"lo": round(b["shares"][s][0], 6),
                                                     "hi": round(b["shares"][s][1], 6)}
                                            for s in SHARES}}
                                for b in bands]
                    lv_t["|".join(str(x) for x in key)] = entry
                if lv_t:
                    ax_t[lv] = lv_t
            table["table"][axis] = ax_t
        Path(args.out).write_text(json.dumps(table, indent=1))
        print(f"serving table ({passing}) -> {args.out}")
    elif args.out:
        print("no axis passed gate 1: no serving table written (the honest outcome)")


if __name__ == "__main__":
    main()
