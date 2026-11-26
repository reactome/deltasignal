"""specs/173 amendment 1: the binned measured call and the split, on synthetic samples."""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lun_a_173 as M  # noqa: E402
from common_173 import split_of  # noqa: E402

HEADER = "experiment,marker,perturbation,row_col,timepoint,origin,target,bin,counts,fit_var,median,overall_var,var_ratio\n"


def gene_in(half):
    return next(g for g in (f"G{i}" for i in range(1000)) if split_of(g) == half)


def rows(exp, marker, well, target, values, counts=50, tp="0"):
    return [f"{exp},{marker},EGF,{well},{tp},GFP,{target},{b},{counts},0,{v},0,0\n" for b, v in enumerate(values)]


def flat(level):
    return [level] * 10


def ramp(lo, hi):
    return [lo + (hi - lo) * b / 9 for b in range(10)]


def write(path, lines):
    path.write_text(HEADER + "".join(lines))


def test_dose_slope_high_minus_low_and_cell_floor():
    bins = {b: (50, v) for b, v in enumerate(ramp(1.0, 2.0))}
    assert M.dose_slope(bins) == pytest.approx(ramp(1.0, 2.0)[8] - ramp(1.0, 2.0)[1])
    thin = {b: (5 if b >= 7 else 50, v) for b, (c, v) in bins.items()}
    assert M.dose_slope(thin) is None


def test_split_routes_controls_to_both_and_unmapped_to_b(tmp_path):
    a, b = gene_in("A"), gene_in("B")
    genes = {"MA": a, "MB": b}
    src = tmp_path / "bin.csv"
    write(src, rows("E1", "FLAG-GFP1", "A1", "p-AKT", flat(1)) + rows("E1", "MA", "A2", "p-AKT", flat(1))
          + rows("E1", "MB", "A3", "p-AKT", flat(1)) + rows("E1", "MX", "A4", "p-AKT", flat(1)))
    counts = M.split(src, tmp_path / "a.csv", tmp_path / "b.csv", genes)
    assert counts == {"control": 10, "A": 10, "B": 10, "B_unmapped": 10}
    seen_a = {r["marker"] for r in csv.DictReader(open(tmp_path / "a.csv"))}
    seen_b = {r["marker"] for r in csv.DictReader(open(tmp_path / "b.csv"))}
    assert seen_a == {"FLAG-GFP1", "MA"} and seen_b == {"FLAG-GFP1", "MB", "MX"}


def test_reader_refuses_half_b(tmp_path):
    b = gene_in("B")
    f = tmp_path / "a.csv"
    write(f, rows("E1", "MB", "A3", "p-AKT", flat(1)))
    with pytest.raises(SystemExit):
        M.load_a(f, {"MB": b})


def test_calls_against_plate_controls_and_threshold(tmp_path):
    a = gene_in("A")
    lines = []
    for i in range(12):  # control slopes spread 0..0.11 on one plate
        lines += rows("E1", f"FLAG-GFP{i}", f"C{i}", "p-AKT", ramp(1.0, 1.0 + 0.01 * i))
    lines += rows("E1", "MA", "X1", "p-AKT", ramp(1.0, 3.0))       # strong rise with dose
    lines += rows("E1", "MA", "X1", "p-ERK1-2", ramp(1.0, 1.0))    # no ERK controls: no threshold, not called
    lines += rows("E1", "MA", "X2", "p-AKT", ramp(1.0, 3.0), tp="10")  # EGF timepoint ignored
    f = tmp_path / "a.csv"
    write(f, lines)
    sl = M.slopes(M.load_a(f, {"MA": a}))
    delta = M.thresholds(sl)
    assert set(delta) == {"p-AKT"} and 0 < delta["p-AKT"] < 0.2
    calls = M.measured_calls(sl, {"MA": a}, delta)
    assert calls[(a, "p-AKT")][0] == "up" and (a, "p-ERK1/2") not in calls
    down = M.measured_calls({**sl, ("E1", "MA", "X1"): {"p-AKT": -1.0}}, {"MA": a}, delta)
    assert down[(a, "p-AKT")][0] == "down"


def test_enrichment_and_lift():
    cases = ([{"gene": "g", "reach": "up", "frozen": "up", "measured": "up"}] * 3
             + [{"gene": "g", "reach": "up", "frozen": "not", "measured": "not"}] * 1
             + [{"gene": "h", "reach": "not", "frozen": "not", "measured": "up"}] * 1
             + [{"gene": "h", "reach": "not", "frozen": "not", "measured": "not"}] * 3)
    assert M.enrichment(cases) == pytest.approx(0.75 / 0.25)
    assert M.lift(cases) == pytest.approx(1.0 - 0.75)
