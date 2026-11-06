"""specs/121 amendment 6: the final scorer's pieces (values-free fixtures)."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_121_final as sf  # noqa: E402

UP = math.log2(2)


def test_parse_minutes():
    assert sf.parse_minutes("5min") == 5
    assert sf.parse_minutes("30 min") == 30
    assert sf.parse_minutes("0.5hr") == 30
    assert sf.parse_minutes("4hr") == 240
    assert sf.parse_minutes("24 hrs") == 1440
    assert sf.parse_minutes("3d") == 3 * 1440
    assert sf.parse_minutes("2 days") == 2 * 1440
    assert sf.parse_minutes("overnight") is None
    assert sf.parse_minutes("") is None
    assert sf.parse_minutes(None) is None


def test_gated_arm_falls_back_to_frozen():
    preds = {"fold": {("D", "x"): 1.0, ("D", "y"): -1.0},
             "reach": {("D", "x"): -1.0, ("D", "y"): 0.0},
             "gated": ({("D", "c1", "x"): 0.25}, {("D", "x"): 1.0, ("D", "y"): -1.0})}
    s = {"drug": "D", "cell_line": "c1", "values": {"x": 0.5, "y": -0.5}}
    assert sf.sample_pairs(s, preds, "gated", {"x", "y"}) == [(0.25, 0.5), (-1.0, -0.5)]
    s2 = dict(s, cell_line="c2")
    assert sf.sample_pairs(s2, preds, "gated", {"x", "y"}) == [(1.0, 0.5), (-1.0, -0.5)]


def test_collect_needs_every_arm_scorable():
    # the reach arm has only 2 of the 3 antibodies, so the sample is left out of every arm (shared denominators)
    preds = {"fold": {("D", a): 1.0 for a in "xyz"},
             "reach": {("D", "x"): 0.0, ("D", "y"): 0.0},
             "gated": ({}, {("D", a): 1.0 for a in "xyz"})}
    s = {"drug": "D", "cell_line": "c", "time": "4hr", "values": {"x": 0.1, "y": 0.2, "z": 0.3}}
    assert sf.collect([s], preds, set("xyz")) == []
    preds["reach"][("D", "z")] = 0.0
    recs = sf.collect([s], preds, set("xyz"))
    assert len(recs) == 1 and recs[0]["minutes"] == 240


def test_primary_bootstraps_the_portion_unit():
    # two cell lines with opposite correlations: resampling cell lines must reach both extremes
    recs = [{"drug": "D", "cell_line": c, "minutes": None, "corr": {"fold": v},
             "calls": {"fold": []}} for c, v in (("c1", 1.0), ("c2", -1.0))]
    point, ci, n_drugs, n_units = sf.primary(recs, "fold", "cell_line")
    assert point == 0.0 and n_drugs == 1 and n_units == 2
    assert ci[0] == -1.0 and ci[1] == 1.0
    # by drug there is one unit, so every resample is the point estimate
    _, ci_d, _, n_units_d = sf.primary(recs, "fold", "drug")
    assert n_units_d == 1 and ci_d == [0.0, 0.0]


def test_balanced_accuracy_three_class():
    recs = [{"drug": "D", "cell_line": "c", "minutes": None, "corr": {},
             "calls": {"fold": [(-1, -1), (0, 0), (1, 0), (1, 1)]}}]
    assert sf.balanced(recs, "fold") == (1.0 + 1.0 + 0.5) / 3
