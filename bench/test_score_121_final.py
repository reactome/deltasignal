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


def _drug_records():
    # four drugs, unequal sample counts and spread, so a regrouping bug changes the interval
    vals = {"D1": [0.9, 0.7], "D2": [-0.4], "D3": [0.1, 0.2, 0.3], "D4": [-0.8, -0.6]}
    return [{"drug": d, "cell_line": f"c{i}", "minutes": None, "corr": {"fold": v}, "calls": {"fold": []}}
            for d, vs in vals.items() for i, v in enumerate(vs)]


def test_drug_bootstrap_matches_direct_reference():
    # review 2026-11-20: resampling drugs must weight a drug by how many times it was drawn. The reference draws
    # drug means directly with the same generator: the textbook bootstrap over drugs.
    import random
    recs = _drug_records()
    point, ci, n_drugs, n_units = sf.primary(recs, "fold", "drug")
    means = sf.per_drug(recs, "fold")
    drugs = sorted(means)
    rng = random.Random(sf.SEED)
    draws = sorted(sum(means[rng.choice(drugs)] for _ in drugs) / len(drugs) for _ in range(sf.N_BOOT))
    ref = [draws[int(0.025 * sf.N_BOOT)], draws[int(0.975 * sf.N_BOOT) - 1]]
    assert n_drugs == n_units == 4
    assert math.isclose(point, sum(means.values()) / 4)
    assert all(math.isclose(a, b, abs_tol=1e-12) for a, b in zip(ci, ref))
    assert math.isclose(ci[1] - ci[0], ref[1] - ref[0], abs_tol=1e-12)
    assert ci[1] - ci[0] > 0.5        # four drugs this far apart cannot give a narrow interval


def test_drug_drawn_twice_counts_twice(monkeypatch):
    # every resample draws D1, D1, D2: the mean over drugs is (m1 + m1 + m2) / 3, not (m1 + m2) / 2
    recs = [r for r in _drug_records() if r["drug"] in ("D1", "D2")] + \
           [{"drug": "D3", "cell_line": "c0", "minutes": None, "corr": {"fold": 0.0}, "calls": {"fold": []}}]
    picks = iter(["D1", "D1", "D2"] * sf.N_BOOT)

    class Fixed:
        def __init__(self, seed):
            pass

        def choice(self, seq):
            return next(picks)

    monkeypatch.setattr(sf.random, "Random", Fixed)
    _, ci, _, _ = sf.primary(recs, "fold", "drug")
    expect = (0.8 + 0.8 - 0.4) / 3
    assert math.isclose(ci[0], expect) and math.isclose(ci[1], expect)
