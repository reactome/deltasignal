"""specs/129: the grid choice and cross-validation in score_129.py on small hand-made data."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_129 as s  # noqa: E402


def test_elasticity():
    assert s.elasticity(math.inf, 1.0) == 1.0
    assert s.elasticity(math.inf, 0.5) == 0.5
    assert abs(s.elasticity(9, 1.0) - 0.9) < 1e-12
    assert abs(s.elasticity(3, 0.5) - 0.5 * 3 ** 0.5 / (1 + 3 ** 0.5)) < 1e-12


def test_read_table(tmp_path):
    p = tmp_path / "cf_K9_h075.tsv"
    p.write_text("# grid\nclass\tK\th\ninput\t9\t0.75\ncatalyst\t9\t0.75\nregulator_pos\t9\t0.75\n")
    assert s.read_table(p) == (9.0, 0.75)
    p.write_text("class\tK\th\ninput\t9\t1\ncatalyst\t3\t1\n")
    try:
        s.read_table(p)
    except ValueError:
        pass
    else:
        raise AssertionError("a table with two curves must be refused")


def test_balanced_accuracy_ignores_absent_classes():
    assert s.balanced_accuracy([("up", "up"), ("not", "up")]) == 0.5
    assert s.balanced_accuracy([("not", "not")] * 9 + [("up", "not")]) == 0.5


def test_folds_are_deterministic_and_balanced():
    a = s.assign_folds(list("abcdefghij") * 3)
    assert a == s.assign_folds(list("jihgfedcba"))
    assert sorted(list(a.values()).count(f) for f in range(5)) == [2] * 5


def _data():
    # Two datasets, ten perturbations each; the identity misses every "up", the contraction gets them all.
    truth, ident, contr, worse = {}, {}, {}, {}
    for ds in ("experimental", "korkut"):
        truth[ds], ident[ds], contr[ds], worse[ds] = {}, {}, {}, {}
        for g in range(10):
            for j, m in enumerate(("up", "not", "down")):
                c = (g, j)
                truth[ds][c] = (f"g{g}", m)
                ident[ds][c] = "not" if m == "up" else m
                contr[ds][c] = m
                worse[ds][c] = "not"
    return truth, {"cf_off": ident, "cf_K9_h1": contr, "cf_K3_h05": worse}


POINTS = {"cf_off": (math.inf, 1.0), "cf_K9_h1": (9.0, 1.0), "cf_K3_h05": (3.0, 0.5)}


def test_choose_and_tie_break():
    truth, preds = _data()
    assert s.choose(POINTS, preds, truth, lambda ds, g: True) == "cf_K9_h1"
    tied = dict(preds, cf_K9_h1=preds["cf_off"])
    assert s.choose(POINTS, tied, truth, lambda ds, g: True) == "cf_off"


def test_cross_validation_and_gain():
    truth, preds = _data()
    cv, chosen, full = s.cross_validate(POINTS, preds, truth)
    assert chosen == ["cf_K9_h1"] * 5 and full == "cf_K9_h1"
    assert cv == preds["cf_K9_h1"]
    point, low, high = s.bootstrap_gain(truth, cv, preds["cf_off"], n=200)
    assert abs(point - 1 / 3) < 1e-12 and abs(low - 1 / 3) < 1e-12 and abs(high - 1 / 3) < 1e-12
    point, low, high = s.bootstrap_gain(truth, preds["cf_off"], preds["cf_off"], n=200)
    assert point == low == high == 0.0


def test_bootstrap_point_matches_balanced_accuracy():
    truth = {"d": {i: (f"g{i % 4}", ("up", "not", "down")[i % 3]) for i in range(24)}}
    a = {"d": {i: ("up", "not")[i % 2] for i in range(24)}}
    b = {"d": {i: "not" for i in range(24)}}
    point, low, high = s.bootstrap_gain(truth, a, b, n=500)
    exp = (s.balanced_accuracy([(m, a["d"][c]) for c, (_, m) in truth["d"].items()])
           - s.balanced_accuracy([(m, "not") for _, m in truth["d"].values()]))
    assert abs(point - exp) < 1e-12 and low <= point <= high
