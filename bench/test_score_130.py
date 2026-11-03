"""specs/130: the readout-baseline transform and the drug cross-validation in score_130.py, on hand-made data."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_130 as s  # noqa: E402


def test_measured_change_limits():
    # p0 far above background: w -> 1, the prediction is the frozen fold
    assert abs(s.measured_change(3.0, 20.0, 0.0) - 3.0) < 1e-5
    # at background (w = 1/2): a knockout halves the signal, a 4x rise reads 2.5x
    assert s.measured_change(-math.inf, 0.0, 0.0) == -1.0
    assert abs(s.measured_change(2.0, 0.0, 0.0) - math.log2(2.5)) < 1e-12
    # no change stays no change at any baseline
    assert s.measured_change(0.0, -3.0, 1.0) == 0.0
    # the lower the resting level, the smaller the predicted fall
    assert s.measured_change(-2.0, -1.0, 0.0) > s.measured_change(-2.0, 1.0, 0.0)


def test_quantile_and_backgrounds():
    assert s.quantile([3, 1, 2], 0.5) == 2
    assert s.quantile([0, 10], 0.25) == 2.5
    bl = {("set1", "a"): {"X": 0.0, "Y": 1.0}, ("set1", "b"): {"X": 2.0}, ("set2", "a"): {"X": 4.0}}
    assert s.backgrounds(bl, 0.5) == {"X": 2.0, "Y": 1.0}


def test_sample_predictions_identity_and_missing_baseline():
    sample = {"set": "set1", "cell_line": "a", "drug": "d", "values": {"X": 0.1, "Y": 0.2, "Z": 0.3}}
    fold = {("d", "X"): -2.0, ("d", "Y"): 1.0}
    bl = {("set1", "a"): {"X": 0.0}}
    assert s.sample_predictions(sample, fold, bl, None) == {"X": -2.0, "Y": 1.0}
    p = s.sample_predictions(sample, fold, bl, {"X": 0.0, "Y": 0.0})
    assert p["Y"] == 1.0                                   # no baseline for Y in this cell: kept
    assert abs(p["X"] - math.log2(0.5 * 0.25 + 0.5)) < 1e-12


def test_classes():
    t = s.THRESHOLD
    assert (s.cls(t), s.cls(-t), s.cls(t * 0.99)) == (1, -1, 0)


def test_choose_cross_validate_and_gain():
    drugs = [f"d{i}" for i in range(10)]
    pd = {None: {d: 0.10 for d in drugs}, 0.25: {d: 0.20 for d in drugs}, 0.5: {d: 0.05 for d in drugs}}
    cv, chosen, full = s.cross_validate(pd, drugs)
    assert chosen == [0.25] * 5 and full == 0.25 and cv == pd[0.25]
    point, low, high = s.bootstrap_gain(cv, pd[None], drugs, n=200)
    assert abs(point - 0.1) < 1e-12 and abs(low - 0.1) < 1e-12 and abs(high - 0.1) < 1e-12
    tie = {None: {d: 0.1 for d in drugs}, 0.25: {d: 0.1 for d in drugs}}
    assert s.choose(tie, drugs) is None


def test_balanced_accuracy():
    scores = [("d1", 0.5, [(1, 1), (0, 0), (0, 1)]), ("d2", 0.1, [(-1, 0)])]
    assert s.balanced_accuracy(scores, {"d1"}) == (1 + 0.5) / 2
    assert s.balanced_accuracy(scores, {"d1", "d2"}) == (1 + 0.5 + 0) / 3
