import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_123  # noqa: E402


def test_paired_difference_and_counts():
    dc = {"ko": {"a": 0.0, "b": 0.1, "c": 0.2}, "mech": {"a": 0.1, "b": 0.1, "c": 0.5}}
    m, ci = score_123.paired(dc, "ko", "mech", ["a", "b", "c"])
    assert abs(m - (0.1 + 0.0 + 0.3) / 3) < 1e-12 and ci[0] <= m <= ci[1] and ci[0] >= 0.0
    out = score_123.report(dc, ["ko", "mech"], ["a", "b", "c"])
    assert out["mech-ko"]["better"] == 2 and out["mech-ko"]["worse"] == 0
    assert abs(out["mean"]["ko"] - 0.1) < 1e-12


def test_identical_arms_give_zero():
    dc = {"x": {"a": 0.3, "b": -0.1}, "y": {"a": 0.3, "b": -0.1}}
    assert score_123.paired(dc, "x", "y", ["a", "b"]) == (0.0, [0.0, 0.0])
