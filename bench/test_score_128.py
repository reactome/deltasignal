"""specs/128: the korkut comparison in score_128.py on small hand-made tables."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import score_128 as s  # noqa: E402


def row(pred, meas):
    return {"predicted": pred, "measured": meas, "status": "scored"}


def test_korkut_counts_and_mcnemar():
    ctrl = {1: row("up", "not"), 2: row("down", "down"), 3: row("up", "up"), 4: row("not", "not")}
    arm = {1: row("not", "not"), 2: row("down", "down"), 3: row("not", "up"), 4: row("not", "not"), 5: row("up", "up")}
    r = s.korkut(ctrl, arm)
    assert r["cases"] == 4 and r["always_not_right"] == 2
    assert (r["control_right"], r["arm_right"]) == (3, 3)
    assert (r["control_change_calls"], r["arm_change_calls"]) == (3, 1)
    assert (r["fixed"], r["broke"]) == (1, 1) and r["mcnemar_p"] == 1.0


def test_macro_f1():
    assert s.macro_f1([("up", "up"), ("down", "down"), ("not", "not")]) == 1.0
    assert s.macro_f1([("up", "not"), ("not", "not")]) == (0 + 0 + 2 / 3) / 3
