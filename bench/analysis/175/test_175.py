"""specs/175: the control comparison and the sign test, on synthetic rows."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_175 as S  # noqa: E402


def row(pred="up", reach="up", ui="2.0", status="scored"):
    return {"status": status, "predicted": pred, "reach_predicted": reach, "pred_ui": ui}


def test_control_ignores_totals_and_catches_phospho_changes():
    old = {("P", "G", "p-AKT"): row(), ("P", "G", "beta-catenin"): row(), ("P", "H", "p-ERK1/2"): row(ui="3.0")}
    new = {("P", "G", "p-AKT"): row(ui="2.0000000001"), ("P", "G", "beta-catenin"): row(pred="not", ui="1.0"),
           ("P", "H", "p-ERK1/2"): row(ui="3.1"), ("P", "X", "p-JNK"): row()}
    bad = S.control(old, new, {"beta-catenin"})
    assert sorted(bad) == [(("P", "H", "p-ERK1/2"), "pred_ui"), (("P", "X", "p-JNK"), "extra in new")]


def test_sign_test():
    assert S.sign_test(0, 0) == 1.0
    assert abs(S.sign_test(5, 0) - 1 / 32) < 1e-12
    assert S.sign_test(3, 3) > 0.5
