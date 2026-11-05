"""specs/138 amendment 1: korkut reported by target engagement."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "138"))
import engagement_report as er  # noqa: E402

ENG = {"A": "yes", "B": "no", "C": "no marker"}
MARK = {"A": [("Ap", -1)], "B": [("Bp", -1)], "C": []}


def row(ds, ro, pred, meas, reach="not"):
    return ("P", ds, ro), {"status": "scored", "predicted": pred, "measured": meas, "reach_predicted": reach}


def test_split_of():
    assert er.split_of("A", ENG) == "every drug engaged"
    assert er.split_of("A+B", ENG) == "some engaged"
    assert er.split_of("B+C", ENG) == "none engaged"


def test_markers_left_out_of_splits_but_not_raw():
    cases = dict([row("A", "Ap", "down", "down"), row("A", "X", "down", "not"),
                  row("B", "Bp", "down", "not"), row("A+B", "X", "not", "not")])
    out = er.report(cases, ENG, MARK)
    assert out["all (raw)"]["model"]["n"] == 4
    assert out["every drug engaged"]["model"]["n"] == 1        # Ap left out
    assert out["none engaged"]["model"]["n"] == 0              # Bp left out
    assert out["some engaged"]["model"]["n"] == 1
    assert out["every drug engaged"]["always_not"]["accuracy"] == 1.0
