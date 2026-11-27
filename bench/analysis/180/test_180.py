import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lun_report as R  # noqa: E402


def test_split_of():
    eng = {"SRC": True, "CSK": False}
    assert [R.split_of(g, eng) for g in ("SRC", "CSK", "XYZ")] == list(R.SPLITS)


def test_report_leaves_direct_pairs_out_of_split_rows_only():
    # (E, delta): E > delta reads UP
    m = {("SRC", "p-FAK"): (2.0, 1.0), ("SRC", "p-ERK1/2"): (2.0, 1.0),
         ("CSK", "p-SRC"): (0.0, 1.0), ("CSK", "p-ERK1/2"): (0.0, 1.0), ("XYZ", "p-ERK1/2"): (0.0, 1.0)}
    pred = {k: ("up", "up") for k in m}
    direct = {("SRC", "p-FAK"), ("CSK", "p-SRC")}
    out = R.report(m, pred, direct, tk={"SRC"}, tyr={"p-FAK", "p-ERK1/2"})
    assert out["all (raw)"]["pairs"] == 5
    assert [out[s]["pairs"] for s in R.SPLITS] == [1, 1, 1]
    assert out["engaged"]["model"]["precision_exact"] == 1.0
    assert out["not engaged"]["model"]["precision_exact"] == 0.0


def test_empty_split_is_none():
    m = {("XYZ", "p-ERK1/2"): (0.0, 1.0)}
    out = R.report(m, {k: ("not", "not") for k in m}, set(), set(), set())
    assert out["engaged"] is None and out["no measured direct pair"]["pairs"] == 1
