import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import train_117 as t  # noqa: E402


def test_resp_bin():
    assert [t.resp_bin(c, n) for c, n in ((0, 0), (0, 5), (1, 10), (2, 10), (5, 10), (6, 10))] == \
        ["unseen", 0, 1, 2, 3, 4]


def test_responsiveness_counts_a_measurement_once():
    # the same (group, readout) scored in three pathways is one measurement
    rows = [{"group": "A", "readout": ("p",), "measured": "down"}] * 3 + \
        [{"group": g, "readout": ("p",), "measured": "not"} for g in ("B", "C")]
    assert t.responsiveness(rows) == {("p",): (3, "down")}     # 1 of 3 changed: (0.25, 0.5]


def test_fallback_to_shorter_cell_then_not():
    train = [({"bin": 1, "dir": "up", "model": "up"}, "up")] * 10 + [({"bin": 1, "dir": "up", "model": "not"}, "not")] * 5 \
        + [({"bin": 0, "dir": "tie", "model": "not"}, "not")] * 20
    tables = t.fit(train, ("bin", "dir", "model"))
    assert t.predict(tables, {"bin": 1, "dir": "up", "model": "up"}) == "up"
    assert t.predict(tables, {"bin": 1, "dir": "up", "model": "not"}) == "up"    # 5 < MIN_CELL: (bin, dir) cell
    assert t.predict(tables, {"bin": 4, "dir": "down", "model": "down"}) == "not"


def test_lodo_responsiveness_never_sees_the_held_out_drug():
    rows = [{"group": g, "readout": ("p",), "measured": m, "model": "not", "reach": "not"}
            for g, m in (("A", "up"), ("B", "not"), ("A+B", "up"), ("C", "not"))]
    seen = []
    orig = t.responsiveness
    t.responsiveness = lambda tr: seen.append({d for r in tr for d in r["drugs"]}) or orig(tr)
    try:
        preds = t.lodo(rows, "korkut")
    finally:
        t.responsiveness = orig
    assert seen == [{"B", "C"}, {"A", "C"}, {"A", "B"}]
    assert len(preds) == 5            # A+B is tested once for A and once for B
