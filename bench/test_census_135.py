"""specs/135: the census's pure parts."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import census_135 as c  # noqa: E402


def _edge(s, t, kind):
    return {"source_id": s, "target_id": t, "edge_type": kind}


def test_copy_producers_starting_made_and_depletion_only():
    rows = [_edge("g", "rx", "input"), _edge("rx", "P", "output"), _edge("D", "S", "depletion"),
            _edge("X", "C", "assembly")]
    assert c.copy_producers(["P", "S", "C", "Z"], rows) == {
        "P": ("made", ["rx"]), "S": ("starting", []), "C": ("other", []), "Z": ("starting", [])}


def test_stock_fraction_and_what_if():
    assert math.isclose(c.stock_fraction(46, 46), 0.5)
    f = c.stock_fraction(24, 46)
    assert math.isclose(f, 1 - 2 ** (-24 / 46)) and 0.30 < f < 0.31
    assert c.what_if(0.0) == 0.0 and c.what_if(0.0, basal=True) == 0.5
    assert math.isclose(c.what_if(0.0, f=f), 1 - f)
    assert math.isclose(c.what_if(0.0, basal=True, f=f), 1 - f / 2)       # 0.848: just under the 0.85 cutoff
    assert c.call(1 - f / 2) == "down" and c.call(0.85) == "not"
