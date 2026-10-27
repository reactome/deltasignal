"""specs/114 part 2: summary_114.py figures on a hand-made case log."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import summary_114 as s  # noqa: E402


def row(drug, readout, pid, measured, pred, reach, status="scored"):
    return {"drug_set": drug, "readout": readout, "pathway_id": pid, "status": status,
            "measured": measured, "predicted": pred, "reach_predicted": reach}


def test_korkut_summary():
    rows = [row("d1", "A", "P1", "down", "down", "down"),
            row("d1", "A", "P2", "down", "not", "down"),       # same unit, other network, other call
            row("d1", "B", "P1", "not", "not", "down"),
            row("d2", "B", "P1", "up", "up", "down"),
            row("d2", "C", "P1", "conflict", "down", "down"),   # excluded
            row("d2", "D", "P1", "down", "", "", status="trivial")]
    out = s.summarize("korkut", rows)
    assert out["scored"] == 4
    assert out["accuracy"] == {"model": 0.75, "reach": 0.5, "always_not": 0.25}
    assert out["regulated"] == {"n": 3, "model": 2, "reach": 2}
    assert out["mcnemar"]["model_only"] == 2 and out["mcnemar"]["reach_only"] == 1
    assert out["units"] == 3 and out["units_with_disagreeing_networks"] == 1
