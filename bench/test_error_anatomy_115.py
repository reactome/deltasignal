import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
from error_anatomy_115 import anatomy, band  # noqa: E402


def row(m, p, reach, ds="d", ro="r", status="scored"):
    return {"pathway_id": "P", "drug_set": ds, "readout": ro, "status": status,
            "measured": m, "predicted": p, "reach_predicted": reach}


def test_classes():
    rows = [row("down", "down", "down"), row("down", "up", "down"), row("down", "not", "down"),
            row("up", "not", "not"), row("not", "down", "down"), row("not", "not", "not"),
            row("conflict", "down", "down"), row("down", "down", "down", status="trivial")]
    assert anatomy(rows)["classes"] == {"changed / right": 1, "changed / wrong direction": 1,
                                        "changed / model NOT": 1, "changed / no route": 1,
                                        "unchanged / false change call": 1, "unchanged / right": 1}


def test_bands_and_direction():
    assert [band(x) for x in ("", "1", "5", "6", "20", "21", "41")] == \
        ["no route", "1-5", "1-5", "6-10", "11-20", "21-40", ">40"]
    rows = [row("down", "down", "down", ro="a"), row("down", "up", "down", ro="b"), row("not", "up", "down", ro="c")]
    hops = {("P", "d", "a"): "2", ("P", "d", "b"): "2", ("P", "d", "c"): "30"}
    b = anatomy(rows, hops, ("pathway_id", "drug_set", "readout"))["by_hops"]
    assert b["1-5"] == {"cases": 2, "changed": 2, "model_calls": 2, "model_calls_right": 1,
                        "both_changed": 2, "direction_agrees": 1}
    assert b["21-40"]["model_calls"] == 1 and b["21-40"]["model_calls_right"] == 0
