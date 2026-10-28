import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import depth_120 as t  # noqa: E402


def row(group, measured, calls):
    return {"group": group, "readout": ("p",), "measured": measured, "model": calls[0.0], "reach": "not",
            "calls": calls}


def test_merge_keeps_cases_scored_in_every_arm():
    a = {("P", "A", "x"): {"group": "A", "readout": ("x",), "measured": "up", "model": "up", "reach": "up"},
         ("P", "B", "x"): {"group": "B", "readout": ("x",), "measured": "not", "model": "up", "reach": "up"}}
    b = {("P", "A", "x"): dict(a[("P", "A", "x")], model="not")}
    rows, dropped = t.merge({0.0: a, 0.5: b})
    assert len(rows) == 1 and rows[0]["calls"] == {0.0: "up", 0.5: "not"}
    assert dropped == {0.0: 1, 0.5: 0}


def test_choose_prefers_the_smaller_depth_on_a_tie():
    train = [row("A", "up", {0.0: "up", 0.5: "up"}), row("B", "not", {0.0: "not", 0.5: "not"})]
    assert t.choose(train, [0.0, 0.5]) == 0.0
    train.append(row("C", "not", {0.0: "down", 0.5: "not"}))
    assert t.choose(train, [0.0, 0.5]) == 0.5


def test_lodo_chooses_the_depth_without_the_held_out_drug():
    # B's cases favour depth 0.5; A's favour 0. Holding out A, B decides, and vice versa.
    rows = [row("A", "up", {0.0: "up", 0.5: "not"}), row("A", "not", {0.0: "not", 0.5: "not"}),
            row("B", "not", {0.0: "down", 0.5: "not"}), row("B", "down", {0.0: "down", 0.5: "down"})]
    preds, chosen = t.lodo(rows, "korkut")
    assert chosen == {"A": 0.5, "B": 0.0}
    assert [p["fitted"] for p in preds if p["drug"] == "A"] == ["not", "not"]
    assert [p["frozen"] for p in preds if p["drug"] == "A"] == ["up", "not"]


def test_false_change():
    assert t.false_change([("up", "not"), ("down", "down"), ("not", "up")]) == 0.5
