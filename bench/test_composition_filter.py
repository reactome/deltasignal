"""Pins for cycle_closing_composition_pairs (specs/029)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import cycle_closing_composition_pairs  # noqa: E402


def _write(d, rows):
    (d / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type\n"
        + "".join(f"{s},{t},pos,and,{et}\n" for s, t, et in rows))


def test_a_bridging_edge_is_kept_and_a_closing_edge_dropped(tmp_path):
    # C -> R -> P ; composition A -> C bridges (acyclic);
    # composition P -> C would close C -> R -> P -> C.
    _write(tmp_path, [("C", "R", "input"), ("R", "P", "output"),
                      ("A", "C", "composition"), ("P", "C", "composition")])
    assert cycle_closing_composition_pairs(tmp_path) == {("P", "C")}


def test_two_edges_that_close_a_cycle_together_are_caught(tmp_path):
    # Each alone is acyclic; together X -> Y -> X. The second in sorted order is dropped.
    _write(tmp_path, [("X", "Y", "composition"), ("Y", "X", "composition")])
    assert cycle_closing_composition_pairs(tmp_path) == {("Y", "X")}


def test_no_composition_edges_nothing_dropped(tmp_path):
    _write(tmp_path, [("A", "B", "input"), ("B", "A", "output")])
    assert cycle_closing_composition_pairs(tmp_path) == set()
