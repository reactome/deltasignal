"""specs/052: the report's Python mirror of the rule, on the Julia suite's fixture."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import zero_supply_report as zr  # noqa: E402

KIND = {"r1": "reaction", "r2": "reaction", "r3": "reaction", "r4": "reaction", "rs": "reaction",
        "rk": "reaction", "C1": "simple_complex", "C2": "simple_complex", "P": "variant_pool"}
LEAVES = {"E": {"LE"}, "X": {"LX"}, "Y": {"LY"}, "C1": {"LE", "LX"}, "C2": {"LE", "LX", "LY"},
          "T": {"LT"}, "M": {"LM"}, "E2": {"LE"}, "K": {"LE"}, "Z": {"LZ"}, "P": {"LE", "LQ"}}
EDGES = [("E", "r1", "input"), ("X", "r1", "input"), ("r1", "C1", "output"),
         ("C1", "r2", "input"), ("Y", "r2", "input"), ("r2", "C2", "output"),
         ("C2", "r3", "input"), ("r3", "C1", "output"), ("C2", "r4", "input"), ("r4", "T", "output")]


def nodes(ids):
    return {u: (KIND.get(u, "simple_entity"), LEAVES.get(u, set())) for u in ids}


BASE = ["E", "X", "Y", "C1", "C2", "T", "r1", "r2", "r3", "r4"]


def test_mirror_matches_the_julia_fixture():
    assert zr.held_by_rule(nodes(BASE), EDGES, {"E": 0.0, "Y": 10.0}) == {"C1", "C2"}
    # synthesis elsewhere keeps supply
    syn = EDGES + [("M", "rs", "input"), ("rs", "E2", "output"), ("E2", "r1", "input")]
    assert zr.held_by_rule(nodes(BASE + ["M", "rs", "E2"]), syn, {"E": 0.0}) == set()
    # knockdown, overexpression, no zero pin: nothing
    for pins in ({"E": 0.5}, {"E": 2.0}, {"Y": 10.0}):
        assert zr.held_by_rule(nodes(BASE), EDGES, pins) == set()
    # a node pinned above 0 supplies what it carries
    assert zr.held_by_rule(nodes(BASE), EDGES, {"E": 0.0, "C1": 1.0}) == set()
    # catalysts do not carry; pools are not held
    cat = EDGES + [("K", "rk", "catalyst"), ("X", "rk", "input"), ("rk", "Z", "output"),
                   ("C1", "P", "variant_pool")]
    assert zr.held_by_rule(nodes(BASE + ["K", "rk", "Z", "P"]), cat, {"E": 0.0, "K": 0.0}) == {"C1", "C2"}


def test_compare_counts_and_by_construction():
    def row(p, e):
        return {"predicted": p, "expected": e}
    ctrl = {("P", "G", "0", "1"): row("2", "0"), ("P", "G", "0", "2"): row("0", "0"),
            ("Q", "H", "2", "3"): row("1", "2")}
    arm = {("P", "G", "0", "1"): row("0", "0"), ("P", "G", "0", "2"): row("1", "0"),
           ("Q", "H", "2", "3"): row("2", "2")}
    m = zr.compare(ctrl, arm, sorted(ctrl), bc={("P", "G", "0", "1")})
    assert (m["fixed"], m["broken"], m["net"]) == (2, 1, 1)
    assert m["by_construction"] == {"fixed": 1, "broken": 0} and m["net_without_by_construction"] == 0
    assert (m["ko_up_control"], m["ko_up_wrong_fixed"], m["ko_up_right_broken"]) == (1, 1, 0)
    assert m["pathways_moved"] == 2 and m["distinct_readouts"] == 3
