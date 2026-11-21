"""specs/166 follow-up: route attribution of change calls to literature edges."""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis" / "166"))
import attribute_166 as a  # noqa: E402


def graph(edges):
    adj, radj = defaultdict(list), defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        radj[v].append(u)
    return adj, radj


def test_an_edge_is_on_route_only_between_source_and_readout():
    # G -> L1 -> R (L1 a literature edge into the route); X -L2-> R (source not reached); G -L3-> Y (dead end)
    adj, radj = graph([("G", "a"), ("a", "b"), ("b", "R"), ("X", "R"), ("G", "Y")])
    lit = [("a", "b", "L1"), ("X", "R", "L2"), ("G", "Y", "L3")]
    assert a.on_route(lit, adj, radj, ["G"], ["R"]) == {"L1"}


def test_classify_new_and_kept():
    assert a.classify("not", "not", "down") is None
    assert a.classify("down", "not", "not") == "new_wrong"
    assert a.classify("down", "not", "down") == "new_right"
    assert a.classify("up", "down", "down") == "kept_wrong"
    assert a.classify("down", "up", "down") == "kept_right"
