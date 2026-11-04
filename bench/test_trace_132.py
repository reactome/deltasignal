import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import trace_132 as tr  # noqa: E402


def test_on_route_keeps_bridges_the_sources_reach_and_that_reach_the_readout():
    fwd = {"g": ["q1"], "q1": ["p1"], "q2": ["p2"], "p1": ["r"], "x": ["p3"]}
    br = [("q1", "p1"), ("q2", "p2"), ("x", "p3")]
    # q2 -> p2 is not reached from g; x -> p3 does not reach r
    assert tr.on_route(br, fwd, tr.reverse(fwd), ["g"], ["r"]) == [("q1", "p1")]


def test_name_bridge_matches_the_head_entity_or_a_set_holding_it():
    census = [{"s": "S1", "r": "R1", "q": "Q", "entity": "E", "role": "input"},
              {"s": "S2", "r": "R2", "q": "Q", "entity": "SET", "role": "catalyst"},
              {"s": "S3", "r": "R3", "q": "Q", "entity": "OTHER", "role": "input"}]
    stids_of = {"p1": {"E"}, "p2": {"M"}}
    members = {"SET": {"M", "N"}}
    assert tr.name_bridge("p1", census, stids_of, members) == [("S1", "R1", "Q", "E", "input")]
    assert tr.name_bridge("p2", census, stids_of, members) == [("S2", "R2", "Q", "SET", "catalyst")]
    assert tr.name_bridge("p9", census, stids_of, members) == []
