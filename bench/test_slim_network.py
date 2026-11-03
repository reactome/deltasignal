import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slim_network as sn  # noqa: E402


def test_slim_keeps_routes_from_sources_to_readouts():
    # T -> a -> R ; b -> a (b not reachable from T: boundary) ; a -> z (z reaches no readout: dropped) ; q -> R2
    E = [("T", "a"), ("a", "R"), ("b", "a"), ("a", "z"), ("q", "R2")]
    parsed = {"nodes": [{"uuid": u, "baseline": 0.01} for u in "T a R b z q R2".split()],
              "edges": [{"parent_uuid": p, "child_uuid": c} for p, c in E]}
    nodes, edges, boundary = sn.slim(parsed, {"T"}, {"R", "R2"})
    assert {n["uuid"] for n in nodes} == {"T", "a", "R", "b"}
    assert {(e["parent_uuid"], e["child_uuid"]) for e in edges} == {("T", "a"), ("a", "R"), ("b", "a")}
    assert boundary == {"b"}
    assert sn.boundary_observations(parsed, boundary) == {"b": [1.0, 1.0]}
