"""specs/050 step 2: the motif census's pure counting functions."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import motif_census as mc  # noqa: E402


def test_largest_scc():
    assert mc.largest_scc([]) == 0
    assert mc.largest_scc([("a", "b"), ("b", "c")]) == 1           # acyclic: singletons
    assert mc.largest_scc([("a", "b"), ("b", "a"), ("b", "c")]) == 2
    # two cycles joined into one component, plus a separate 2-cycle
    edges = [("a", "b"), ("b", "c"), ("c", "a"), ("c", "d"), ("d", "c"), ("x", "y"), ("y", "x")]
    assert mc.largest_scc(edges) == 4
    # deep chain into a cycle: iterative, no recursion limit
    chain = [(str(i), str(i + 1)) for i in range(5000)] + [("5000", "4990")]
    assert mc.largest_scc(chain) == 11


def e(s, rxn):
    return {"source_id": s, "edge_type": "input", "edge_reaction_id": rxn}


def test_competition_sites():
    edges = [e("S", "R-1"), e("S", "R-2"), e("T", "R-1"), e("U", "R-3"), e("U", "R-4")]
    assert mc.competition_sites(edges, set(), {"R-1"}) == 1          # S only: T has one consumer
    assert mc.competition_sites(edges, set(), {"R-1", "R-3"}) == 2   # S and U
    assert mc.competition_sites(edges, {"S"}, {"R-1"}) == 0          # carrier / pool node excluded
    assert mc.competition_sites(edges, set(), set()) == 0            # no GO table: nothing
