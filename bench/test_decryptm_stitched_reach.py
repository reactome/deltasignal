"""specs/051: the stitched-reachability diagnostic on a two-pathway fixture."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
import decryptm_stitched_reach as sr  # noqa: E402

NODES = "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
NET = "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"


def _pathway(d, nodes, edges):
    d.mkdir(parents=True)
    (d / "nodes.csv").write_text(NODES + "".join(f"{u},simple_entity,{st},,{st},,\n" for u, st in nodes))
    (d / "logic_network.csv").write_text(NET + "".join(f"{a},{b},pos,and,input,1,\n" for a, b in edges))


def test_boundary_gain_and_split(tmp_path, capsys):
    b = tmp_path / "build"
    # P1: target T -> hub H. P2: hub H -> site S (different uuids, same identity for H)
    _pathway(b / "R-HSA-1", [("t", "R-T"), ("h", "R-H")], [("t", "h")])
    _pathway(b / "R-HSA-2", [("h2", "R-H"), ("s", "R-S"), ("t2", "R-T"), ("x", "R-X")], [("h2", "s")])
    (tmp_path / "t.csv").write_text('"R-T", "TGT"\n')
    (tmp_path / "p.csv").write_text('"R-S", "P1", "G", "5", "O-phospho-L-serine"\n')
    cases = tmp_path / "c.tsv"
    head = "status\tsplit\tmeasured\ttargets\tuniprot\tposition\tpathway_id\n"
    cases.write_text(head + "scored\theld-out\tdown\tTGT\tP1\t5\tR-HSA-2\n"
                     "scored\theld-out\tnot\tTGT\tP1\t5\tR-HSA-1\n"
                     "trivial\theld-out\tdown\tTGT\tP1\t5\tR-HSA-2\n")
    t = sr.main(["x", str(b), str(cases), str(tmp_path / "t.csv"), str(tmp_path / "p.csv")])
    # in R-HSA-2 the target copy t2 has no edge to s; stitched, t -> H -> s
    assert t[("held-out", "regulated", False, True)] == 1
    # in R-HSA-1 the site is not in the pathway at all; stitched it is reachable
    assert t[("held-out", "unchanged", False, True)] == 1
    assert sum(t.values()) == 2               # the trivial case is not counted
    out = capsys.readouterr().out
    assert "only when stitched 1" in out
    # keyed own pathway: R-HSA-2's own copy of T has no edge, so merging copies inside it adds nothing
    assert "keyed own pathway     0" in out
