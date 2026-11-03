import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drug_activity as da  # noqa: E402

NODES = ("uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
         "mek,simple_entity,MEK,,LMEK,,\nrafmek,simple_complex,RAFMEK,,LRAF|LMEK,,\nmekerk,simple_complex,MEKERK,,LMEK|LERK,,\n"
         "r1,reaction,R1,,,,\nr2,reaction,R2,,,,\nr3,reaction,R3,,,,\nerk,simple_entity,ERK,,LERK,,\n")
EDGES = ("source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
         "mek,r1,pos,and,catalyst,1,R1\nrafmek,r2,pos,and,catalyst,1,R2\nmekerk,r3,pos,and,catalyst,1,R3\n"
         "mek,r3,pos,and,input,1,R3\nerk,r1,pos,and,input,1,R1\n")


def _net(tmp_path):
    (tmp_path / "nodes.csv").write_text(NODES)
    (tmp_path / "logic_network.csv").write_text(EDGES)


def test_cut_respects_active_units(tmp_path):
    _net(tmp_path)
    units = {("R1", "MEK"): None,                      # no active unit: the holder's activity is cut
             ("R2", "RAFMEK"): {"LRAF"},               # RAF is the active unit: a MEK drug leaves it
             ("R3", "MEKERK"): {"LMEK"}}               # MEK is the active unit
    assert da.cut_edges(tmp_path, {"LMEK"}, units) == {("mek", "r1"), ("mekerk", "r3")}   # inputs never cut


def test_redirect_keeps_the_edge_and_holds_its_source_at_zero():
    parsed = {"nodes": [{"uuid": "mek", "baseline": 0.01, "entity_type": "x"}, {"uuid": "r1", "baseline": 0.01}],
              "edges": [{"parent_uuid": "mek", "child_uuid": "r1", "edge_type": "catalyst"},
                        {"parent_uuid": "mek", "child_uuid": "r9", "edge_type": "input"}],
              "pathways": []}
    payload, obs = da.redirected(parsed, {("mek", "r1")}, {})
    assert [e["parent_uuid"] for e in payload["edges"]] == ["drugblock:mek:r1", "mek"]
    assert obs == {"drugblock:mek:r1": [0.0, 1.0]}
    assert len(payload["nodes"]) == 3 and parsed["edges"][0]["parent_uuid"] == "mek"     # the original is untouched
