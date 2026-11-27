"""specs/179: engagement, the placebo, the gate and the direct-pair table, on synthetic inputs."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import direct_179 as D  # noqa: E402
import engage_179 as E  # noqa: E402


def test_engagement_and_placebo():
    m = {("K", "s1"): (0.5, 0.2), ("K", "s2"): (0.0, 0.2), ("L", "s1"): (-0.5, 0.2), ("L", "x"): (0.5, 0.2),
         ("M", "x"): (0.5, 0.2), ("K", "x"): (0.0, 0.2)}
    direct = {("K", "s1"), ("K", "s2"), ("L", "s1"), ("N", "s1")}
    assert E.engagement(m, direct) == {"K": True, "L": False}          # L's substrate fell; M and N have none measured
    reach = {("K", "x"): "not", ("L", "x"): "not", ("M", "x"): "up", ("K", "s1"): "not"}
    assert E.placebo(m, reach, direct) == {"K": False, "L": True}      # M's readout has a route; direct pairs left out


def test_gate_counts():
    cases = ([{"gene": "a", "engaged": False, "model": "up", "reach": "up", "measured": "not"}] * 6
             + [{"gene": "b", "engaged": False, "model": "up", "reach": "up", "measured": "up"}]
             + [{"gene": "c", "engaged": True, "model": "up", "reach": "up", "measured": "up"}] * 200)
    g = E.gate([dict(c) for c in cases], "model")
    assert (g["withheld"], g["fixed"], g["broken"]) == (7, 6, 1)
    assert g["recall"] == 1.0 and abs(g["recall_gated"] - 200 / 201) < 1e-12
    assert g["passes"] is False                                         # 6 vs 1 is p = 0.0625
    g = E.gate([dict(c) for c in cases[:6] + cases[7:]], "model")
    assert g["passes"] is True                                          # 6 vs 0, recall untouched


def test_direct_pairs(tmp_path):
    p = tmp_path / "R-HSA-1"
    p.mkdir()
    (p / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
        "s,r1,pos,and,input,1,R1\nk,r1,pos,,catalyst,1,R1\nr1,sp,pos,and,output,1,R1\n"
        "sp,r2,pos,and,input,1,R2\nq,r2,pos,,catalyst,1,R2\nr2,spc,pos,and,output,1,R2\n")
    (p / "nodes.csv").write_text(
        "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
        "s,simple_entity,S,,L-S,,\nk,simple_entity,K,,L-K,,\nr1,reaction,R1,,,,\nsp,simple_entity,SP,,L-SP,,\n"
        "q,simple_entity,Q,,L-Q,,\nr2,reaction,R2,,,,\nspc,simple_complex,SPC,,L-SP|L-Q,,\n")
    (p / "leaf_refs.csv").write_text(
        "leaf_stid,reference_stid,modifications\nL-S,uniprot:PS,\nL-K,uniprot:PK,\nL-Q,uniprot:PQ,\n"
        "L-SP,uniprot:PS,O-phospho-L-serine at 9\n")
    site, own = {("PS", 9): "p-S"}, {"p-S": {"S"}}
    pairs = D.direct_pairs(tmp_path, {"PK": "K", "PQ": "Q", "PS": "S"}, site, own)
    assert pairs == {("K", "p-S"): ["R-HSA-1"]}       # Q only binds the already-phosphorylated form
