"""specs/176: the budget's scoring helpers, the curated enzyme classes and the network reach, on synthetic inputs."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import budget_176 as B  # noqa: E402
import followup_176 as F  # noqa: E402
import routes_176 as R  # noqa: E402


def test_call_and_auroc():
    assert B.call_of(0.3, 0.2) == "up" and B.call_of(-0.3, 0.2) == "down" and B.call_of(0.2, 0.2) == "not"
    assert B.auroc([0.1, 0.9, 0.8, 0.2], [0, 1, 1, 0]) == 1.0
    assert B.auroc([1, 1, 1, 1], [0, 1, 1, 0]) == 0.5            # all tied
    assert B.auroc([0, 1, 1, 1], [0, 0, 1, 1]) == pytest.approx(0.75)
    assert B.auroc([1, 2], [1, 1]) is None


def test_scores():
    cases = ([{"measured": "up", "frozen": "up", "reach": "up"}] * 2 + [{"measured": "down", "frozen": "up", "reach": "up"}]
             + [{"measured": "not", "frozen": "up", "reach": "up"}] + [{"measured": "up", "frozen": "not", "reach": "not"}]
             + [{"measured": "not", "frozen": "not", "reach": "not"}] * 3)
    s = B.scores(cases)
    assert s["changes"] == 4 and s["calls"] == 4
    assert s["recall"] == pytest.approx(2 / 4) and s["precision_exact"] == pytest.approx(2 / 4)
    assert s["precision_any"] == pytest.approx(3 / 4)
    assert s["enrichment"] == pytest.approx((3 / 4) / (1 / 4))


def write(d, name, text):
    (d / name).write_text(text)


def toy(d):
    """K:S -> (catalyst K:S) -> K:pY-S, an autophosphorylation in trans is absent; K phosphorylates S on tyrosine.
    A (isoform reference) autophosphorylates on serine. B sits in the catalyst complex with K."""
    p = d / "R-HSA-1"
    p.mkdir()
    write(p, "logic_network.csv", "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
          "s,r1,pos,and,input,1,R1\nkb,r1,pos,,catalyst,1,R1\nr1,sp,pos,and,output,1,R1\n"
          "a,r2,pos,and,input,1,R2\na,r2,pos,,catalyst,1,R2\nr2,ap,pos,and,output,1,R2\n"
          "sp,r3,pos,and,input,1,R3\nr3,t,pos,and,output,1,R3\n")
    write(p, "nodes.csv", "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
          "s,simple_entity,S,,L-S,,\nkb,simple_complex,KB,,L-K|L-B,,\nr1,reaction,R1,,,,\nsp,simple_entity,SP,,L-SP,,\n"
          "a,simple_entity,A,,L-A,,\nr2,reaction,R2,,,,\nap,simple_entity,AP,,L-AP,,\nr3,reaction,R3,,,,\n"
          "t,simple_entity,T,,L-T,,\n")
    write(p, "leaf_refs.csv", "leaf_stid,reference_stid,modifications\nL-S,uniprot:PS,\nL-K,uniprot:PK,\n"
          "L-B,uniprot:PB,\nL-SP,uniprot:PS,O4'-phospho-L-tyrosine at 10\nL-A,uniprot:PA-2,\n"
          "L-AP,uniprot:PA-2,O-phospho-L-serine at 5\nL-T,uniprot:PT,\n")
    return p


def test_enzyme_classes(tmp_path):
    toy(tmp_path)
    cls = F.enzyme_classes(tmp_path, {"PK", "PB", "PA", "PS"})
    assert cls["PK"] == {"catalyst", "tyr_kinase"} and cls["PB"] == {"catalyst", "tyr_kinase"}   # complex members both
    assert cls["PA"] == {"catalyst", "self_phos", "st_kinase"}                                   # isoform suffix dropped
    assert "PS" not in cls


def test_net_reach_and_sites(tmp_path):
    net = R.Net(toy(tmp_path))
    assert net.by_site[("PS", 10)] == {"sp"} and net.by_site[("PA", 5)] == {"ap"}
    assert net.unmodified["PS"] == {"s"}
    assert "t" in net.reach({"s"}) and "t" not in net.reach({"a"})
