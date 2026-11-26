"""specs/174: the registered case classes, readout types, groups and the network tracing, on synthetic inputs."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import k562_174 as K  # noqa: E402
import trace_174 as T  # noqa: E402


def test_case_classes():
    assert T.case_class("up", "up", "up") == "R1"
    assert T.case_class("up", "up", "down") == "R2"
    assert T.case_class("up", "down", "not") == "R3"
    assert T.case_class("up", "not", "up") == "R4"
    assert T.case_class("up", "conflict", "not") == "R5"      # a conflict counts as no call
    assert T.case_class("not", "up", "up") == "N1"
    assert T.case_class("not", "not", "not") == "N2"


def cases(readout, on, on_changed, off, off_changed):
    return ([{"readout": readout, "reach": "up", "measured": "up" if i < on_changed else "not"} for i in range(on)]
            + [{"readout": readout, "reach": "not", "measured": "up" if i < off_changed else "not"} for i in range(off)])


def test_readout_types():
    t = T.readout_types(cases("bg", 10, 5, 20, 5) + cases("rs", 10, 4, 20, 2) + cases("few", 9, 9, 20, 0)
                        + cases("low", 10, 1, 20, 0))
    assert t["bg"]["type"] == "background"          # 5/20 = 0.25 without a route
    assert t["rs"]["type"] == "route-specific"      # 0.4 vs 0.1
    assert t["few"]["type"] == "quiet"              # under 10 pairs with a route
    assert t["low"]["type"] == "quiet"              # 0.1 with a route is under the 0.2 floor
    assert T.readout_types(cases("ratio", 10, 3, 20, 4))["ratio"]["type"] == "quiet"   # 0.3 < 2 x 0.2


def test_groups_cover_half_then_stop():
    cs = ([{"cls": "R3", "readout": "a", "rtype": "quiet"}] * 6 + [{"cls": "R3", "readout": "b", "rtype": "quiet"}] * 3
          + [{"cls": "R3", "readout": "c", "rtype": "quiet"}] * 1
          + [{"cls": "R1", "readout": "x", "rtype": "background"}] * 5)
    g = T.groups(cs)
    assert g["R3"]["groups"] == [("a", 6)] and g["R3"]["covered"] == 6
    assert g["R1"]["cases"] == 0                     # R1 groups only on route-specific readouts


def write(d, name, text):
    (d / name).write_text(text)


def toy_pathway(d):
    """P1 (protein A) -> r1 -> complex A:B -| r2 -> C-phospho; A is a root."""
    write(d, "logic_network.csv", "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
          "a,r1,pos,and,input,1,R1\nr1,ab,pos,and,output,1,R1\nab,r2,neg,,regulator,1,R2\n"
          "c0,r2,pos,and,input,1,R2\nr2,cp,pos,and,output,1,R2\n")
    write(d, "nodes.csv", "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
          "a,simple_entity,S-A,,L-A,,\nr1,reaction,R1,,,,\nab,simple_complex,S-AB,,L-A|L-B,,\n"
          "r2,reaction,R2,,,,\nc0,simple_entity,S-C,,L-C,,\ncp,simple_entity,S-CP,,L-CP,,\n")
    write(d, "leaf_refs.csv", "leaf_stid,reference_stid,modifications\nL-A,uniprot:PA,\nL-B,uniprot:PB,\n"
          "L-C,uniprot:PC,\nL-CP,uniprot:PC,O-phospho-L-serine at 15\n")


def test_trace_signed_path(tmp_path):
    toy_pathway(tmp_path)
    net = T.Net(tmp_path)
    assert net.nodes_with("PC", site=15) == ["cp"] and net.nodes_with("PC", unmodified=True) == ["c0"]
    assert sorted(net.nodes_with("PA")) == ["a", "ab"]
    start, steps = net.path(["a"], ["cp"])
    assert start == "a" and [s[0] for s in steps] == ["r1", "ab", "r2", "cp"]
    assert steps[-1][1] == -1                        # one negative edge on the way


def test_k562_regulator_pairs(tmp_path):
    p = tmp_path / "R-HSA-1"
    p.mkdir()
    write(p, "logic_network.csv", "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
          "gene,rx,pos,and,input,1,RX\nrx,prot,pos,and,output,1,RX\ntf,rx,pos,,regulator,1,RX\n"
          "x,ry,pos,and,input,1,RY\nry,y,pos,and,output,1,RY\ntf,ry,pos,,regulator,1,RY\n")
    write(p, "nodes.csv", "uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
          "gene,simple_entity,G,,L-G,,\nrx,reaction,RX,,,,\nprot,simple_entity,P,,L-T,,\n"
          "tf,simple_complex,TF,,L-TF|L-CO,,\nx,simple_entity,X,,L-X,,\nry,reaction,RY,,,,\ny,simple_entity,Y,,L-Y,,\n")
    write(p, "leaf_refs.csv", "leaf_stid,reference_stid,modifications\nL-G,ensembl:ENSG1,\nL-T,uniprot:PT,\n"
          "L-TF,uniprot:PTF,\nL-CO,uniprot:PCO,\nL-X,uniprot:PX,\nL-Y,uniprot:PY,\n")
    # RY has no gene or transcript leaf: not an expression reaction, so TF -> Y is not a pair
    assert K.regulator_pairs(tmp_path) == {("PTF", "PT"), ("PCO", "PT")}
