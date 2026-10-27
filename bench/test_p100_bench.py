"""specs/113: the reachability arms on a hand-built network and interaction table."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p100_bench as p  # noqa: E402


def _net(root):
    d = root / "R-HSA-1"
    d.mkdir()
    (d / "leaf_refs.csv").write_text("leaf_stid,reference_stid,modifications\n"
                                     "L1,uniprot:P1,\nL2,uniprot:P2-2,\nL3,uniprot:P3,\nLc,chebi:1,\n")
    # n1 (P1) -> r -> n2 (P2:P3 complex); n3 (P3) is a separate copy reaching nothing
    (d / "nodes.csv").write_text("uuid,node_kind,diagram_entity_id,compartment,member_leaves,source_sets,chosen_members\n"
                                 "n1,entity,,,L1,,\nr,reaction,,,,,\nn2,entity,,,L2|L3|Lc,,\nn3,entity,,,L3,,\n")
    (d / "logic_network.csv").write_text("source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
                                         "n1,r,pos,and,input,1,R\nr,n2,pos,and,output,1,R\n")


def test_network_reach(tmp_path):
    _net(tmp_path)
    reach = p.network_reach(tmp_path)
    assert reach["P1"] == {"P1", "P2", "P3"}          # length 0 included; isoform suffix dropped
    assert reach["P3"] == {"P2", "P3"}                # the complex holding P3 reaches itself
    assert "P9" not in reach


def test_partners_human_only(tmp_path):
    f = tmp_path / "i.csv"
    f.write_text('i.dbId, a.identifier, s.taxId\n1, "P1", "9606"\n1, "Q7", "9606"\n'
                 '2, "P1", "9606"\n2, "W1", "6239"\n3, "Q7", "9606"\n3, "P2-1", "9606"\n')
    assert p.load_partners(f) == {"P1": {"Q7"}, "Q7": {"P1", "P2"}, "P2": {"Q7"}}


def test_arms(tmp_path):
    _net(tmp_path)
    reach = p.network_reach(tmp_path)
    partners = {"X": {"P1"}, "P1": {"X"}, "Y": {"P2", "Z"}, "P2": {"Y"}, "Z": {"Y"}}
    assert p.score_a_b({"P1"}, "P2", reach, partners) == (1, 2)
    assert p.score_a_b({"X"}, "P2", reach, partners) == (0, 1)      # target widened by its partner
    assert p.score_a_b({"P1"}, "Y", reach, partners) == (0, 1)      # probe widened by its partner
    assert p.score_a_b({"P2"}, "P1", reach, partners) == (0, 0)     # direction counts
    assert p.score_c({"X"}, "P1", partners) == 2
    assert p.score_c({"P2"}, "Z", partners) == 1
    assert p.score_c({"X"}, "Z", partners) == 0


def test_permuted_keeps_degrees():
    import numpy as np
    partners = {"a": {"b", "c"}, "b": {"a"}, "c": {"a"}, "d": set()}
    pp = p.permuted(partners, np.random.default_rng(1))
    assert sorted(len(v) for v in pp.values()) == sorted(len(v) for v in partners.values())
    assert all(a in pp[b] for a, bs in pp.items() for b in bs)       # still symmetric


def test_truth_against_dmso(tmp_path):
    import gzip
    import pytest
    h5py = pytest.importorskip("h5py")
    import numpy as np
    probes, cols = ["pA", "pB", "pC"], [f"s{i}" for i in range(7)]
    # four DMSO wells (MAD of pA = 0.1, pB = 0, pC = 0.1), two wells of D at dose 10, one at dose 1
    m = np.array([[0.0, 1, 0.0], [0.1, 1, 0.1], [-0.1, 1, -0.1], [0.2, 1, np.nan],
                  [0.5, 5, 0.3], [0.7, 5, 0.3], [9.0, 9, 9.0]])
    with h5py.File(tmp_path / "m.gctx", "w") as h:
        h["0/DATA/0/matrix"] = m
        h["0/META/ROW/id"] = np.array([p.encode() for p in probes])
        h["0/META/COL/id"] = np.array([c.encode() for c in cols])
    rows = ["inst_id\tcell_id\tpert_dose\tpert_id\tpert_type"]
    rows += [f"s{i}\tX\t0\tDMSO\tctl_vehicle" for i in range(4)]
    rows += ["s4\tX\t10\tD\ttrt_cp", "s5\tX\t10\tD\ttrt_cp", "s6\tX\t1\tD\ttrt_cp"]
    with gzip.open(tmp_path / "inst.gz", "wt") as f:
        f.write("\n".join(rows) + "\n")
    tv = p.truth(tmp_path / "inst.gz", tmp_path / "m.gctx")
    # pA: |0.6 - 0.05| = 0.55 >= 4.45 * 0.1; pB dropped (MAD 0); pC: 3 DMSO values, |0.3 - 0| < 4.45 * 0.1
    assert tv == {("D", "X"): {"pA": 1, "pC": 0}}
