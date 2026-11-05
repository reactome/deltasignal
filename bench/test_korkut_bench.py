"""specs/062: the korkut2015 bench's pure parts (no Neo4j, no API, no values)."""
import math
from pathlib import Path

import pytest

import korkut_bench as kb


def test_call_thresholds():
    t = math.log2(1.5)
    assert kb.call(t, t) == "up" and kb.call(-t, t) == "down"
    assert kb.call(t - 1e-9, t) == "not" and kb.call(0.0, t) == "not"


def test_drug_set_call_agree_conflict_none():
    t = math.log2(1.5)
    assert kb.drug_set_call([1.0, 0.7], t) == "up"
    assert kb.drug_set_call([1.0, 0.1], t) == "up"          # one dose regulated, one not
    assert kb.drug_set_call([1.0, -1.0], t) == "conflict"
    assert kb.drug_set_call([0.1, -0.2], t) == "not"


def test_conditions_drop_doses_and_sort(tmp_path):
    p = tmp_path / "expr_index.txt"
    p.write_text("901|1.5        1   0\nPLX|60,AK|5 1 1\nAK|10 1 0\n")
    assert kb.conditions(p) == [(0, ("901",)), (1, ("AK", "PLX")), (2, ("AK",))]


def test_readout_stids_total_and_phospho():
    readouts = [{"readout": "AKT", "kind": "total", "sites": "AKT1,AKT2"},
                {"readout": "AKTpS473", "kind": "phospho", "sites": "AKT1:473,AKT2:474"}]
    g2s = {"AKT1": ["R-1", "R-1p"], "AKT2": ["R-2"]}
    s2s = {("AKT1", 473): ["R-1p"], ("AKT2", 474): ["R-2p"]}
    out = kb.readout_stids(readouts, g2s, s2s)
    assert out == {"AKT": ["R-1", "R-1p", "R-2"], "AKTpS473": ["R-1p", "R-2p"]}


def test_maps_are_well_formed():
    drugs = kb.read_tsv(kb.MAPS / "drugs.tsv")
    assert len(drugs) == 12 and len({d["code"] for d in drugs}) == 12
    ro = kb.read_tsv(kb.MAPS / "readouts.tsv")
    assert len(ro) == 82 and {r["kind"] for r in ro} == {"phospho", "total", "exclude"}
    for r in ro:
        if r["kind"] == "phospho":
            assert all(":" in s and s.split(":")[1].isdigit() for s in r["sites"].split(","))
        if r["kind"] == "exclude":
            assert r["note"]


@pytest.mark.skipif(not (kb.DATA / "expr_index.txt").exists(), reason="korkut2015 design not on this machine")
def test_every_design_code_is_mapped():
    codes = {c for _, ds in kb.conditions() for c in ds}
    assert codes == {d["code"] for d in kb.read_tsv(kb.MAPS / "drugs.tsv")}


def test_measured_calls_refuses_misaligned_columns(tmp_path, monkeypatch):
    (tmp_path / "node_Index.csv").write_text("p53\np21\n")
    (tmp_path / "expr.csv").write_text("0.1,0.9\n7,0.1,0.9\n")
    monkeypatch.setattr(kb, "DATA", tmp_path)
    design = [(0, ("NT",)), (1, ("AK",))]
    with pytest.raises(SystemExit, match="misalign"):
        kb.measured_calls(kb.THRESHOLD, design)
    (tmp_path / "expr.csv").write_text("0.1,0.9\n-0.7,0.1\n")
    m = kb.measured_calls(kb.THRESHOLD, design)
    assert m[(("NT",), "p21")] == "up" and m[(("AK",), "p53")] == "down"


def test_readout_nodes_reads_only_own_nodes_in_a_context_build():
    r_stids = {"AKT_pS473": ["R-HSA-1"], "MAPK_pT202": ["R-HSA-2"]}
    s2u = {"R-HSA-1": ["p1", "q1", "p1"], "R-HSA-2": ["q2"]}
    assert kb.readout_nodes(r_stids, s2u) == {"AKT_pS473": ["p1", "q1"], "MAPK_pT202": ["q2"]}
    # a readout only a neighbour pathway holds is not P's readout
    assert kb.readout_nodes(r_stids, s2u, own={"p1"}) == {"AKT_pS473": ["p1"]}


def _edge(s, t, kind):
    return {"source_id": s, "target_id": t, "edge_type": kind}


def test_abundance_copies_keep_starting_and_translated_copies_only():
    # A (start) binds B in r1 -> A:B; A is phosphorylated in r2 -> pA; mRNA m is translated in r3 -> A2 (a new copy);
    # pA is degraded by D (depletion: not a source); a set variant V feeds pool Q (identity edge: Q is made from V)
    rows = [_edge("A", "r1", "input"), _edge("B", "r1", "input"), _edge("r1", "AB", "output"),
            _edge("A", "r2", "input"), _edge("K", "r2", "catalyst"), _edge("r2", "pA", "output"),
            _edge("m", "r3", "input"), _edge("r3", "A2", "output"),
            _edge("D", "pA", "depletion"), _edge("V", "Q", "variant_pool")]
    nodes = ["A", "AB", "pA", "A2", "V", "Q"]
    assert kb.abundance_copies(nodes, rows) == ["A", "A2", "V"]


def test_abundance_copies_catalyst_in_set_is_not_a_source():
    # the protein catalysing a reaction that makes another of its copies does not make that copy from itself
    rows = [_edge("X", "r", "input"), _edge("A", "r", "catalyst"), _edge("r", "A3", "output")]
    assert kb.abundance_copies(["A", "A3"], rows) == ["A", "A3"]


def test_abundance_copies_assembled_and_released_copies_are_made():
    # specs/134 amendment 1: a starting complex A:P built from A (assembly), A released from the bound complex A:Q
    # (dissociation) and a set pool S of A (set_member) are made from A; a complex built only from partners is not
    rows = [_edge("A", "AP", "assembly"), _edge("P", "AP", "assembly"), _edge("AQ", "A4", "dissociation"),
            _edge("A", "S", "set_member"), _edge("P", "PR", "assembly")]
    assert kb.abundance_copies(["A", "AP", "AQ", "A4", "S"], rows) == ["A", "AQ"]


def test_abundance_copies_pool_of_the_protein_carries_it():
    # specs/134 amendment 2: complex C is made from pool Q, which pools A's copies (a set input), so C is a bound form;
    # D is made from pool R of partners only, and a partner B released from A:B does not carry A
    rows = [_edge("A", "Q", "variant_pool"), _edge("Q", "r1", "input"), _edge("r1", "C", "output"),
            _edge("B", "R", "variant_pool"), _edge("R", "r2", "input"), _edge("r2", "D", "output"),
            _edge("AB", "B", "dissociation"), _edge("B", "r3", "input"), _edge("r3", "E", "output")]
    assert kb.abundance_copies(["A", "C", "D", "AB", "E"], rows) == ["A", "D", "AB", "E"]


def test_abundance_copies_loop_without_start_has_none():
    rows = [_edge("pA", "r1", "input"), _edge("r1", "A", "output"), _edge("A", "r2", "input"), _edge("r2", "pA", "output")]
    assert kb.abundance_copies(["A", "pA"], rows) == []


def test_total_readout_nodes_leaves_forms_and_phospho_alone(tmp_path, monkeypatch):
    ro = {"AKT": ["A", "AB"], "AKTpS473": ["pA"]}
    monkeypatch.setattr(kb, "TOTAL_READOUT", "forms")
    assert kb.total_readout_nodes(ro, {"AKT"}, tmp_path) is ro
    monkeypatch.setattr(kb, "TOTAL_READOUT", "abundance")
    (tmp_path / "logic_network.csv").write_text(
        "source_id,target_id,pos_neg,and_or,edge_type,stoichiometry,edge_reaction_id\n"
        "A,r1,pos,and,input,1,R1\nB,r1,pos,and,input,1,R1\nr1,AB,pos,and,output,1,R1\n")
    out = kb.total_readout_nodes(ro, {"AKT"}, tmp_path)
    assert out == {"AKT": ["A"], "AKTpS473": ["pA"]}
