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
