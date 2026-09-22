"""Tests for empirical_axis_feasibility.py.

The first version of this script resolved a readout's node_kind catalog-wide
with first-match-wins. 2,785 of 15,616 stable ids in the 92-pathway catalog
appear as `simple_entity` in one pathway and `dissociation_sink` in another, so
the answer depended on directory iteration order: two runs over the same data
reported 1,691 and 1,812 triples. These tests pin the fix.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from empirical_axis_feasibility import (  # noqa: E402
    BARE_SYMBOL,
    SINGLE_ENTITY_KINDS,
    load_pathway_ids,
    load_readout_kinds,
    triples,
)


def _pathway(root, pid, rows):
    """rows: (uuid, node_kind, stable_id)"""
    d = os.path.join(root, pid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "nodes.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["uuid", "node_kind"])
        for uuid, kind, _ in rows:
            w.writerow([uuid, kind])
    with open(os.path.join(d, "stid_to_uuid_mapping.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["uuid", "stable_id"])
        for uuid, _, sid in rows:
            w.writerow([uuid, sid])


def test_kind_is_resolved_per_pathway_not_catalog_wide(tmp_path):
    root = str(tmp_path)
    # The SAME stable id is an entity in one pathway and a sink in another.
    _pathway(root, "R-HSA-111", [("u1", "simple_entity", "R-HSA-999")])
    _pathway(root, "R-HSA-222", [("u2", "dissociation_sink", "R-HSA-999")])
    kinds = load_readout_kinds(root)
    assert kinds[("R-HSA-111", "999")] == "simple_entity"
    assert kinds[("R-HSA-222", "999")] == "dissociation_sink"
    # and crucially it is keyed by pathway, so neither can shadow the other
    assert len([k for k in kinds if k[1] == "999"]) == 2


def test_within_one_pathway_entity_beats_sink_regardless_of_row_order(tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    _pathway(a, "R-HSA-111", [("u1", "dissociation_sink", "R-HSA-999"),
                              ("u2", "simple_entity", "R-HSA-999")])
    _pathway(b, "R-HSA-111", [("u2", "simple_entity", "R-HSA-999"),
                              ("u1", "dissociation_sink", "R-HSA-999")])
    assert load_readout_kinds(a)[("R-HSA-111", "999")] == "simple_entity"
    assert load_readout_kinds(b)[("R-HSA-111", "999")] == "simple_entity"


def test_a_sink_still_counts_as_a_single_entity():
    # A dissociation sink is a handle for a released subunit -- the same gene
    # product -- so it must not be excluded from the measurable set.
    assert "dissociation_sink" in SINGLE_ENTITY_KINDS
    assert "simple_entity" in SINGLE_ENTITY_KINDS
    assert "simple_complex" not in SINGLE_ENTITY_KINDS
    assert "reaction" not in SINGLE_ENTITY_KINDS


def test_bare_symbol_accepts_real_gene_symbols():
    for s in ("BAX", "CDKN1A", "EGFR", "MT2A", "HLA-A", "TP53", "IFITM3"):
        assert BARE_SYMBOL.match(s), s


def test_bare_symbol_rejects_everything_a_transcript_assay_cannot_measure():
    for s in (
        "p-S133-CREB1",              # phospho-form: a state, not an abundance
        "p-S317_S345-CHEK1",
        "K6PolyUb-p-BRCA1",          # ubiquitin-form
        "4xPalmC-CD36",              # lipid modification
        "ACIN1(1-1093)",             # proteolytic fragment
        "APC(778-2843)",
        "5649637",                   # a bare numeric id, i.e. unnamed
        "",                          # no name at all
        "EXO1,DNA2:BLM,WRN",         # a set or complex expression
        "beta-catenin",              # lowercase: not an HGNC symbol
    ):
        assert not BARE_SYMBOL.match(s), s


def test_pathway_list_maps_names_to_ids_and_tolerates_a_bare_number(tmp_path):
    p = tmp_path / "list.tsv"
    p.write_text("id\tpathway_name\nR-HSA-73894\tDNA_Repair\n453279\tMitotic_G1\n")
    m = load_pathway_ids(str(p))
    assert m["DNA_Repair"] == "R-HSA-73894"
    assert m["Mitotic_G1"] == "R-HSA-453279"


def test_triples_deduplicates_on_pathway_gene_readout():
    rows = [
        {"pathway": "P", "gene": "G", "key_output": "1"},
        {"pathway": "P", "gene": "G", "key_output": "1"},   # duplicate
        {"pathway": "P", "gene": "G", "key_output": "2"},
        {"pathway": "Q", "gene": "G", "key_output": "1"},
    ]
    assert len(triples(rows)) == 3
