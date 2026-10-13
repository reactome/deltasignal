"""specs/051 R4: the drug -> target map follows the registered rule and amendment A1."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import drug_targets as dt  # noqa: E402


def _a(mol, tid, rel, v, typ="Kd", units="nM"):
    return {"molecule_chembl_id": mol, "target_chembl_id": tid, "standard_relation": rel,
            "standard_value": str(v), "standard_type": typ, "standard_units": units}


T = {"t1": ("P1", "EGFR"), "t2": ("P2", "MET"), "t3": ("P3", "KIT"), "tc": None, "tx": ("P9", None)}


def test_cutoff_relation_and_fallbacks():
    acts = [_a("m1", "t1", "=", 2), _a("m1", "t2", "=", 1500), _a("m1", "t1", "=", 9),
            _a("m1", "tc", "=", 1), _a("m1", "tx", "=", 1),
            _a("m2", "t3", ">", 30000),
            _a("m3", "t1", "=", 5000)]
    mols = {"DrugA": "m1", "DrugB": "m2", "DrugC": "m3", "Antibody": ""}
    hand = {"DrugA": ["EGFR", "ERBB2", "KIT"], "DrugB": ["KIT"], "DrugC": ["EGFR"], "Antibody": ["ERBB2"]}
    rows, tally = dt.build(acts, T, mols, hand, cutoff_nm=1000)
    got = {(d, g): (kd, src) for d, g, _, kd, src in rows}
    # the lowest Kd per gene; MET above the cutoff is out; complex and non-Reactome targets skipped
    assert got[("DrugA", "EGFR")] == (2.0, "kinobeads") and ("DrugA", "MET") not in got
    assert tally["not_single_protein"] == 1 and tally["not_in_reactome"] == 1
    # A1: ERBB2 was never measured for any drug -> kept; KIT was measured (for DrugB) -> decided by measurement
    assert got[("DrugA", "ERBB2")] == ("", "hand_unmeasured") and ("DrugA", "KIT") not in got
    # no target at the cutoff (inactive, or too weak) -> the hand table
    assert got[("DrugB", "KIT")] == ("", "hand") and got[("DrugC", "EGFR")] == ("", "hand")
    assert got[("Antibody", "ERBB2")] == ("", "hand")
    assert tally["drugs_kinobeads"] == 1 and tally["drugs_hand"] == 3
    rows10k, _ = dt.build(acts, T, mols, hand, cutoff_nm=10000)
    assert ("DrugA", "MET", "P2", 1500.0, "kinobeads") in rows10k
    assert ("DrugC", "EGFR", "P1", 5000.0, "kinobeads") in rows10k


def test_built_map_is_well_formed():
    """The map lives outside the repo (CC BY-SA source); checked where it exists."""
    import pytest
    path = Path.home() / "deltasignal-catalogs/external/klaeger2017_chembl37/drug_targets_1uM.tsv"
    if not path.exists():
        pytest.skip("drug->target map not built on this machine")
    rows = dt.read_tsv(path)
    assert {r["source"] for r in rows} <= {"kinobeads", "hand", "hand_unmeasured"}
    assert len({(r["drug"], r["gene"]) for r in rows}) == len(rows)
    assert {r["drug"] for r in rows} == set(dt.load_molecules())
    for r in rows:
        if r["source"] == "kinobeads":
            assert 0 < float(r["kd_app_nM"]) <= dt.PRIMARY_CUTOFF_NM
