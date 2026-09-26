"""Pins for release_skew_reason (specs/025): which unscorable cases cannot
measure the generator or the solver, and so are left out of the stats."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import GENE_NAME_CORRECTIONS, release_skew_reason  # noqa: E402

PID = "1"


def cache(genes=(), participants=(), dbids=()):
    c = {("genes", PID): set(genes), ("participants", PID): set(participants)}
    for d, exists in dbids:
        c[("dbid", d)] = exists
    return c


def test_unknown_gene_name_is_unresolved():
    assert release_skew_reason(PID, "P", "CACNAD1", [], {"CACNAD1": []}, "9", ["o"], {}, cache()) \
        == "gene_name_unresolved"


def test_gene_absent_from_this_release_pathway():
    got = release_skew_reason(PID, "P", "COL1A1", [], {"COL1A1": ["R-HSA-1"]}, "9", ["o"], {},
                              cache(genes={"OTHER"}))
    assert got == "gene_not_in_pathway"


def test_gene_in_pathway_but_unpinnable_is_KEPT():
    # MIR675: in Reactome's pathway and our network, no root form -- ours to fix.
    got = release_skew_reason(PID, "P", "MIR675", [], {"MIR675": ["R-HSA-2"]}, "9", ["o"], {},
                              cache(genes={"MIR675"}))
    assert got == ""


def test_readout_gone_from_the_release():
    got = release_skew_reason(PID, "P", "G", ["u"], {"G": ["s"]}, "69589", [], {},
                              cache(dbids=[("69589", False)]))
    assert got == "readout_not_in_release"


def test_readout_exists_but_left_the_pathway():
    got = release_skew_reason(PID, "P", "G", ["u"], {"G": ["s"]}, "7", [], {"7": "R-HSA-7"},
                              cache(participants={"R-HSA-8"}))
    assert got == "readout_not_in_pathway"


def test_readout_in_pathway_but_unresolved_is_KEPT():
    got = release_skew_reason(PID, "P", "G", ["u"], {"G": ["s"]}, "7", [], {"7": "R-HSA-7"},
                              cache(participants={"R-HSA-7"}))
    assert got == ""


def test_corrections_only_name_genes_that_exist():
    assert GENE_NAME_CORRECTIONS[("DNA_Double-Strand_Break_Repair", "PRKDC1")] == "PRKDC"
    assert all(v.isupper() for v in GENE_NAME_CORRECTIONS.values())
