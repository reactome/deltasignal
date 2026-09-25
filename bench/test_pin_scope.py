"""Pins for entry_occurrences (specs/023): pin where the gene enters the
network, compute everything downstream."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from benchmark_vs_mpbiopath import entry_occurrences  # noqa: E402


def test_root_entity_is_pinned_and_its_complexes_are_not():
    # A (the gene) + B -> AB complex -> reaction R -> AB:C
    adj = {"A": ["AB"], "B": ["AB"], "AB": ["R"], "R": ["ABC"]}
    assert entry_occurrences(["A", "AB", "ABC"], adj) == ["A"]


def test_independent_occurrences_are_each_an_entry():
    # Two branches that never meet: both are where the gene enters.
    adj = {"A1": ["X"], "A2": ["Y"]}
    assert entry_occurrences(["A1", "A2"], adj) == ["A1", "A2"]


def test_first_mid_pathway_occurrence_is_pinned_when_there_is_no_root():
    # The gene is only produced by a reaction (transcription) and flows on.
    adj = {"DNA": ["TX"], "TX": ["P"], "P": ["PC"]}
    assert entry_occurrences(["P", "PC"], adj) == ["P"]


def test_members_that_only_reach_each_other_fall_back_to_all():
    adj = {"A": ["B"], "B": ["A"]}
    assert sorted(entry_occurrences(["A", "B"], adj)) == ["A", "B"]


def test_a_cycle_downstream_of_the_root_is_not_pinned():
    adj = {"A": ["C1"], "C1": ["C2"], "C2": ["C1"]}
    assert entry_occurrences(["A", "C1", "C2"], adj) == ["A"]


from benchmark_vs_mpbiopath import root_occurrences  # noqa: E402


def test_root_scope_keeps_only_nodes_with_no_incoming_edge():
    indeg = {"AB": 2, "R": 1, "ABC": 1}           # A and B are roots
    assert root_occurrences(["A", "AB", "ABC"], indeg) == ["A"]


def test_a_root_complex_containing_the_gene_is_pinned():
    # A root complex that is not decomposed still counts: it contains the gene.
    assert root_occurrences(["AB", "ABC"], {"ABC": 1}) == ["AB"]


def test_root_scope_does_not_invent_a_pin_for_a_gene_with_no_root_form():
    # Unlike `entry`, a gene only produced mid-pathway is NOT perturbed.
    assert root_occurrences(["P", "PC"], {"P": 1, "PC": 1}) == []


def test_root_scope_drops_duplicates():
    assert root_occurrences(["A", "A", "B"], {}) == ["A", "B"]
