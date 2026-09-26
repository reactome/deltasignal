"""Pins for case_oracle.py's pure functions."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from case_oracle import classify, signed_parities  # noqa: E402


def test_parities_track_inhibition():
    adj = {"G": {("R", 1)}, "R": {("X", 1)}, "X": {("T", -1)}}
    assert signed_parities(adj, {"G"}, {"T"}) == {-1}
    adj["G"].add(("T", 1))
    assert signed_parities(adj, {"G"}, {"T"}) == {1, -1}


def test_classify_strict_lenient_none_and_sign():
    strict = {"G": {("R", 1)}, "R": {("T", 1)}}
    assert classify(strict, strict, {"G"}, {"T"}, "0", "0") == ("strict", "matches")      # KD -> DOWN
    assert classify(strict, strict, {"G"}, {"T"}, "2", "0") == ("strict", "opposite_only")  # KD -> UP needs -1
    assert classify(strict, strict, {"G"}, {"T"}, "1", "0") == ("strict", "n/a")
    lenient = {"G": {("C", 1)}, "C": {("T", 1)}}
    assert classify({}, lenient, {"G"}, {"T"}, "0", "0") == ("lenient", "n/a")
    assert classify({}, {}, {"G"}, {"T"}, "0", "0") == ("none", "n/a")


from case_oracle import USE, MADE, oracle_graph  # noqa: E402

H = lambda *xs: ["R-HSA-" + x for x in xs]


def test_no_hop_between_siblings_of_a_set():
    S, P, T, R1, OUT = H("S", "PMAIP1", "tBID", "R1", "OUT")
    adj = oracle_graph([(R1, "Reaction", S, OUT, "", "", "")], [], [(S, P), (S, T)])
    assert signed_parities(adj, {P}, {R1}) == {1}          # a member acts as the set
    assert signed_parities(adj, {P}, {T}) == set()         # but never reaches a sibling


def test_a_produced_set_hands_on_to_its_members():
    A, S, M1, M2, Z, R1, R2 = H("A", "S", "M1", "M2", "Z", "R1", "R2")
    rows = [(R1, "Reaction", A, S, "", "", ""), (R2, "Reaction", M1, Z, "", "", "")]
    adj = oracle_graph(rows, [], [(S, M1), (S, M2)])
    assert signed_parities(adj, {A}, {Z}) == {1}


def test_assembly_one_level_and_lenient_release():
    A, B, AB, P, R1 = H("A", "B", "AB", "P", "R1")
    rows = [(R1, "Reaction", AB, P, "", "", "")]
    comp = [(AB, A), (AB, B)]
    assert signed_parities(oracle_graph(rows, comp, []), {A}, {P}) == {1}
    assert signed_parities(oracle_graph(rows, comp, []), {AB}, {B}) == set()
    assert signed_parities(oracle_graph(rows, comp, [], lenient=True), {AB}, {B}) == {1}


def test_a_set_component_is_reached_through_its_use_node():
    # PDGF shape: member M -> set S (a component) -> complex C -> reaction R1.
    M, S, C, R1, OUT = H("M", "S", "C", "R1", "OUT")
    adj = oracle_graph([(R1, "Reaction", C, OUT, "", "", "")], [(C, S)], [(S, M)])
    assert signed_parities(adj, {M}, {OUT}) == {1}


def test_a_small_molecule_can_be_a_readout_but_not_a_carrier():
    A, R1, R2, Z = H("A", "R1", "R2", "Z")
    SM = "R-ALL-113592"
    rows = [(R1, "Reaction", A, SM, "", "", ""), (R2, "Reaction", SM, Z, "", "", "")]
    adj = oracle_graph(rows, [], [])
    assert signed_parities(adj, {A}, {SM}) == {1}       # reached as a readout
    assert signed_parities(adj, {A}, {Z}) == set()      # but carries nothing on


def test_both_parities_are_reported_as_such():
    adj = {"G": {("R", 1)}, "R": {("T", 1), ("X", -1)}, "X": {("T", 1)}}
    assert classify(adj, adj, {"G"}, {"T"}, "0", "0") == ("strict", "both_parities")


from route_breaks import first_break  # noqa: E402


def test_first_break_types_the_missing_step():
    klass = {"P": "EntityWithAccessionedSequence", "C": "Complex", "R": "Reaction"}.get
    expand = lambda s: {s}
    assert first_break(["P", "C", "R"], {"P"}, expand, {"P", "C"}, klass)[0] == \
        "component -> complex | node exists, edge missing"
    assert first_break(["P", "C", "R"], {"P"}, expand, {"P"}, klass)[0] == \
        "component -> complex | node absent"
    assert first_break(["P", "C"], {"P", "C"}, expand, set(), klass)[0] == "whole route reached"
    assert first_break(["P", "C"], set(), expand, set(), klass)[0] == \
        "route starts from another form of the gene"


def test_release_of_a_set_component_does_not_reach_a_sibling():
    A, B, S, C, R1, OUT = H("A", "B", "S", "C", "R1", "OUT")
    # A stands in for S, S is a component of C; lenient releases S from C.
    adj = oracle_graph([(R1, "Reaction", C, OUT, "", "", "")], [(C, S)], [(S, A), (S, B)], lenient=True)
    assert signed_parities(adj, {A}, {B}) == set()


def test_shortest_path_ties_do_not_follow_set_order():
    # Two equal-length routes S->A->G and S->B->G. The route must not depend on
    # the iteration order of the adjacency sets (route_breaks.py read 186 vs
    # 166 "edge missing" on identical input before this was sorted).
    from curator_oracle import shortest_path
    for order in (["A", "B"], ["B", "A"]):
        adj = {"S": dict.fromkeys(order).keys(), "A": {"G"}, "B": {"G"}}
        assert shortest_path(adj, {"S"}, {"G"}) == ["S", "A", "G"]
    assert shortest_path({"S": {"A"}}, {"S", "T"}, {"A"}) == ["S", "A"]


def test_faithful_comp_pairs_follows_what_the_generator_decomposes():
    from case_oracle import faithful_comp_pairs
    # K is a root (used, never produced) and contains nested N (unproduced) and
    # Q (produced). P is produced. N2 is unproduced but nested only inside P, and
    # no reaction uses it, so the generator never decomposes it.
    rows = [("R1", "Reaction", "K", "P", "", "", ""),
            ("R2", "Reaction", "X", "Q", "", "", ""),
            ("R3", "Reaction", "P", "Z", "", "", "")]
    comp = {("K", "N"), ("N", "a"), ("K", "Q"), ("Q", "q"), ("K", "b"),
            ("P", "c"), ("P", "N2"), ("N2", "d")}
    got = faithful_comp_pairs(rows, comp)
    assert got == {("K", "N"), ("N", "a"), ("K", "Q"), ("K", "b")}
    # (K, Q) is a join onto the produced Q; Q's own pairs are not representable.
    assert ("Q", "q") not in got and ("N2", "d") not in got and ("P", "c") not in got
