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
