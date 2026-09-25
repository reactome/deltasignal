"""Pins for arm_compare.py: it must pair only cases valid in BOTH arms, count
fixed/broke correctly, and compute macro-F1 over DOWN/NORMAL/UP."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arm_compare import compare, macro_f1  # noqa: E402


def row(pred, exp, valid="1"):
    return {"predicted": pred, "expected": exp, "valid": valid}


def test_only_cases_valid_in_both_arms_are_paired():
    a = {("P", "G", "0", "1"): row("0", "0"), ("P", "G", "0", "2"): row("1", "0")}
    b = {("P", "G", "0", "1"): row("1", "0"), ("P", "G", "0", "2"): row("0", "0", valid="0")}
    r = compare(a, b, lambda k: True)
    assert r["shared"] == 1 and r["valid_base"] == 2 and r["valid_arm"] == 1
    assert (r["fixed"], r["broke"]) == (0, 1)       # the invalid case is not scored as a gain


def test_fixed_broke_and_pathway_direction():
    a = {("P", "G", "2", str(i)): row("1", "2") for i in range(3)}
    a[("Q", "H", "0", "9")] = row("0", "0")
    b = {k: row("2", "2") for k in a if k[0] == "P"}
    b[("Q", "H", "0", "9")] = row("1", "0")
    r = compare(a, b, lambda k: True)
    assert (r["fixed"], r["broke"]) == (3, 1)
    assert r["per_pathway"] == {"P": 3, "Q": -1}
    assert (r["up"], r["down"], r["perturbations"]) == (1, 1, 2)


def test_macro_f1_is_the_mean_over_three_classes():
    assert macro_f1([("0", "0"), ("1", "1"), ("2", "2")]) == 1.0
    # all predicted NORMAL: only the NORMAL class scores
    got = macro_f1([("1", "0"), ("1", "1"), ("1", "2")])
    assert abs(got - (2 * 1 / (2 * 1 + 2 + 0)) / 3) < 1e-12
