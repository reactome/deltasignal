"""experimental_gap.py: the join must be exact and the quadrants right."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from experimental_gap import report, split  # noqa: E402


def case(p, g, d, k, pred, exp, valid="1", cat="pass"):
    return {"pathway": p, "gene": g, "direction": d, "key_output": k, "predicted": pred,
            "expected": exp, "valid": valid, "category": cat}


MPB = {("P", "A", "0", "1"): ("0", "0"), ("P", "A", "0", "2"): ("2", "1"),
       ("P", "B", "2", "1"): ("2", "2"), ("Q", "C", "0", "9"): ("1", "2")}


def test_quadrants_and_invalid_counts_as_wrong():
    rows = split(MPB, [case("P", "A", "0", "1", "2", "0"),            # only MPB
                       case("P", "A", "0", "2", "2", "2"),            # only DS
                       case("P", "B", "2", "1", "1", "2", valid="0"),  # invalid -> only MPB
                       case("Q", "C", "0", "9", "0", "1")])           # neither
    assert [r["quadrant"] for r in rows] == ["MPB/-", "-/DS", "MPB/-", "-/-"]
    text = report(rows)
    assert "MP-BioPath 2, DeltaSignal 1, gap +1" in text
    assert "either 3" in text


def test_a_label_disagreement_is_an_error():
    with pytest.raises(SystemExit):
        split(MPB, [case("P", "A", "0", "1", "0", "2")])


def test_a_case_missing_on_either_side_is_an_error():
    with pytest.raises(SystemExit):
        split(MPB, [case("P", "A", "0", "1", "0", "0")])               # 3 MPB cases unmatched
    with pytest.raises(SystemExit):
        split({}, [case("P", "A", "0", "1", "0", "0")])
