import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
from korkut_hops import check_pins, hops_from  # noqa: E402


def test_check_pins_refuses_a_different_perturbation():
    check_pins(3, 3, "P", "AK")
    with pytest.raises(SystemExit):
        check_pins(5, 3, "P", "AK")      # e.g. DS_PIN_SCOPE=all exported in the shell


def test_hops_from():
    assert hops_from({"a": ["b"], "b": ["c"]}, {"a"}) == {"a": 0, "b": 1, "c": 2}
