"""specs/051: the register parses and the guard refuses what it should."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dataset_roles as dr  # noqa: E402


def _reg(tmp_path, body):
    p = tmp_path / "d.tsv"
    p.write_text("id\tslice\trole\tsince\tlicence\tsource\tnote\n" + body)
    return p


def test_committed_register_parses():
    reg = dr.load()
    assert reg["mpbiopath_curator_heldout"]["role"] == "SELECT"
    assert reg["decryptm1"]["role"] == "FIT"


def test_guard(tmp_path):
    p = _reg(tmp_path, "a\tx\tFIT\td\tl\ts\tn\nb\tx\tSELECT\td\tl\ts\tn\nc\tx\tTEST\td\tl\ts\tn\n")
    assert dr.require_role("a", "fit", p) == "FIT"
    assert dr.require_role("b", "select", p) == "SELECT"
    assert dr.require_role("c", "score", p) == "TEST"
    for ds, purpose in (("b", "fit"), ("c", "fit"), ("c", "select")):
        with pytest.raises(SystemExit):
            dr.require_role(ds, purpose, p)
    with pytest.raises(SystemExit):
        dr.require_role("unregistered", "score", p)
    with pytest.raises(ValueError):
        dr.require_role("a", "tune", p)


def test_bad_register(tmp_path):
    with pytest.raises(ValueError):
        dr.load(_reg(tmp_path, "a\tx\tTRAIN\td\tl\ts\tn\n"))
    with pytest.raises(ValueError):
        dr.load(_reg(tmp_path, "a\tx\tFIT\td\tl\ts\tn\na\ty\tFIT\td\tl\ts\tn\n"))
