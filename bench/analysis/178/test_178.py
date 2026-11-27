"""specs/178: the release comparison on synthetic pair tables."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_178 as C  # noqa: E402


def test_compare():
    old = {("G", "a"): ("up", "up"), ("G", "b"): ("not", "not"), ("H", "a"): ("not", "up")}
    new = {("G", "a"): ("not", "up"), ("G", "b"): ("not", "up"), ("K", "a"): ("up", "up")}
    r = C.compare(old, new)
    assert r["pairs_both"] == 2 and r["only_old"] == 1 and r["only_new"] == 1
    assert r["model_calls_moved"] == 1 and r["model_transitions"] == {"up->not": 1}
    assert r["reach_calls_moved"] == 1 and r["reach_transitions"] == {"not->up": 1}
    assert r["only_new_by_readout"] == {"a": 1}
