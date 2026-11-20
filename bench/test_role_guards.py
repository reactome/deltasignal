"""specs/051 guards at runtime (review 2026-11-20): each fitting or selection script refuses before it reads
anything when its dataset is registered TEST. The arguments name files that do not exist, so a script that read
data before calling require_role would fail with an IO error instead of the guard's SystemExit."""
import importlib.util
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH))
import dataset_roles as dr  # noqa: E402

MISSING = "/nonexistent/specs051"
SCRIPTS = [
    ("analysis/146/fit_confidence.py", [MISSING]),
    ("analysis/144/engagement_dial.py", [MISSING]),
    ("analysis/162/engagement_162.py", [MISSING] * 5),
    ("analysis/162/arm_162.py", [MISSING] * 6),
    ("analysis/162/score_162.py", [MISSING] * 8),
    ("analysis/163/score_163.py", ["cppa"] + [MISSING] * 6),
    ("analysis/164/score_164.py", [MISSING] * 6),
]


def _load(rel):
    spec = importlib.util.spec_from_file_location(f"guarded_{Path(rel).stem}", BENCH / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("rel,argv", SCRIPTS, ids=[s[0] for s in SCRIPTS])
def test_script_refuses_a_test_dataset(rel, argv, monkeypatch):
    try:
        mod = _load(rel)
    except ImportError as e:            # a script whose heavy imports are absent here cannot be exercised
        pytest.skip(f"{rel}: {e}")
    real = dr.load()
    monkeypatch.setattr(dr, "load", lambda path=None: {k: dict(v, role="TEST") for k, v in real.items()})
    monkeypatch.setattr(sys, "argv", [rel] + argv)
    with pytest.raises(SystemExit, match="registered TEST"):
        mod.main()
