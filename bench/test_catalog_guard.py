"""Runtime test of scripts/catalog.sh's LNG_* leak guard (feedback: a guard
needs a runtime test). The guard once killed every build silently when nothing
leaked (grep exit 1 under set -euo pipefail), so both directions are checked by
running the script itself. Nothing is built: the catalog home does not exist
and the generator checkout is missing, so a clean run stops right after the
guard with "no generator checkout"."""
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "catalog.sh"


def run(extra_env):
    env = {k: v for k, v in os.environ.items() if not k.startswith("LNG_")}
    env.update(DS_CATALOG_HOME=str(REPO / ".no-such-catalog-home"),
               LOGIC_NETWORK_GENERATOR=str(REPO / ".no-such-generator"), **extra_env)
    return subprocess.run(["bash", str(SCRIPT), "build"], env=env, capture_output=True, text=True)


def test_clean_environment_passes_the_guard():
    r = run({})
    assert r.returncode == 1 and "no generator checkout" in r.stderr, r.stderr
    assert "exports" not in r.stderr


def test_lng_python_is_exempt():
    r = run({"LNG_PYTHON": "/usr/bin/python3"})
    assert "no generator checkout" in r.stderr and "exports" not in r.stderr, r.stderr


@pytest.mark.parametrize("name", ["LNG_BOUNDARY_HIERARCHY", "LNG_foo", "LNG_PYTHONX"])
def test_a_leaked_generator_flag_is_refused(name):
    r = run({name: "1"})
    assert r.returncode == 1 and "exports" in r.stderr and name in r.stderr, r.stderr


def test_a_multiline_value_does_not_fake_a_leak():
    r = run({"SOME_VAR": "a\nLNG_FAKE=1"})
    assert "exports" not in r.stderr and "no generator checkout" in r.stderr, r.stderr
