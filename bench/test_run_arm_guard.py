"""run_arm.sh refuses an arm NAME that starts with '-'.

`scripts/run_arm.sh --help` used to start a real arm named "--help" against
the current catalog (2026-10-09): the name check allowed a leading dash. The
guard must fire before anything is started, so this runs the script itself.
"""
import pathlib
import subprocess

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "run_arm.sh"


@pytest.mark.parametrize("name", ["--help", "-h", "-x"])
def test_a_leading_dash_is_refused_before_anything_runs(name):
    r = subprocess.run(["bash", str(SCRIPT), name], capture_output=True, text=True, timeout=30)
    assert r.returncode == 1
    assert "not start with '-'" in r.stderr
    assert "usage:" in r.stderr


def test_decryptm_axis_needs_its_inputs_before_anything_runs():
    r = subprocess.run(["bash", str(SCRIPT), "probe", "--axis", "decryptm"],
                       capture_output=True, text=True, timeout=30)
    assert r.returncode == 1
    assert "--axis decryptm needs --decryptm-cases and --decryptm-phospho" in r.stderr


def test_decryptm_inputs_must_exist():
    r = subprocess.run(["bash", str(SCRIPT), "probe", "--axis", "decryptm", "--decryptm-cases", "/nonexistent"],
                       capture_output=True, text=True, timeout=30)
    assert r.returncode == 1
    assert "--decryptm-cases needs an existing file" in r.stderr
