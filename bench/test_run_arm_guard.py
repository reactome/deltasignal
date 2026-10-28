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


def _run(*args, **kw):
    return subprocess.run(["bash", str(SCRIPT), "probe", *args], capture_output=True, text=True, timeout=30, **kw)


def test_limit_is_refused_on_the_drug_axes(tmp_path):
    f = tmp_path / "x.tsv"
    f.write_text("")
    for extra in (["--limit", "2"], ["--max-edges", "100"]):
        r = _run("--axis", "decryptm", "--decryptm-cases", str(f), "--decryptm-phospho", str(f), *extra)
        assert r.returncode == 1
        assert "apply only to the curator and experimental axes" in r.stderr


def test_decryptm_inputs_without_the_axis_are_refused(tmp_path):
    f = tmp_path / "x.tsv"
    f.write_text("")
    r = _run("--axis", "korkut", "--decryptm-cases", str(f))
    assert r.returncode == 1
    assert "without --axis decryptm" in r.stderr


def test_drug_axis_allowlist_names_the_shared_knobs():
    text = SCRIPT.read_text()
    line = next(l for l in text.splitlines() if l.strip().startswith("DRUG_KNOBS="))
    for knob in ("DS_DOWN_CUTOFF", "DS_UP_CUTOFF", "DS_PIN_NETWORK_ROOT", "DS_SOLVE_TIMEOUT", "DS_PERTURB_UI_DOWN"):
        assert knob in line
