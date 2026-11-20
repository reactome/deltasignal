"""Runtime test of `scripts/catalog.sh bench`'s served-build refusal (review 2026-11-20). On the dev box the API
serves PATHWAY_CATALOG from the untracked .env while `current` is the frozen build, so a plain `bench` always
refused with only "recreate it". The refusal must name the cause and the --build fix, and `--build` must select a
build by directory while comparing the manifest's build id (a copy such as *-plus-leads keeps the id).

Nothing is benchmarked: every case stops at the served-build check, before any benchmark starts. A fake API
answers /api/health. The catalog home lives under the repository, not /tmp, because catalog.sh refuses a tmpfs."""
import http.server
import json
import os
import shutil
import subprocess
import sys
import threading
import uuid
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "catalog.sh"
SERVED = "20261113-1311_baf608b_r98fix"
FROZEN = "20261026-1349_8b6527e"


@pytest.fixture
def home():
    h = REPO / f".catalog-test-{uuid.uuid4().hex[:8]}"
    for d, bid in ((FROZEN, FROZEN), ("20261113-r98fix-plus-leads", SERVED), ("other-dir", "x")):
        (h / "builds" / d).mkdir(parents=True)
        (h / "builds" / d / "BUILD.json").write_text(json.dumps({"build_id": bid, "status": "complete"}))
    (h / "current").symlink_to(f"builds/{FROZEN}")
    # a docker that knows no container, so a run that passes the served-build check stops at the container check
    # and never starts a benchmark, on a dev box with a live julia-api as in CI
    (h / "bin").mkdir()
    (h / "bin" / "docker").write_text("#!/bin/sh\nexit 1\n")
    (h / "bin" / "docker").chmod(0o755)
    yield h
    shutil.rmtree(h)


@pytest.fixture
def api():
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            body = json.dumps({"catalog": {"build_id": SERVED}}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass

    srv = http.server.HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_port}"
    srv.shutdown()


def run(home, api, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("DS_", "LNG_"))}
    env.update(DS_CATALOG_HOME=str(home), DELTASIGNAL_BASE=api, LNG_PYTHON=sys.executable,
               PATH=f"{home / 'bin'}{os.pathsep}{os.environ['PATH']}")
    r = subprocess.run(["bash", str(SCRIPT), "bench", *args], env=env, capture_output=True, text=True)
    if "uncommitted change" in r.stderr:
        pytest.skip("src/ or bench/ is uncommitted here; the bench refuses before the served-build check")
    return r


def test_plain_bench_names_the_cause_and_the_fix(home, api):
    r = run(home, api)
    assert r.returncode != 0
    assert f"serves '{SERVED}'" in r.stderr and f"scoring '{FROZEN}'" in r.stderr
    assert "PATHWAY_CATALOG" in r.stderr and "--build" in r.stderr
    assert "20261113-r98fix-plus-leads" in r.stderr      # the matching directory is listed
    assert "other-dir" not in r.stderr


def test_build_option_compares_the_manifest_id(home, api):
    # the served build's directory passes the served-build check and stops at the (fake) container check
    r = run(home, api, "--build", "20261113-r98fix-plus-leads")
    assert r.returncode != 0 and "serves '" not in r.stderr
    assert "julia-api (dev compose) is not running" in r.stderr
    r = run(home, api, "--build", "other-dir")
    assert r.returncode != 0 and "scoring 'x' (builds/other-dir)" in r.stderr


@pytest.mark.parametrize("bad", ["../builds", "a/b", ".."])
def test_build_option_refuses_a_path(home, api, bad):
    r = run(home, api, "--build", bad)
    assert r.returncode != 0 and "not a path" in r.stderr


def test_build_option_refuses_an_unknown_build(home, api):
    r = run(home, api, "--build", "nope")
    assert r.returncode != 0 and "no build nope" in r.stderr
