#!/usr/bin/env python3
"""Validate inhibition formulas against synthetic test pathways with known expected
qualitative behavior. Runs each DS_INHIBITION_MODE against each case and reports
which formulas predict the right direction.

Test pathways live in the logic-network-generator's output/ directory:
  TEST_simple_inhibit_R-HSA-90000001    A → R ⊣ C → B
  TEST_chain_inhibit_R-HSA-90000002     A → R1 ⊣ C → X → R2 → B
  TEST_negfeedback_loop_R-HSA-90000003  A → R1 → B, B ⊣ R1

Created during the 2026-05-28 inhibition-formula investigation; rebuilt by
the corresponding `mkdir + cat` block in the deltasignal commit history if
those directories disappear.

The API container must be restarted with the desired DS_INHIBITION_MODE before
each pass. This driver does that via docker compose then calls /api/parse and
/api/solve for each scenario."""

import json
import os
import subprocess
import sys
import time
from urllib.request import Request, urlopen

DS = "http://127.0.0.1:8080"

# UUIDs (kept short readable)
A1, C1, R1, B1 = (
    "00000000-0000-0000-0000-000000000001",
    "00000000-0000-0000-0000-000000000002",
    "00000000-0000-0000-0000-000000000003",
    "00000000-0000-0000-0000-000000000004",
)
A2, C2, R1a, X2, R2a, B2 = (
    "00000000-0000-0000-0000-000000000011",
    "00000000-0000-0000-0000-000000000012",
    "00000000-0000-0000-0000-000000000013",
    "00000000-0000-0000-0000-000000000014",
    "00000000-0000-0000-0000-000000000015",
    "00000000-0000-0000-0000-000000000016",
)
A3, R3, B3 = (
    "00000000-0000-0000-0000-000000000021",
    "00000000-0000-0000-0000-000000000023",
    "00000000-0000-0000-0000-000000000024",
)

# Each scenario: (pathway_dir, label, observations dict, readout_uuid, expected_direction)
SCENARIOS = [
    # ---- simple inhibition: A → R ⊣ C → B ----
    ("TEST_simple_inhibit_R-HSA-90000001",   "simple/baseline",     {},                       B1, "NORMAL"),
    ("TEST_simple_inhibit_R-HSA-90000001",   "simple/C-knockout",   {C1: [0.0, 1.0]},         B1, "UP"),
    ("TEST_simple_inhibit_R-HSA-90000001",   "simple/C-upregulate", {C1: [80.0, 1.0]},        B1, "DOWN"),
    ("TEST_simple_inhibit_R-HSA-90000001",   "simple/A-knockout",   {A1: [0.0, 1.0]},         B1, "DOWN"),
    ("TEST_simple_inhibit_R-HSA-90000001",   "simple/A-upregulate", {A1: [80.0, 1.0]},        B1, "UP"),
    # ---- chain: A → R1 ⊣ C → X → R2 → B ----
    ("TEST_chain_inhibit_R-HSA-90000002",    "chain/baseline",      {},                       B2, "NORMAL"),
    ("TEST_chain_inhibit_R-HSA-90000002",    "chain/C-knockout",    {C2: [0.0, 1.0]},         B2, "UP"),
    ("TEST_chain_inhibit_R-HSA-90000002",    "chain/C-upregulate",  {C2: [80.0, 1.0]},        B2, "DOWN"),
    # ---- neg-feedback loop: A → R1 → B, B ⊣ R1 ----
    ("TEST_negfeedback_loop_R-HSA-90000003", "loop/baseline",       {},                       B3, "NORMAL"),
    ("TEST_negfeedback_loop_R-HSA-90000003", "loop/A-knockout",     {A3: [0.0, 1.0]},         B3, "DOWN"),
    ("TEST_negfeedback_loop_R-HSA-90000003", "loop/A-upregulate",   {A3: [80.0, 1.0]},        B3, "NORMAL"),
    # ^ at UP A, the loop self-regulates: B rises, B inhibits R1, R1 dampens. Net: B settles
    # back toward baseline-ish (true negative feedback). NORMAL is the "correct" loop behavior.
]

# Discretization thresholds (UI scale; internal × 100).
DOWN, NORMAL, UP = 0.5, 1.5, 2.0  # UI < DOWN → DOWN; UI ≥ UP → UP; in between → NORMAL


def api(path, body):
    r = Request(f"{DS}{path}", data=json.dumps(body).encode(),
                headers={"Content-Type": "application/json"}, method="POST")
    with urlopen(r, timeout=120) as x:
        return json.loads(x.read())


def classify(ui_value):
    if ui_value < DOWN:
        return "DOWN"
    if ui_value >= UP:
        return "UP"
    return "NORMAL"


def restart_api_with_mode(mode, **kwargs):
    env = {"DS_AND_MODE": "signed", "DS_INHIBITION_MODE": mode, "DS_OR_MODE": "max", **kwargs}
    env_str = " ".join(f"{k}={v}" for k, v in env.items())
    subprocess.run(
        f"cd /home/awright/gitroot/deltasignal && {env_str} "
        "docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api",
        shell=True, capture_output=True, timeout=60,
    )
    # Wait for health
    for _ in range(40):
        try:
            with urlopen(f"{DS}/api/health", timeout=2):
                time.sleep(2)
                return
        except Exception:
            time.sleep(2)
    raise RuntimeError("API never came up")


def run_pass(mode, mode_label, **env_extra):
    print(f"\n=== {mode_label} ===")
    restart_api_with_mode(mode, **env_extra)
    parsed_cache = {}
    rows = []
    for pathway, label, obs, readout, expected in SCENARIOS:
        if pathway not in parsed_cache:
            parsed_cache[pathway] = api("/api/parse", {"pathway_id": pathway})
        parsed = parsed_cache[pathway]
        net = {"nodes": parsed["nodes"], "edges": parsed["edges"], "pathways": parsed["pathways"]}
        res = api("/api/solve", {"network": net, "observations": obs})
        v_internal = res["node_activities"].get(readout, 0.01)
        v_ui = v_internal * 100.0
        pred = classify(v_ui)
        mark = "✓" if pred == expected else "✗"
        rows.append((label, expected, pred, v_ui, mark))
        print(f"  {mark} {label:30s} expected={expected:6s} got={pred:6s}  B_UI={v_ui:.3f}")
    n_ok = sum(1 for r in rows if r[-1] == "✓")
    print(f"  --- {n_ok}/{len(rows)} correct ---")
    return rows


if __name__ == "__main__":
    print("Synthetic test pathways for inhibition formula validation\n")
    all_results = {}
    for mode_label, env in [
        ("spec, β=0 (current baseline)",  {"DS_INHIBITION_MODE": "spec"}),
        ("spec, β=1",                     {"DS_INHIBITION_MODE": "spec",     "DS_INHIBITOR_BETA": "1"}),
        ("krep, K=0.1",                   {"DS_INHIBITION_MODE": "krep",     "DS_INHIBITOR_K":    "0.1"}),
        ("krep, K=0.05",                  {"DS_INHIBITION_MODE": "krep",     "DS_INHIBITOR_K":    "0.05"}),
        ("divide, ε=0.001",               {"DS_INHIBITION_MODE": "divide",   "DS_INHIBITOR_EPS":  "0.001"}),
        ("divide, ε=0.005",               {"DS_INHIBITION_MODE": "divide",   "DS_INHIBITOR_EPS":  "0.005"}),
        ("devspec, β=2",                  {"DS_INHIBITION_MODE": "devspec",  "DS_INHIBITOR_BETA": "2"}),
    ]:
        mode = env.pop("DS_INHIBITION_MODE")
        all_results[mode_label] = run_pass(mode, mode_label, **env)

    # Summary
    print("\n\n=== SUMMARY ===")
    print(f"{'mode':36s} {'simple':>8s} {'chain':>6s} {'loop':>6s} {'total':>6s}")
    for mode_label, rows in all_results.items():
        n_simple = sum(1 for r in rows if r[0].startswith("simple") and r[-1] == "✓")
        d_simple = sum(1 for r in rows if r[0].startswith("simple"))
        n_chain  = sum(1 for r in rows if r[0].startswith("chain")  and r[-1] == "✓")
        d_chain  = sum(1 for r in rows if r[0].startswith("chain"))
        n_loop   = sum(1 for r in rows if r[0].startswith("loop")   and r[-1] == "✓")
        d_loop   = sum(1 for r in rows if r[0].startswith("loop"))
        n_total  = sum(1 for r in rows if r[-1] == "✓")
        print(f"{mode_label:36s} {n_simple}/{d_simple:<4d}   {n_chain}/{d_chain:<2d}    {n_loop}/{d_loop:<2d}    {n_total}/{len(rows):<2d}")
