# specs/154 tasks

## US1+US2 MVP (solver + API + tests, behind the parity gate)

- [ ] T1. Solver: `detect_switch_components` post-pass at the end of
  `solve_scc_ordered!` — per cyclic component (size > 1, not fully pinned,
  canonical residual below tolerance): re-sweep from high (1.0) and low (0.0)
  basin inits with pins held and outside frozen, same damped gauss-seidel +
  `update_pools!` sequence, `DS_SWITCH_MAX_SWEEPS` cap; accept a second point
  on convergence + distinctness (> `DS_SWITCH_FLOOR` = 1e-3 internal on any
  member); propagate the flip through downstream acyclic updates and
  downstream components' feed-forward recompute (single pass, values only);
  emit sparse alternates (nodes differing > floor) + diagnostics
  (`switch_detect`, `_searched`, `_found`, `_sweeps_capped`).
- [ ] T2. Env/config: `DS_SWITCH_DETECT` (on|off, default on),
  `DS_SWITCH_MAX_SWEEPS` (int, default 250), `DS_SWITCH_FLOOR` (float,
  default 1e-3) — validated like every DS_* flag (config validation suite).
- [ ] T3. `/api/solve`: serve `switch_components` (sparse, UI scale) when
  detection ran and the solve converged; `[]` vs absent per the contract.
- [ ] T4. `/api/readouts`: accept optional `switch_components` (validate:
  400 on malformed); value each resolved readout under canonical and each
  overlay; serve `switch`/`state`/`other_state_fold` by the call-cutoff rule;
  protocol block additions.
- [ ] T5. `test/test_switch_states.jl`: bistable fixture (two-node mutual
  activation above the knife-edge; both basins; selection under a pin to each
  side), monostable control, floor rejection, API presence/absence, env
  validation; register in `test/asserting_suites.txt` + CLAUDE.md table.
- [ ] T6. Parity replays (FR-001): specs/141 readout dump and specs/142 rails
  dump through the branch — zero differences required.

## US3 (the serving bar)

- [ ] T7. `bench/analysis/154/dump_switches.jl` + `report_154.py`: detection
  pass over the benchmark case list, join against the 3,922 railed folds;
  RAF + PTEN/MKRN1 traced-case checks; label-invariance re-run on the
  identical-content rebuild.

## Serving

- [ ] T8. docs/API.md blocks; scope statement `switch_pathways` limitation
  gains the served-fields pointer; PR; live smoke (quickstart §3).
