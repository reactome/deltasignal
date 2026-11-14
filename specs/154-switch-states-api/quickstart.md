# specs/154 quickstart: validating switch-state serving

## Prerequisites

The dev API on the served catalog (`docker compose -f docker-compose.dev.yml
up -d julia-api`; health reports the current build).

## 1. Unit fixtures (FR-007)

```bash
docker compose -f docker-compose.dev.yml run --rm julia-api \
  julia --project=/app /app/test/test_switch_states.jl
```

Expected: the bistable fixture finds both points from both basins and selects
correctly under a pin to each side; the monostable control claims nothing;
the floor-rejection case serves `switch: false`; API field presence/absence
per the contract.

## 2. Parity replays (FR-001, the additive proof)

Replay the specs/141 readout dump and the specs/142 rails dump through the
feature branch (the specs/152 g4 pattern: dev container, RW gate mount,
`compare_141.py` / `compare_142.py`). Expected: zero differences.

## 3. The live scenario (US1)

HRAS over-expression in RAF (`R-HSA-5673001`), pins from `/api/pins`, solve,
readouts with the returned `switch_components` passed back. Expected: the
railed MAPK readouts serve `switch: true`, `state: "high"`,
`other_state_fold` near the low state; folds/calls/railed unchanged from the
recorded values.

## 4. The coverage report (US3, the serving bar)

```bash
# detection pass over the benchmark record (dev container, gate mount)
docker compose -f docker-compose.dev.yml run --rm \
  -v $PWD/specs/154-switch-states-api/report:/gate julia-api \
  julia --project=/app /app/bench/analysis/154/dump_switches.jl /gate
python3 bench/analysis/154/report_154.py specs/154-switch-states-api/report
```

Expected: the rails-coverage share printed; the RAF and PTEN/MKRN1 checks
PASS (the bar); the label-invariance re-run on the identical-content rebuild
(`20261113-1235_baf608b_hgncfix97`) flags the same set.
