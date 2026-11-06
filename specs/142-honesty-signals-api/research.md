# specs/142 research: decisions and the parity result

## Decisions

- **D1 — the flag is `score_124.railed`, nothing else.** Every recorded rails%
  figure (specs/107 behaviour checks, specs/124 comparisons) is that one
  function: `fold >= 0.99 * ceiling or fold <= 0.01`, boundaries railed. The
  engine's `readout_railed` is a line-for-line port; the gate holds them equal
  on every scored case fold of the frozen record.
- **D2 — stops served, not hardcoded.** `protocol.rail_low` / `rail_high` carry
  the resolved values (0.01 and 0.99 × ceiling, ceiling = 1 / `DS_BASELINE` as
  specs/141 resolves baseline), the same pattern as the served cutoffs.
- **D3 — `rail_low` stays 0.01 under any baseline.** It is a fraction of rest
  (1% of resting level), matching the bench's constant `0.01`; only the high
  stop follows the ceiling. Documented in the spec's assumptions and pinned by
  a unit test at baseline 0.001 (stops (0.01, 990.0)).
- **D4 — valued-but-valueless readouts carry `railed: null`,** mirroring
  `value`/`fold`/`call` exactly (basis `none`, abundance `trivial`/`no_copy`);
  resolution-only queries carry no `railed` key at all, as they carry no
  `value` key.
- **D5 — the gate runs offline on recorded folds.** `railed` is a pure
  function of the fold, so the frozen build's recorded `pred_ui` values are
  the complete input; no re-solve, no new catalog. The gate reads the frozen
  scoring run (build `20261026-1349_8b6527e`, solver `05dc867` — the
  RESULTS.json freeze record: curator 21502/24100, experimental 643/849), not
  the specs/141 dump, which is resolution-only and carries no folds (checked
  before registration).
- **D6 — `/api/solve` untouched.** It already serves `converged`, `residual`
  and the `nonconv_*` diagnostics (specs/066/073/076); the UI half of gap 3
  consumes those as served.

## Parity result (the FR-003 gate)

```
$ python3 bench/analysis/142/dump_rails.py ~/deltasignal-catalogs/builds/20261026-1349_8b6527e/results/05dc867
24356 valid case folds dumped (3922 railed) -> bench/analysis/142/bench_rails.json
$ docker compose -f docker-compose.dev.yml run --rm -v $PWD/bench/analysis/142:/gate julia-api \
    julia --project=/app /app/bench/analysis/142/rails_142.jl /gate/bench_rails.json /gate/engine_rails.json
24356 folds flagged (stops (0.01, 99.0)) -> /gate/engine_rails.json
$ python3 bench/analysis/142/compare_142.py
24356 case folds compared; agreeing: 24356 (3922 railed) / differences: 0
no differences — the gate passes
```

24,356 = every valid case of the frozen record (curator 23,511, experimental
845); 3,922 railed (curator 3,706 = 15.8%, experimental 216 = 25.6%) —
the flags agree case for case, and the engine's stops equal the bench's
ceiling convention. The bench side imports `railed` from
`bench/analysis/score_124.py`; the engine side calls the served
`readout_railed` with its env-resolved default stops.

## Unit coverage

- `test/test_readout_protocol.jl` (+15): railed high through a computed
  entity-mean reaching the ceiling, railed low through a knocked-out readout,
  mid-range un-railed, boundary folds exactly at both stops (railed, as the
  bench's `>=`/`<=` say), a baseline-read group folding to 1.0 (not railed),
  and the stops following `DS_BASELINE` (0.001 → (0.01, 990.0)).
- `test/test_api_errors.jl` (+9): `protocol.rail_low`/`rail_high` served;
  `railed` false on a mid fold, true at fold 100 and fold 0, `null` on basis
  `none` and abundance `trivial`, false on a scored abundance fold of 50, and
  absent on resolution-only queries.
