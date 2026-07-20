# attic/ — quarantined dead code

These files were part of `src/` and `test/` but a call-graph audit (2026-07)
confirmed **none of them are reachable from the live parse/solve path** — not
from the CLI (`cli/deltasignal.jl` `parse`/`solve`), not from the HTTP API
(`src/api/server.jl` `/api/parse`, `/api/solve`). They were only ever exercised
by their own unit tests, which are also moved here.

They are kept (not deleted) so the work is easy to revive if the time-dynamic,
parameter-learning, or biological-realism modes are ever wired to a real entry
point. They are plain module-internal source: they assume they are `include`d
inside `module DeltaSignal` and are not standalone-runnable as-is.

## Why each is dead

| file | reason |
|---|---|
| `src/core/biological_realism.jl` | entry fns called only by other attic modules + `test/test_enhanced_biological_realism.jl` |
| `src/core/compartmentalization.jl` | callers only in `test_enhanced_biological_realism.jl` |
| `src/core/experimental_constraints.jl` | zero callers anywhere (not even a test) |
| `src/core/temporal_dynamics.jl` | callers only in `test_enhanced_biological_realism.jl` |
| `src/core/stochastic_effects.jl` | callers only in `test_enhanced_biological_realism.jl` |
| `src/core/feedback_enhancements.jl` | used only by `enhanced_steady_state.jl` (itself dead) + `test_enhanced_negative_feedback.jl` |
| `src/solvers/enhanced_steady_state.jl` | `solve_steady_state_enhanced` — CLI/API call plain `solve_steady_state`; only tests call this |
| `src/solvers/time_dynamic.jl` | `rollout_time_dynamic` — the CLI `rollout` command is a "not yet implemented" stub |
| `src/learning/parameter_learning.jl` | `learn_parameters` — the CLI `train` command is a "not yet implemented" stub |

## The live path (kept in src/)

`cli/deltasignal.jl` → `io/tsv_parser.jl` → `solvers/steady_state.jl` →
`core/reaction_model.jl` → `core/{aggregators,hill_functions,sensitivity}.jl`.
`core/sensitivity.jl` in particular is LIVE (`apply_sensitivity_transform`) —
it stayed in `src/`.

## Reviving

Move the file back under `src/`, re-add its `include(...)` and the relevant
`export`s in `src/DeltaSignal.jl`, and move its test back to `test/`.
