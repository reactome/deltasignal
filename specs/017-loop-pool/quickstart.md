# Quickstart: validating the loop pool

Prerequisites: the dev image `deltasignal-julia-api`, the `deltasignal_julia-cache`
volume, and the production catalog mounted as in `bench/` scripts.

```bash
REPO=$(git rev-parse --show-toplevel)
run() { docker run --rm -v $REPO/src:/app/src -v $REPO/test:/app/test -v $REPO/examples:/app/examples \
  -v $REPO/Project.toml:/app/Project.toml -v deltasignal_julia-cache:/root/.julia deltasignal-julia-api "$@"; }

# 1. Unit acceptance (US1–US3): the new assertion file
run julia --project=/app /app/test/test_loop_pool.jl
#    expected: every @test passes; the file prints pooled/iterated counts per fixture

# 2. Default untouched (US2): the existing suites, unchanged
for f in test_propagator_invariants test_config_validation test_loop_elasticity test_solver_determinism \
         test_and_curves test_cycle_handling test_worked_example test_observation_pinning; do
  run julia --project=/app /app/test/$f.jl | tail -1; done

# 3. Differential check against the pre-feature tree (constitution IV):
#    solve the six cyclic fixtures + examples/ under default env on both trees and diff the reprs
#    (the reviewer's harness: scratchpad/diffsuite.jl) — expected: empty diff

# 4. Guard rail
DS_SCC_METHOD=pool_parityy run julia --project=/app -e 'using DeltaSignal' # -> ArgumentError naming the four values

# 5. Measurement (US4): pre-registration committed FIRST, then
#    bench arm scripts with -e DS_SCC_METHOD=pool and =pool_parity on cat_os, production defaults;
#    the arm's `verified:` line must echo METHOD and a probe solve's `scc.pooled` must be > 0;
#    holdout_report.py --cases ab_onesided.tsv --compare ab_pool*.tsv; experimental axis likewise;
#    TP53 AKT1/AKT2-KO cases scored directly (102 cases, 96 truth UP).
```
