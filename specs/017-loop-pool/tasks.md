# Tasks: Loop as a conserved pool

**Input**: `specs/017-loop-pool/{spec,plan,research,data-model}.md`
**Tests**: requested by the spec (assertions in a new file; the repo's rule).

## Phase 1: Setup

- [ ] T001 Create `test/test_loop_pool.jl` with the fixture helpers (`N`, `E`, `mk`, `P`), reusing `posloop()`'s shape from `test/test_loop_elasticity.jl`; one outer `@testset "loop pool"` so failures accumulate.

## Phase 2: Foundational (blocks all stories)

- [ ] T002 Extend the `DS_SCC_METHOD` allow-list in `src/solvers/steady_state.jl` (line ~651) to `fixed_point | minimize | pool | pool_parity`; keep the error message naming all four; add `pooling = scc_method in ("pool","pool_parity")`, `parity = scc_method == "pool_parity"`.
- [ ] T003 Make `solve_scc_ordered!` return a third value `(pooled, iterated, fallback_negative, fallback_inconsistent, pooled_nodes)` and update the call site (`steady_state.jl` line ~925) to store them in `SolverResult.diagnostics` as `scc_method`, `scc_pooled`, `scc_iterated`, `scc_fallback_negative`, `scc_fallback_inconsistent`, `scc_pooled_nodes`; the existing iterated path counts every cyclic component as `iterated`.

## Phase 3: US2 — default untouched, fallback never silent (P1) + US1 — pool on positive loops (P1) [MVP]

**Goal**: under the default nothing changes; under `pool` a positive loop reads its external supply.
**Independent test**: `test/test_loop_pool.jl` scenarios 1–5 of US1 and 1–3 of US2 pass; the eight existing assertion files pass unchanged.

- [ ] T004 [US1] In the per-component branch of `solve_scc_ordered!` (line ~761), before the iteration: if `pooling` and the component has no internal negative edge (`comp_neg_frac[c] == 0`, R2), build `xpool` (R1: copy of `x`, non-pinned member nodes set to `baseline_vec[i]`), compute `e_r` for every `r in rs`, collect entries (`e_r != 1`), sort entry folds by value, multiply with 0-absorb and `clamp(·, 0, 100)`, write `x[m] = clamp(baseline_vec[m] * pool, 0, 1)` for every non-pinned member node, bump `pooled` and `pooled_nodes`, and skip the iteration.
- [ ] T005 [US2] Under `pool`, a component with `comp_neg_frac[c] > 0` bumps `fallback_negative` and takes the existing iteration path unchanged.
- [ ] T006 [US1] Tests in `test/test_loop_pool.jl`: `posloop` at U=2 → A and B at fold 2 ± 1e-9 (both modes); U=0 → both 0; four co-required entries (1,1,3,2) → fold 6; two OR-alternative entries (fold 3 and 1) → the OR combination (mean 2), not 3; a closed loop with no external input → baseline.
- [ ] T007 [US2] Tests: default (`DS_SCC_METHOD` unset) equals `fixed_point` bit-for-bit on `posloop`, the parity fixture and `cyclic_net()` from `test/test_solver_determinism.jl`; a three-component network (two positive, one with an internal inhibitor) reports `scc_pooled == 2, scc_iterated == 1, scc_fallback_negative == 1` under `pool`; `DS_SCC_METHOD=pool_parityy` throws `ArgumentError`.
- [ ] T008 [US2] Relabel invariance: reuse `relabel` + `RELABELLINGS` from `test/test_solver_determinism.jl` (extract to a small shared include or duplicate the 10-line helper) and assert bit-identical activities under `pool` and `pool_parity` on every fixture; assert edge-order invariance by reversing edges.
- [ ] T009 [US2] Additive API field: `src/api/server.jl` solve response gains `"scc" => Dict(...)` from diagnostics; document in `docs/API.md` per `specs/017-loop-pool/contracts/scc-diagnostics.md`; add `scc_pooled`/`scc_iterated` to `solve_provenance` in `cli/deltasignal.jl`; extend `test/test_api_errors.jl` or `test/test_cli_observations.jl` with one assertion that the field is present and integer-valued.

## Phase 4: US3 — internal negatives by parity (P2)

**Goal**: `pool_parity` reads members downstream of an odd number of internal negatives with the inverse fold; inconsistent signs fall back.
**Independent test**: the X → M, M ⊣ T fixture reads T = 2x at X = 0.5 under `pool_parity` and falls back under `pool`.

- [ ] T010 [US3] Parity BFS in `solve_scc_ordered!` (R3): internal signed adjacency over member nodes from `rs` (activators +1, inhibitors/depleters −1, only indices with `comp_id == c`); for each entry node BFS assigning signs; conflict ⇒ `fallback_inconsistent += 1` and iterate; else member value `baseline[m] * prod_e fold_e^(s_e[m])` with the entry product taken in sorted order.
- [ ] T011 [US3] Tests: X → r1 → M, M ⊣ r2, r2 → T, T → r3 → M (so M and T share a component) with X = 0.5: `pool_parity` gives M ≈ 0.5, T ≈ 2 (± 1e-6) and `scc_pooled == 1`; `pool` gives `scc_fallback_negative == 1` and activities equal to `fixed_point`'s; an odd negative cycle (A ⊣ B, B → A) with an external drive reports `scc_fallback_inconsistent == 1` under both modes; a two-entry parity fixture (one entry at 0.5 upstream of the inhibitor, one at 2 downstream) reads T = 2 × 2 = 4.

## Phase 5: US4 — measurement, prediction first (P2)

- [ ] T012 [US4] Write the pre-registration in `specs/017-loop-pool/research.md` §Results (predictions P1–P5 with the decision rule from SC-004, the arms, catalog `cat_os`, baseline `ab_onesided.tsv`, both axes) and **commit it** before any arm runs.
- [ ] T013 [US4] Arm scripts from the `floor_ab.sh` pattern: `-e DS_SCC_METHOD=pool` and `=pool_parity` on `cat_os`, production defaults; the `verified:` line echoes `METHOD`; a probe solve prints the `scc` object and the arm aborts if `pooled == 0`; curator arms then experimental arms (`--ground-truth experimental`).
- [ ] T014 [US4] Score: `holdout_report.py` (held-out/tuning, concentration columns), `el_analysis.py` (false change, per-pathway net), McNemar p, the 102 TP53 AKT1/AKT2-KO cases directly, relabel churn on one cyclic pathway (SC-002); record each prediction held/failed in `research.md` §Results and commit.

## Phase 6: Polish

- [ ] T015 Differential check: default env on the pre-feature tree vs this tree on the six cyclic fixtures + `examples/` (the reviewer's `diffsuite.jl`); record "empty diff" in `research.md`.
- [ ] T016 CLAUDE.md: add `DS_SCC_METHOD` values and the new test file to the assertion table; `docs/RESULTS.md` only if a variant is adopted.

## Dependencies

T001 → T002 → T003 → {T004, T005} → {T006, T007, T008, T009} → T010 → T011 → T012 → T013 → T014 → {T015, T016}.
US1 and US2 are one increment (MVP); US3 depends on the MVP; US4 depends on US3 (both variants are measured together).

## Parallel opportunities

T006/T007/T008/T009 are independent files/sections once T004–T005 land; T015 and T016 are independent.
