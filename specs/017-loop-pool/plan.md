# Implementation Plan: Loop as a conserved pool

**Branch**: `feat/017-loop-pool` | **Date**: 2026-09-19 | **Spec**: [spec.md](spec.md)

**Input**: Feature specification from `specs/017-loop-pool/spec.md`

## Summary

Replace the per-component damped iteration with a pool rule for cyclic
components under two new `DS_SCC_METHOD` values: each member reaction's
external fold is computed with in-component inputs at baseline, the entry
folds multiply into one pool fold, and every non-pinned member reads its
baseline × that fold (`pool`; falls back on any internal negative edge) or ×
the parity-signed fold (`pool_parity`; falls back on inconsistent signs).
Pooled/iterated counts are reported in diagnostics and the API. Default
`fixed_point` is bit-identical. Both variants are measured on the production
catalog, both axes, predictions committed first.

## Technical Context

**Language/Version**: Julia 1.10 (LTS)
**Primary Dependencies**: none new; `solve_scc_ordered!`, `compute_reaction_output_vec`, Tarjan components — all existing
**Storage**: N/A
**Testing**: `test/test_loop_pool.jl` (new, assertions) + the eight existing assertion files; differential run vs the pre-feature tree
**Target Platform**: the `deltasignal-julia-api` container
**Project Type**: library + HTTP API
**Performance Goals**: no measurable change on the default path (the pool path does one extra evaluation per member reaction per component — cheaper than iteration)
**Constraints**: default-OFF; label- and edge-order-invariant; additive API change only
**Scale/Scope**: 251 cyclic components across 92 pathways; ~30 large ones

## Constitution Check

I ✓ (solver only) · II ✓ (wide set, both axes, p, held-out; pre-registration committed before the arm) · III ✓ (SC-004 names the negative outcome) · IV ✓ (default unchanged; differential check is a task) · V ✓ (fallback counts reported) · VI ✓ (additive `scc` object). No violations.

## Project Structure

### Documentation (this feature)

```text
specs/017-loop-pool/
├── plan.md
├── research.md          # R1–R5 resolved from code
├── data-model.md
├── quickstart.md
├── contracts/scc-diagnostics.md
└── tasks.md             # /speckit-tasks
```

### Source Code

```text
src/solvers/steady_state.jl   # solve_scc_ordered!: pool branch, parity BFS, counts; call site stores diagnostics
src/core/reaction_model.jl    # (no change needed: DS_SCC_METHOD is validated in steady_state.jl; extend that list)
src/api/server.jl             # additive "scc" object in the solve response
docs/API.md                   # document it
cli/deltasignal.jl            # provenance already records scc_method; add pooled counts
test/test_loop_pool.jl        # new assertion file
specs/017-loop-pool/research.md  # §Results: pre-registration, then scores
```

## Design (from research)

1. **Pool branch** in the per-component `else` (line ~761): if `pooling` and the component is eligible, compute `e_r` for all `rs` via `xpool` (R1), collect entries, sort their folds, multiply with 0-absorb and cap 100 (data-model), assign members (R3 for parity), bump counts, `continue` — no iteration.
2. **Eligibility**: `pool`: `comp_neg_frac[c] == 0` (R2) else `fallback_negative += 1` and iterate. `pool_parity`: parity BFS (R3); inconsistent ⇒ `fallback_inconsistent += 1` and iterate.
3. **Counts** returned as a third value; call site → diagnostics; API → `scc` object (R4).
4. **Validation**: extend the `DS_SCC_METHOD` allow-list at line 651; guard-rail test iterates it.
5. **Sequence**: US2+US1 (MVP: pool on positive loops, counts, default identical) → US3 (parity) → US4 (arms).

## Complexity Tracking

None: no new dependencies, no new modules, one additive API field.
