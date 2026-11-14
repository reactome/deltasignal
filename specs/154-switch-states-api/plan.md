# Implementation Plan: Serve bistable components as switches

**Branch**: `154-switch-states-api` | **Date**: 2026-11-13 | **Spec**: [spec.md](spec.md)

**Input**: Feature specification from `/specs/154-switch-states-api/spec.md`

## Summary

After the canonical solve (untouched), detect per cyclic component whether a
second stable fixed point exists by re-sweeping the component from a high and
a low basin with everything outside it frozen at canonical values; propagate a
found flip through the component's feed-forward downstream to get a sparse
consistent alternate; `/api/solve` serves the alternates sparsely, and
`/api/readouts` resolves each readout a second time through the identical
specs/141 machinery to serve `switch`/`state`/`other_state_fold`. Additive
only, parity-replayed; the coverage report (bench/analysis/154) is the serving
bar, anchored on the two traced cases (RAF specs/054, PTEN/MKRN1 specs/082).

## Technical Context

**Language/Version**: Julia 1.10 (engine), Python 3 (bench report scripts)

**Primary Dependencies**: existing solver (`src/solvers/steady_state.jl`),
reaction model (`compute_reaction_output_vec`), pool machinery
(`src/core/cycles.jl` `update_pools!`), readout machinery
(`src/api/readouts.jl`), HTTP handlers (`src/api/server.jl`)

**Storage**: none (stateless API; sparse alternates ride the solve response)

**Testing**: new `test/test_switch_states.jl` in `test/asserting_suites.txt`
+ CLAUDE.md count table; parity replays of the specs/141/142 dumps

**Target Platform**: the dev/prod Docker API, unchanged deployment

**Project Type**: engine + HTTP API (UI consumption is a later WA spec)

**Performance Goals**: detection adds at most 2 bounded component re-solves
per cyclic component, env-capped (`DS_SWITCH_MAX_SWEEPS`, default = the
specs/076 restart sweep budget, 250); hitting the cap = honest absence

**Constraints**: FR-001 byte-identical existing fields; label-invariance
(specs/013); no scoring change; no new catalog

**Scale/Scope**: 92 served pathways; the largest component is DSB Repair's
454 nodes — two extra bounded sweeps of one component are well under a second

## Constitution Check

- **I (processing is DeltaSignal's job)**: pure engine capability; no
  generator change. PASS.
- **II (measure, then claim)**: no accuracy claim is made at all; the
  registered coverage report (FR-006) runs before serving and reports
  whatever it finds. PASS.
- **III (negative results are results)**: if detection fails the traced-case
  bar, that is recorded in research.md and the feature does not serve. PASS.
- **IV (default-off for behaviour changes)**: the feature changes no
  behaviour — canonical predictions byte-identical (FR-001 parity proof);
  precedent: specs/142's `railed` shipped default-on because additive.
  `DS_SWITCH_DETECT=on` default, `off` escape. PASS with the parity proof as
  the justification.
- **V (honest reporting)**: honest absence is a named FR (FR-004). PASS.

Post-design re-check: unchanged — the design adds fields and a post-pass, no
existing path is modified except to RECORD (not alter) per-component
bookkeeping at the end of the solve.

## Phase 0 findings

See [research.md](research.md). Every NEEDS CLARIFICATION is resolved; the
key design decisions: the alternate is served SPARSELY per bistable component
(only nodes that differ beyond the floor), downstream feed-forward re-solved;
components are re-swept with the SAME mutated `IndexedReaction`s and the same
per-sweep sequence (reaction sweep + `update_pools!`), so every adopted rule
applies identically in both basins; pool-managed components are searched, not
excluded (RAF's component carries pools and is the serving bar).

## Phase 1 artifacts

- [data-model.md](data-model.md) — the SwitchComponent diagnostic, the sparse
  alternate, the readout switch fields and their presence rules.
- [contracts/switch-fields.md](contracts/switch-fields.md) — the additive
  /api/solve and /api/readouts contract changes, with presence/absence rules
  and the protocol block.
- [quickstart.md](quickstart.md) — end-to-end validation: unit fixtures, the
  HRAS/RAF live scenario, the parity replays, the coverage report.

## Sequence

1. **US1+US2 (MVP)**: solver post-pass + sparse alternates + solve/readouts
   fields + `test/test_switch_states.jl`, behind the FR-001 parity replays.
2. **US3**: `bench/analysis/154/` coverage report (the serving bar: both
   traced cases flagged with matching states; rails coverage share reported).
3. Serving PR: docs/API.md, CLAUDE.md table, scope statement's
   `switch_pathways` limitation gains the pointer to the served fields.
