# Implementation Plan: Handle loops properly, and compare honestly

**Branch**: `004-loop-handling` | **Date**: 2026-09-10 | **Spec**: [spec.md](./spec.md)

**Input**: Feature specification from `/specs/004-loop-handling/spec.md`

## Summary

MP-BioPath's published networks are acyclic; ours are not. Four of the nine
published benchmark networks have no cycle at all and the largest component
anywhere is eleven nodes, against 836, 643, 465 and 443 in ours. The
365-versus-407 headline has therefore never compared like with like.

The plan attributes the deficit before trying to close it. First run the
DeltaSignal propagator on MP-BioPath's own networks — never done, and it
splits the 42-case gap into a propagator share and a network share. In
parallel, fix the reporting so no comparison is stated without the
acyclicity difference beside it. Only then intervene on loops, ordered by
how much curated structure each intervention destroys, with a
feedback-arc-set DAGification last as an unshippable upper bound.

Phase 0 resolved all three unknowns by measurement: the conjunction
convention in MP-BioPath's files (`0` = AND), the identifier mapping (sound
as-is; a 17.9% Reactome-wide divergence that is 100% clean on the entities
in play), and the artifact-versus-feedback discriminator (nodes per distinct
reaction stable id: 20–161 for artifacts, 3–9 for genuine components).

## Technical Context

**Language/Version**: Julia 1.10 (engine), Python 3.11 (benchmark harness)

**Primary Dependencies**: existing only — no new package. SCC machinery
(Tarjan) already exists in `src/solvers/steady_state.jl`; the benchmark
harness already loads networks, maps cases and scores baselines.

**Storage**: filesystem — LNG catalog directories, MP-BioPath's published
four-column TSVs, benchmark output directories.

**Testing**: Julia `@test` in `test/`; assertions must go in a file that
already has them (`test_config_validation.jl`, `test_and_curves.jl`,
`test_worked_example.jl`) or a new one.

**Target Platform**: Linux, local Julia server on a per-arm port.

**Project Type**: research engine plus offline benchmark harness.

**Performance Goals**: none. Each benchmark arm is minutes; the constraint is
memory, not speed — one arm at a time, kill each Julia server before the
next.

**Constraints**: every arm reads the same catalog directory; every arm
reports macro-F1, per-pathway net change, both-arms-converged count and
coverage change. Non-converged solves differ run to run because uuid4 ids set
sweep order, so a raw case delta is not evidence on its own.

**Scale/Scope**: ten pathways, 564 experimental and 3,914 curator scored
cases, 43,968 catalog edges, 4,561 cycle-resident nodes across 34 components.

## Constitution Check

*GATE: checked before Phase 0 and re-checked after Phase 1.*

| Principle | Status | Note |
|---|---|---|
| I. Processing is DeltaSignal's job | **GATE — binding** | Edge-type removal collapses 4,561 cycle-resident nodes to 494, but deleting upstream structure to make the graph easy is exactly what this principle forbids. Interventions that delete curated causality (catalyst, regulator, depletion edges) are admissible **only as diagnostics**. A shippable outcome is either DeltaSignal solving the cyclic structure correctly, or removal of something that is not curated causality. `diagram_bridge` is the one candidate that qualifies: it is a generated connectivity heuristic, not a curator assertion. |
| II. Measure, then claim | PASS | FR-006 requires macro-F1 with per-pathway breakdown and both-arms-converged count for every arm; McNemar per the workflow section. |
| III. Negative results are results | PASS | FR-008 and SC-004 require recording arms that reduce cyclicity without improving accuracy. R2 in research.md is already a recorded non-finding. |
| IV. New behaviour ships default-OFF | PASS | FR-005 requires each intervention selectable and defaulted to current behaviour. |
| V. Honest solver reporting | PASS | Coverage loss must be reported as coverage loss (FR-007), not absorbed into accuracy. |
| VI. API boundary is a trust boundary | N/A | No API surface change. |

**Post-Phase-1 re-check**: unchanged. The design keeps every structure-deleting
arm behind a flag and labels the DAGification arm an upper bound in its own
output, so principle I is respected by construction rather than by discipline.

**Justified deviation**: none required.

## Project Structure

### Documentation (this feature)

```text
specs/004-loop-handling/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 — R1..R5, all unknowns resolved
├── data-model.md        # Phase 1 — entities and classification rules
├── quickstart.md        # Phase 1 — reproduce every figure
├── checklists/
│   └── requirements.md
└── tasks.md             # Phase 2 (/speckit-tasks)
```

### Source Code

```text
bench/
├── benchmark_mpbiopath_cases.py     # gains: cyclicity classification per case,
│                                    #   acyclic/cyclic split reporting
└── analysis/
    ├── cycle_structure.py           # NEW — reproducible SCC report (FR-010)
    ├── mpbiopath_network_adapter.py # NEW — 4-column TSV -> DeltaSignal network
    └── loop_interventions.py        # NEW — produce intervention catalogs

src/                                 # unchanged unless the control implicates
                                     #   the propagator; then feature 005
test/
└── test_cycle_classification.jl     # NEW if classification lands in Julia
```

**Structure decision**: everything lands in `bench/`, not `src/`. This feature
is measurement and attribution; it changes no solver behaviour by default.
If the control shows the propagator is at fault, that is a finding handed to a
new feature, not scope absorbed here.

## Phase sequencing

**US1 first (P1, the control).** It is cheap, has never been run, and its
answer redirects everything after it. If DeltaSignal scores near 407 on
MP-BioPath's acyclic networks, the propagator is exonerated and loops/network
structure own the whole gap — which makes US3 the main event. If it scores
near 365, the networks are exonerated and this feature's remaining arms are
expected to be neutral; the finding is then that the propagator is the
problem and feature 003 becomes the priority.

**US2 in parallel (P1, no experimental dependency).** The honest-reporting
obligation stands whatever the experiments show, so it must not be sequenced
behind them.

**US3 after US1 (P2), ordered by destructiveness:**
1. `diagram_bridge` removal — 99.9% cycle-resident, previously measured
   accuracy-neutral, generated heuristic rather than curated causality. The
   only arm that could ship under principle I.
2. Recycling-artifact breaking — targets components with a high
   nodes-per-reaction ratio. Principled, but it deletes edges derived from
   curated reactions, so shipping requires the orientation evidence
   (`precedingEvent`) rather than the ratio heuristic alone.
3. Combined.

**US4 last (P3), diagnostic only.** Feedback-arc-set DAGification bounds the
prize. Never shippable; reported with its deleted-edge count.

## Risks

- **Coverage confound.** Breaking cycles can disconnect a readout, turning a
  wrong answer into an unscoreable one and inflating accuracy. FR-007 exists
  for this; every arm reports coverage change and it must be checked before
  reading the score.
- **Noise masquerading as effect.** 179 of 564 cases do not converge and
  their values depend on sweep order. A change of a few cases concentrated in
  TP53 with no both-arms-converged support is the known artifact signature,
  seen twice before.
- **The control may not be constructible faithfully.** If MP-BioPath's
  published networks turn out not to be what its published predictions were
  computed on, US1 is invalid. Spot-check before trusting it: recompute
  MP-BioPath's own predictions from its networks for a handful of cases.
- **Attractive wrong conclusion.** "Removing edges improved the score" will
  be true for some arm. Principle I is the guard, and the spec states it as a
  non-goal, but it will need saying again when the number appears.
