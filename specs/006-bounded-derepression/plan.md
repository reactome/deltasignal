# Implementation Plan: Bound de-repression on purpose, not by accident

**Branch**: `006-bounded-derepression` | **Date**: 2026-09-14 | **Spec**: [spec.md](./spec.md)

## Summary

The inhibitor branch of the reaction model lets removing a single inhibitor
raise its target elevenfold, because the constant that prevents division by
zero also sets the ceiling. Depletion edges one branch over already have an
explicit named bound for exactly this, so the fix is a consistency change
with a precedent in the same function rather than a new idea.

Phase 0 established that the guard moves **three** things at once — the
ceiling, the curve shape at intermediate values, and the suppression floor —
and that they push in opposite directions. That makes the measured +15/+59
gain unattributable as it stands, and makes a 2×2 of (guard, ceiling) the
first piece of work rather than the last.

## Technical Context

**Language/Version**: Julia 1.10 (`src/core/reaction_model.jl`); Python 3.11
benchmark harness.

**Primary Dependencies**: none new.

**Storage**: none. This is a per-solve configuration change.

**Testing**: `test/test_config_validation.jl` and `test/test_and_curves.jl`
are the two files that carry real assertions and are where new ones go. The
arithmetic here is exactly the kind that is cheap to pin and expensive to get
wrong silently.

**Target Platform**: Linux, local Julia server per benchmark arm.

**Project Type**: research engine plus offline benchmark.

**Performance Goals**: none — this is one `min` in an inner loop.

**Constraints**: default must reproduce current behaviour exactly until
measured (principle IV); every arm reports macro-F1, per-pathway net,
both-arms-converged count and coverage delta; adoption requires the
92-pathway catalog against both ground truths.

**Scale/Scope**: one function branch, two configuration controls, four
benchmark arms on two ground truths and two catalogs.

## Constitution Check

| Principle | Status | Note |
|---|---|---|
| I. Processing is DeltaSignal's job | PASS | Entirely in the propagator, which is where a modelling assumption belongs. Nothing upstream changes. |
| II. Measure, then claim | **GATE — and it shaped the plan** | The existing +15/+59 result has an unattributed mechanism. US1 exists to fix that before anything is adopted, and FR-004 forbids reporting a cause that has not been isolated. |
| III. Negative results are results | PASS | FR-010 requires recording the outcome where the gain turns out to come from suppression, which would overturn this feature's framing. That must be written down, not quietly reframed. |
| IV. New behaviour ships default-OFF | PASS | Arm A (ε 1e-3, ceiling 11) reproduces today's behaviour and is the default until the evidence lands. |
| V. Honest solver reporting | PASS | The ceiling is reported in the solve provenance alongside the other config, so a result can be traced to the assumption that produced it. |
| VI. API boundary is a trust boundary | PASS | The new control is a `DS_*` numeric and goes through the existing `_float_env` validation, which rejects a non-numeric rather than silently defaulting. |

**Post-Phase-1 re-check**: unchanged.

**Justified deviation**: none.

## Project Structure

```text
src/core/reaction_model.jl        # the inhibitor branch; new bound + guard split
test/test_config_validation.jl    # defaults, validation of the new control
test/test_and_curves.jl           # or a sibling: the arithmetic of the bound
bench/analysis/compare_set_rules.py   # already reports what the arms need
specs/006-bounded-derepression/
├── plan.md  spec.md  research.md  quickstart.md
├── contracts/inhibition-bounds.md
└── checklists/requirements.md
```

**Structure decision**: no new module. The change is a `min` and a config
field in a branch that already has the shape it needs, plus tests. Inventing
an abstraction here would obscure a three-line change.

## Phase sequencing

**US1 first — attribution.** Separate the controls, then run the 2×2 from
research R3:

| arm | ε | ceiling | isolates |
|---|---|---|---|
| A | 1e-3 | 11 | current behaviour |
| B | 1e-9 | 11 | ε's curve-shape and suppression effects |
| C | 1e-9 | 2 | the de-repression ceiling |
| D | 1e-2 | none | the original probe, for continuity |

The arms discriminate because A→B *strengthens* suppression while C weakens
de-repression. If the gain was a suppression effect, B is worse than A and C
does not recover it.

**US2 alongside** — the separation is the same edit. Its distinct deliverable
is the evidence that the guard no longer sets the ceiling: vary ε over three
orders of magnitude with the ceiling fixed and show predictions unchanged
(FR-002). That test is what makes the separation real rather than nominal.

**US3 last** — choose the default on the 92-pathway catalog, both ground
truths, with the 742-case set reported as secondary.

## Risks

- **The gain may be a suppression effect.** Then this feature's headline is
  wrong and the interesting change is elsewhere. FR-010 requires saying so.
  Framing this as a risk rather than an afterthought is the point of US1.
- **The result may be threshold-adjacent.** A ceiling near the 0.85/1.15
  cutoffs could move many cases for reasons unrelated to biology. SC-006
  requires measuring against the cutoffs, because a gain that vanishes when
  the cutoffs shift slightly is an artifact.
- **Noise.** Only 9 of 29 changed experimental predictions and 34 of 239
  curator ones converged in both arms. Per-arm both-converged counts are
  mandatory, and a small delta without both-converged support is not a
  result.
- **Tuning by another name.** The ceiling is a free parameter and this
  project has a documented history of fitting them to this case set. FR-007
  moves adoption to the 92-pathway catalog; the temptation to sweep finely on
  the 742 remains and should be resisted rather than managed.
- **`DS_DEPLETION_H_MAX` is the same assumption, unexamined.** Its default of
  10 has exactly the arguments against it that this feature makes. It is
  deliberately out of scope; changing it silently alongside would confound
  the attribution this feature exists to establish.
