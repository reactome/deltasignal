# Tasks: Long-range upregulation propagation

**Input**: Design documents from `/specs/002-upregulation-propagation/`
**Prerequisites**: spec.md, plan.md, research.md

**Tests**: Included. The change alters default predictions for every user of
the tool, and three parameters have already each been wrong in a different
region of the range while looking correct near baseline. Curve behaviour is
pinned by assertion so the next reader cannot repeat that.

**Organization**: Grouped by user story from spec.md.

## Format: `[ID] [P?] [Story] Description`

---

## Phase 1: Setup

- [ ] T001 Confirm a single reproducible catalog build exists for all ten
      evaluation pathways at Release97 on current LNG `main`, and record its
      path in `specs/002-upregulation-propagation/quickstart.md`. Every arm in
      this feature must share one build — see research.md R1.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: pin the curve behaviour before changing any default, so a
regression is caught by assertion rather than by benchmark archaeology.

**⚠️ CRITICAL**: No default may change until T003 passes.

- [ ] T002 [P] Add `test/test_and_curves.jl` asserting the AND design intent
      at the magnitudes signals actually reach, not only near baseline:
      `[0.5,0.5]→0.25`, `[1,1]→1`, `[2,0.5]→1`, `[2,2]→4` (all ±0.01);
      `[50]→50`, `[100]→100`, `[10,10]→100` (±0.1, i.e. constrained to the
      0–100 range rather than 10,000). Assert for `hill_sat` at
      `DS_HILL_SAT_EPS=1e-5`. Include the currently-failing `hill_log` values
      (`[50]→41.4`, `[100]→74.1`) as documentation of why the default changes.
- [ ] T003 [P] Add to `test/test_and_curves.jl` a low-end assertion that
      catches the epsilon defect: at `DS_HILL_SAT_EPS=1e-3`,
      `[0.25,0.25]` yields 0.09 rather than 0.0625, and at `1e-5` it yields
      0.0625 ±0.001. This is the bug that produced a systematically
      upward-biased model (research.md R6).
- [ ] T004 Verify both tests FAIL against current `main` defaults before any
      change, and record the failure output in the task notes. A test that
      passes before the fix is not testing the fix.

**Checkpoint**: curve behaviour pinned; defaults may now change.

---

## Phase 3: User Story 1 — A knocked-in gene raises its downstream readout (P1) 🎯 MVP

**Goal**: an increase applied at a root input survives to a terminal readout.

**Independent Test**: correct upregulation calls rise from 119/246 without
downregulation falling below 169/247, on one shared catalog build.

- [ ] T005 [US1] Change the `DS_HILL_SAT_EPS` default from `0.001` to `1e-5`
      in `src/core/reaction_model.jl`. This is a defect fix independent of the
      rest of the feature: at the current default anyone selecting `hill_sat`
      gets a model that predicts UP on 317 of 564 cases against 246 actually
      UP, because the epsilon exceeds the internal values where knockouts live
      (~0.0006) and acts as a floor.
- [ ] T006 [US1] Change the `DS_AND_MODE` default from `hill_log` to
      `hill_sat` in `src/core/reaction_model.jl`, with a comment recording
      that `hill_log` tanh-squashes the summed log-fold at `z_max=10` and
      `e^10` lies far outside the UI range, so compression is active
      throughout the operating range rather than only at the ceiling.
- [ ] T007 [US1] Change the `DS_ASSEMBLY_LIMITING` default from `1` to `0` in
      `src/core/reaction_model.jl`. Record that `min(elevated, baseline)` is
      exactly baseline, so a complex transmitted scarcity perfectly and
      blocked abundance completely, and that the rule was adopted on a curator
      win measured when perturbations barely had to cross a complex.
- [ ] T008 [US1] Update the Solver Configuration section of `CLAUDE.md` to the
      new defaults, replacing the `DS_AND_MODE=hill_log`,
      `DS_ASSEMBLY_LIMITING=1` line. State the design intent explicitly — AND
      multiplies fold-changes, constrained to 0–100, near-exact multiplication
      around baseline — so the next reader does not have to re-derive it.
- [ ] T009 [US1] Re-run the experimental benchmark on the shared catalog build
      and confirm ≥365/564, macro-F1 ≥0.578, UP ≥168/246, DOWN ≥169/247.
- [ ] T010 [US1] Confirm the result is not a prediction bias: median
      `predicted_ui` must remain 1.000, and predicted-UP must move toward 246
      rather than past it. This check is mandatory — its absence is what
      caused a bias to be reported as a win (research.md R6).

**Checkpoint**: SC3, SC4, SC5, SC6 evaluable; upregulation restored.

---

## Phase 4: User Story 2 — Results are trustworthy against a null model (P2)

**Goal**: DeltaSignal beats the structural baseline, and the curator win the
clamp was adopted for is not silently discarded.

**Independent Test**: DeltaSignal exceeds 393/564 on identical cases; the
curator comparison is measured rather than assumed.

- [ ] T011 [US2] Run the same arm against **curator** ground truth
      (`--ground-truth curator`) on the shared catalog build, both with and
      without the new defaults. The clamp was adopted for a +0.9pp curator
      gain (DSB +185, Interferon α/β −75); dropping it without re-measuring
      that axis trades a known win for an unmeasured one.
- [ ] T012 [US2] If the curator axis regresses, record the trade-off with its
      numbers in `research.md` and raise it for a decision rather than
      resolving it unilaterally. A change that helps experimental and hurts
      curator is a judgement call about which ground truth the tool serves.
- [ ] T013 [US2] Report DeltaSignal against the `shortest_signed_path`
      baseline (393/564) and MP-BioPath (407/564) on identical scored cases,
      with per-pathway attribution and a both-arms-converged count per FR5.

**Checkpoint**: SC2 evaluable; the curator axis is measured, not assumed.

---

## Phase 5: Polish & Cross-Cutting

- [ ] T014 [P] Record the measured outcome in
      `specs/002-upregulation-propagation/research.md`, including any negative
      or neutral result, per constitution principle III.
- [ ] T015 [P] Write `specs/002-upregulation-propagation/quickstart.md` with
      the exact commands to reproduce every figure in the spec, including the
      catalog build and each benchmark arm.
- [ ] T016 Note in `research.md` that non-convergence rises from 122 to 179
      under the new defaults — larger signals reach the cyclic components —
      so feature 003 becomes more important, not less. Do not attempt to fix
      it here; it is out of scope.
- [ ] T017 Re-check the residual gap to MP-BioPath after this feature lands.
      +31 closes part of 73; the remainder (42 wrong-direction upregulation
      cases and whatever else) is the next feature, not a rounding error to
      be absorbed into a success claim.

---

## Dependencies & Execution Order

- **Phase 1 (Setup)**: no dependencies.
- **Phase 2 (Foundational)**: depends on Phase 1. **Blocks all user stories** —
  T004 must show the tests failing before any default changes.
- **Phase 3 (US1)**: depends on Phase 2.
- **Phase 4 (US2)**: depends on Phase 3 — the curator run must exercise the
  new defaults, so it cannot precede them.
- **Phase 5**: depends on Phases 3 and 4.

### Within User Story 1

T005, T006 and T007 are three independent default changes but touch the same
file, so they are **not** parallel. Sequence them and re-run T009 after each,
to attribute the gain per change rather than only in aggregate — the separate
effects (+5, +21, +31 combined) are what showed the two mechanisms were
masking each other.

### Parallel Opportunities

- T002 and T003 (both new assertions in the same new file — write together).
- T014 and T015 (different files).

---

## Implementation Strategy

**MVP is User Story 1.** It delivers the measured +31 cases and is
independently testable. Stop and validate there.

**Do not skip Phase 2.** Every defect in this feature was invisible to a check
performed near baseline: the clamp, the `hill_log` compression, and the
`hill_sat` epsilon each look correct at the 0.85/1.15 cutoffs and each fail
somewhere else in the range. The assertions in T002–T003 exist specifically to
make that class of error impossible to reintroduce silently.

**T010 and T011 are not optional polish.** T010 is the check whose absence
turned a bias into a reported win. T011 is the check that stops us trading an
unmeasured curator regression for a measured experimental gain.
