# Feature Specification: Combination perturbations as a validation axis

**Feature Branch**: `specs/159-combination-axis`

**Created**: 2026-11-16

**Status**: Draft

**Input**: Improvement backlog #256 (the 2026-11-13 approved list, item 6 — the last; items 1–5 shipped end to end).

## The problem, sharpened by the record

The backlog note says no benchmark has scored combinations. The record says otherwise in a useful way: the korkut design IS mostly drug pairs — 66 of its 78 drug sets are two-drug conditions, every pair's both singles are also in the design, and the recorded benchmark already solved each pair as simultaneous pins (the engine's one-observation-set semantics). What has NEVER been measured is whether the pair solve carries information BEYOND the single-drug solves: the recorded pair predictions have only ever been scored as more rows, never against a combination-naive null. The record holds 1,442 pair cases joinable to both singles' cases, of which 41 are combination-specific (the pair's measured call differs from BOTH singles') — a small but real class only a genuine combination model can get right.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - The registered measurement: does the pair solve beat combination-naive nulls? (Priority: P1)

Entirely offline on the recorded korkut TSV (no new solves; the record's pair predictions ARE the pair solve): for each joinable pair case, compare the RECORDED pair prediction's call against two registered combination-naive nulls computed from the RECORDED single predictions — (a) **product**: the pair's predicted fold is the product of the two singles' predicted folds (log-additivity, the no-interaction null), called at the standard cutoffs; (b) **extreme**: the single prediction with the larger |log2 fold| stands for the pair. The gate: the pair solve's accuracy against the measured calls exceeds the PRODUCT null's, with a 95% bootstrap interval over the 66 pairs excluding zero; balanced accuracy and the extreme null reported beside it. Secondary, reported with its (wide) interval whatever it says: on the 41 combination-specific cases, the pair solve's hit rate vs each null's.

**Why this priority**: it is the axis — a pair prediction that never beats multiplying the singles is API convenience, not a modelled combination.

**Independent Test**: the scorer is fixture-checked on computed known answers (a synthetic record where the pair call equals the product null exactly must show delta 0; one where pairs are right and nulls wrong must pass) before it reads the real TSV.

---

### User Story 2 - The record is read, never re-made (Priority: P1)

No solves run; no recorded figure moves; korkut stays FIT and nothing is fitted (both nulls are parameter-free and registered here). The analysis commits whatever it shows.

---

### User Story 3 - Serving: a measured scope line (Priority: P2)

If the gate passes, the served scope statement (specs/143) gains a validated-axis line stating what combinations were measured to do (and the combination-specific hit rate beside it); if it fails, the honest outcome is the recorded negative and no claim — the product keeps accepting simultaneous pins either way, as it always has, with no new accuracy implied.

---

### Edge Cases

- A pair case whose readout is untouched (fold exactly 1.0) in the pair but touched in a single (or vice versa): joinable, scored as recorded — reach differences between pair and singles are part of what a combination changes.
- Non-converged recorded rows: scored as recorded (the record's own convention), counted in the report.
- The product null's fold exceeding the UI ceiling: clamped at the model's own bounds (the null must live on the same scale), stated in the script.

## Requirements *(mandatory)*

- **FR-001 (offline, read-only)**: the measurement reads the recorded 98 korkut TSV; `require_role("korkut2015", "fit")` guards it; no solves.
- **FR-002 (registered nulls)**: product and extreme, exactly as US1 defines them, no parameters.
- **FR-003 (the gate)**: pair vs product-null accuracy, 95% bootstrap over pairs, 10k resamples, seed 159; balanced accuracy and the combination-specific class reported whatever they say.
- **FR-004 (serving only per the gate)**: a scope-statement line on pass; a recorded negative on fail; nothing else changes either way.
- **FR-005 (fixtures that compute)**: the scorer fixture-checked before the real read.

## Success Criteria *(mandatory)*

- **SC-001**: the measurement is on the record with intervals, the nulls' numbers beside the pair solve's, and the combination-specific class broken out.
- **SC-002**: any served claim cites the gate; a fail serves nothing new.
- **SC-003**: no recorded figure, endpoint, or protocol changes.

## Assumptions

- The recorded single predictions are the correct inputs to the nulls (same catalog, same protocol, same scoring population — they are rows of the same TSV).
- WA needs nothing: the panel already takes several perturbations.

## Non-goals

- No synergy modelling, no dose–response surfaces, no new endpoints.
- No new solves and no new machinery — the backlog note's premise, honoured literally.
