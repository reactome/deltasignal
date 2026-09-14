# Feature Specification: Bound de-repression on purpose, not by accident

**Feature Branch**: `006-bounded-derepression`

**Created**: 2026-09-14

**Status**: Draft

**Input**: The accuracy investigation in `specs/004-loop-handling/research.md`
R15–R16 identified de-repression as the single largest systematic error, with
a mechanism confirmed in the code and a fix direction confirmed on two
independent ground truths.

## Why this feature exists

### The model asserts something nobody chose

Inhibition is computed as a ratio of the inhibitor's baseline level to its
current level, floored by a small constant. That constant exists to stop a
division by zero. It also — as a side effect nobody decided — sets the
**maximum de-repression**: how far removing an inhibitor can raise its
target. At the current value that maximum is **eleven-fold from a single
inhibitor**, compounding along a path, which is where observed predictions of
25× and 100× baseline come from.

Nothing measured that eleven. It is the arithmetic consequence of a numerical
guard. The model separately has an explicit control over how far an inhibitor
may **suppress** its target, so the asymmetry is stark: suppression is
bounded by choice, de-repression is bounded by accident.

There is a second inconsistency. At the inhibitor's own baseline the formula
yields exactly "no effect", for every value of the guard. So the model says an
inhibitor sitting at its normal level exerts no inhibition at all — and then
says that removing that same inhibitor produces an eleven-fold increase, out
of a state it just described as uninhibited.

### It is the dominant error, and it is measured

On the 742-case experimental set, of the 142 errors on cases that have a
directed path, **85 (60%) are outright direction reversals**. Accuracy split
by which route signs exist: positive-only 0.896, ambiguous 0.758,
**negative-only 0.320**.

Across all knockout cases:

| prediction | cases | correct | accuracy |
|---|---|---|---|
| knockout → decrease | 241 | 213 | **0.884** |
| **knockout → increase (de-repression)** | **78** | **28** | **0.359** |

MP-BioPath scores **60 of those same 78**, and calls them an increase only 37
times against our 78. The truth on the 78 is 28 increase, 22 no change, 28
decrease — so "never de-repress" is also wrong. **The direction is not
inverted; the magnitude is ungraded.**

This is roughly 32 of the 46-case actionable gap in one error mode, and it
spans six pathways, so it is not an artifact of the two pathways that
dominate the case set.

### The fix direction is confirmed, the attribution is not

Weakening the guard from its current value to one giving 2× maximum
de-repression gives, on identical cases:

| axis | before | after |
|---|---|---|
| experimental (742) | 491 correct, macro-F1 0.5771 | **506, 0.5986** |
| curator (4,346) | 3,016, macro-F1 0.6711 | **3,075, 0.6867** |

Held-out gain exceeds development gain, which is the generalisation
direction. The targeted subgroup behaves exactly as diagnosed: spurious
knockout-to-increase calls fall 78 → 68 while the correct ones hold at 28.

**But the guard changes two things at once.** Raising it weakens
de-repression (11× → 2×) *and* weakens suppression (the strongest an
inhibitor can suppress goes from 0.0110 to 0.0198 of baseline). The measured
gain **cannot be attributed** to either. Resolving that is the first
requirement here, not a footnote — without it we would ship a result we
cannot explain.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Know which half of the change did the work (Priority: P1)

The measured gain came from moving one control that governs two behaviours.
Before anything is adopted, those behaviours must be separately adjustable
and separately measured, so the reported cause matches the actual cause.

**Why this priority**: it is the difference between a finding and a
coincidence. Every downstream decision depends on knowing which effect
mattered, and adopting a default whose mechanism is unattributed is the kind
of result this project has had to retract before.

**Independent Test**: de-repression strength and suppression strength can be
varied independently, and each is benchmarked with the other held fixed.

**Acceptance Scenarios**:

1. **Given** the maximum de-repression is changed and suppression strength is
   held fixed, **When** the benchmark runs, **Then** the change in score is
   attributable to de-repression alone.
2. **Given** suppression strength is changed and maximum de-repression held
   fixed, **Then** likewise for suppression.
3. **Given** both are measured, **When** results are reported, **Then** the
   reported explanation names which effect produced the gain, or states
   plainly that both contribute and in what proportion.

---

### User Story 2 — De-repression is a stated assumption (Priority: P1)

Whatever the maximum de-repression is, it should be something a reader can
find, question and change — not an arithmetic consequence of a division
guard. The guard should do only what its name says.

**Why this priority**: this is the actual defect. Even if the measured gain
turns out to come from the suppression side, the model would still be
asserting an unexamined eleven-fold de-repression, and the next person would
still have no way to see it.

**Independent Test**: the maximum de-repression can be read from
configuration, and changing the division guard across a wide range does not
change it.

**Acceptance Scenarios**:

1. **Given** the configuration, **When** a reader looks for how far removing
   an inhibitor can raise its target, **Then** a single named control states
   it, with its default and the reason for that default recorded.
2. **Given** the division guard is varied across several orders of magnitude,
   **When** predictions are compared, **Then** they do not change except
   where a genuine division by zero was being prevented.
3. **Given** an inhibitor is removed entirely, **Then** the resulting
   increase does not exceed the stated maximum, at a single reaction or
   compounded along a path.

---

### User Story 3 — Decide the default on evidence from a wider sample (Priority: P2)

The sweep that produced the candidate value was run on the same case set the
result is reported against, 60% of which is two pathways. The default should
be chosen against the larger catalog.

**Why this priority**: it guards the decision rather than enabling it, and it
depends on User Stories 1 and 2 existing first.

**Independent Test**: the chosen default is justified by results on the
larger catalog, and the smaller set is reported as a secondary check.

**Acceptance Scenarios**:

1. **Given** candidate values, **When** they are evaluated, **Then** the
   evaluation covers the larger catalog and both ground truths.
2. **Given** a value is adopted, **Then** the record states what it was
   chosen against and what it was not.
3. **Given** a candidate improves the smaller set but not the larger,
   **Then** it is not adopted, and the discrepancy is recorded.

---

### Edge Cases

- **An inhibitor removed when its target has no activator supply.** Bounding
  de-repression by a fold-change still multiplies whatever is there; if the
  target is already at zero, no bound makes it rise. Behaviour must be
  defined rather than emergent.
- **Several inhibitors removed on one path.** Individually bounded increases
  still compound. Whether the bound is per-reaction or end-to-end is a
  modelling decision and must be stated, because the observed 25× and 100×
  predictions came from compounding, not from a single step.
- **An inhibitor raised rather than removed.** The same control governs the
  suppression direction; changing one must not silently change the other,
  which is the defect being fixed.
- **The truth is genuinely an increase.** 28 of the 78 cases are real
  de-repression. A bound tight enough to eliminate all of them would trade
  one error for another; the measure must show both.
- **Interaction with the classification cutoffs.** A bound near the
  change/no-change boundary could move many cases at once for a reason
  unrelated to biology. Sensitivity to that must be checked, not assumed.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The maximum amount by which removing an inhibitor can raise its
  target MUST be expressed as a single named configuration control, with a
  documented default and a recorded reason for that default.
- **FR-002**: The control that prevents division by zero MUST NOT determine
  the maximum de-repression. Varying it across several orders of magnitude
  MUST NOT change predictions other than where a division by zero would
  otherwise occur.
- **FR-003**: De-repression strength and suppression strength MUST be
  independently adjustable, and each MUST be benchmarkable with the other
  held fixed.
- **FR-004**: The reported explanation of any measured gain MUST name which
  effect produced it, or state the proportion attributable to each.
- **FR-005**: The bound MUST hold when several inhibitors are removed along a
  path, and the specification MUST state whether it applies per reaction or
  end to end.
- **FR-006**: Every arm MUST be reported with macro-F1, a per-pathway
  breakdown, the count of changed predictions where both arms converged, and
  the coverage delta.
- **FR-007**: A default MUST NOT be adopted on the strength of the 742-case
  set alone; adoption requires results on the larger catalog against both
  ground truths.
- **FR-008**: Any change MUST default to current behaviour until measured.
- **FR-009**: The effect on the targeted subgroup — knockouts predicted to
  increase — MUST be reported separately from the overall score, showing both
  how many wrong calls were removed and how many correct ones were lost.
- **FR-010**: Negative and neutral results MUST be recorded with their
  numbers, including any arm showing the gain came from suppression rather
  than de-repression.

### Key Entities

- **Maximum de-repression**: how far removing an inhibitor may raise its
  target. Currently eleven-fold and implicit; to become explicit.
- **Division guard**: the constant preventing division by zero. Currently
  also determining the above; to do only its stated job.
- **Suppression bound**: how far an inhibitor may lower its target. Already
  explicit; the model for what de-repression should look like.
- **Targeted subgroup**: knockout cases predicted to increase — 78 cases at
  0.359 accuracy against MP-BioPath's 0.769 on the same cases.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The measured gain is attributed — the report states whether it
  came from weaker de-repression, weaker suppression, or both, with numbers.
- **SC-002**: Maximum de-repression is readable from configuration, and
  varying the division guard over at least three orders of magnitude changes
  no prediction.
- **SC-003**: Accuracy on the targeted subgroup improves from 0.359, with the
  count of correct de-repression calls retained reported alongside.
- **SC-004**: Any adopted default improves macro-F1 on the larger catalog
  against both ground truths, or is not adopted.
- **SC-005**: Results are reported for the larger catalog as primary and the
  742-case set as secondary, with both stated.
- **SC-006**: The relationship between the bound and the classification
  cutoffs is measured, so a gain cannot be a threshold artifact in disguise.
- **SC-007**: Every arm tried is recorded with its numbers, including those
  that fail.

## Assumptions

- The divide form of inhibition is settled and not under review here; only
  its bounds are.
- Baseline activity is uniform across nodes, so a fold-change bound has the
  same meaning everywhere.
- The larger catalog is the 92-pathway build, which now carries the
  entity-to-node resolution exports and can therefore score set-valued
  readouts on the same basis as the smaller set.
- The two ground truths remain the curator and experimental sets. The
  experimental one carries a known ~19% subset with no directed path, which
  this feature does not address and should not be expected to move.
- The ~2× candidate is a starting point from a sweep on the smaller set, not
  a recommendation. It may well not survive the wider evaluation, and the
  record should say so if it does not.
- Some of the measured gain is within the known run-to-run variation from
  non-converged solves; the attribution work is expected to narrow that, not
  eliminate it.
