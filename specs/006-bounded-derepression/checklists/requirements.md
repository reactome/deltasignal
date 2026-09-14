# Specification Quality Checklist: Bound de-repression on purpose, not by accident

**Purpose**: Validate specification completeness and quality before proceeding to planning
**Created**: 2026-09-14
**Feature**: [spec.md](../spec.md)

## Content Quality

- [x] No implementation details (languages, frameworks, APIs)
- [x] Focused on user value and business needs
- [x] Written for non-technical stakeholders
- [x] All mandatory sections completed

## Requirement Completeness

- [x] No [NEEDS CLARIFICATION] markers remain
- [x] Requirements are testable and unambiguous
- [x] Success criteria are measurable
- [x] Success criteria are technology-agnostic (no implementation details)
- [x] All acceptance scenarios are defined
- [x] Edge cases are identified
- [x] Scope is clearly bounded
- [x] Dependencies and assumptions identified

## Feature Readiness

- [x] All functional requirements have clear acceptance criteria
- [x] User scenarios cover primary flows
- [x] Feature meets measurable outcomes defined in Success Criteria
- [x] No implementation details leak into specification

## Notes

Validation notes, iteration 2 (iteration 1 failed one item):

- **Iteration 1 failure — implementation details.** The first draft named
  `DS_INHIBITOR_EPS`, `DS_INHIBITOR_FLOOR`, `reaction_model.jl:1247` and the
  formula `(baseline + eps)/(x + eps)` throughout the requirements. Rewritten
  in domain terms — "the division guard", "the maximum de-repression", "the
  suppression bound" — with the code references left to the plan. The
  arithmetic facts (eleven-fold, 0.0110 vs 0.0198) are kept as *evidence*,
  since removing them would make the requirements unfalsifiable.
- **The first requirement is an attribution requirement, not a change.**
  FR-003 and FR-004 exist because the measured gain moved one control that
  governs two behaviours, so the reported cause is currently unknown. A spec
  that went straight to "bound de-repression" would be assuming the answer to
  the question this feature has to ask first.
- **No [NEEDS CLARIFICATION] markers.** Two candidates went to Assumptions:
  whether the bound applies per reaction or end to end (raised explicitly in
  FR-005 as something the design must state, since the observed 25× and 100×
  came from compounding), and whether the ~2× candidate is right (recorded as
  a starting point that may not survive, not a recommendation).
- **Constitution check.** Principle II drives FR-004 and FR-006 — a gain
  whose mechanism is unattributed is exactly the "measure, then claim"
  failure. Principle III drives FR-010, including the case where the result
  turns out to come from suppression and the headline story changes.
  Principle IV drives FR-008. FR-007 encodes the overfitting audit's own
  conclusion against the smaller case set that produced this candidate.
- **SC-006 guards a specific trap.** A de-repression bound near the
  change/no-change cutoff could move many cases for reasons unrelated to
  biology, which would look like a modelling win and be a threshold artifact.
