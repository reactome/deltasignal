# Specification Quality Checklist: Metabolic cofactors are participants, not conduits

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

## Constitution Check

- [x] **I. Processing is DeltaSignal's job** — this is the principle the
      feature exists to honour. The generator-side filter was abandoned for
      exactly this reason and the negative result is recorded.
- [x] **II. Measure, then claim** — paired A/B on one shared catalog build,
      macro-F1 reported alongside accuracy, both-arms-converged count stated,
      per-pathway attribution given including the regression.
- [x] **III. Negative results are results** — two negative results recorded
      permanently (node deletion −84, `SimpleEntity` class −163), plus the
      neutral ten-pathway reading and three defects found in my own work.

## Notes

Two judgement calls, made rather than escalated:

1. **`inert` as the default, not opt-in.** It is the modelling position the
   feature is for, and it measures positive on the only evidence that has held
   up at scale. Reversible with one environment variable.
2. **A hand-curated list rather than a rule.** Every rule tried (entity class,
   diagram multiplicity, degree) deleted second messengers along with the
   cofactors. The list is release-pinned and expected to be refreshed rather
   than derived.

One residual, recorded and deliberately not fixed: `GPVI-mediated_activation_
cascade` regresses −12. Patching it would be a pathway-specific rule, which the
non-goals forbid.
