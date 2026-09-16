# Specification Quality Checklist: Handle cycles correctly in the propagator

**Purpose**: Validate specification completeness and quality before planning
**Created**: 2026-09-16
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
- [x] Success criteria are technology-agnostic
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

- [x] **I. Processing is DeltaSignal's job** — the non-goals forbid deleting
      edges or breaking cycles in the generator, following the cofactor
      precedent. Cycles are a processing concern.
- [x] **II. Measure, then claim** — FR-008 pins the wide curator set as the
      decision basis and explicitly bars the nine-pathway set, which is what
      caused the anti-collapse machinery to be abandoned.
- [x] **III. Negative results are results** — "Prior Findings This Feature Must
      Not Repeat" records three failed classification formulations with their
      numbers, the in-degree confound, and the four candidate discriminators
      that did not survive.

## Notes

**Deliberately unresolved, and the feature's first task.** The `4.045534594765421e-6`
residual is not explained. The spec states the observation and the eliminations
(multiple reactions per target, catalyst-supply freezing, iteration budget,
damping) without asserting a cause. Planning must start by explaining it, not by
choosing a fix.

**Two candidate directions, not yet chosen**, which is why this is Draft rather
than Planned:

1. *Classify then treat* — reaction-level typing into the four curator types,
   handling each differently. Stronger biologically, but ~98.5% of edges are
   curator-derived so there is little to sort on, and misclassification would
   silently damage real feedback.
2. *Fix the dynamics generically* — exclude the spurious all-zero root and
   repair the convergence measure, without classifying. Simpler, no new
   generator artifact, no misclassification risk; but recycling and causal
   feedback arguably *should* behave differently and this treats them alike.

The measurement in this spec does not settle which. Deciding is a planning
question and needs Adam's steer.

**Two judgement calls made rather than escalated:**

1. The all-zero collapse is treated as numerical, not biological. Stated by
   Adam directly: the math is satisfied by all-zero even when the loop's inputs
   are ≥ 1. An earlier reading of it as a mass-balance question was wrong.
2. US1 and US2 are both P1. The convergence defect could be deferred as
   "only reporting", but the reported non-convergence rate is the evidence
   anyone would use to size this work, so a wrong number there misdirects
   everything downstream.
