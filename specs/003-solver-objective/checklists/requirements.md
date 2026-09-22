# Specification Quality Checklist: Solve the objective the model specifies

**Purpose**: Validate specification completeness and quality before planning
**Created**: 2026-09-10
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

## Notes

Deliberate choices worth recording:

- The spec avoids naming the solver technique, the optimiser, or the objective
  formula in the requirements, stating them as observable behaviour instead
  ("returns a result", "reports how far from self-consistency"). The formula
  belongs in the plan.
- **SC1 says exceed, not match.** Matching a comparator is a floor. The lead's
  framing was explicit and the criterion reflects it.
- **SC5 protects the stable subset.** The temptation with a large solver
  change is to celebrate fixing 122 broken cases while quietly regressing the
  442 working ones. The criterion forbids that trade.
- **SC7 bounds runtime.** A minimisation can be arbitrarily slower than an
  iteration; without a bound, "correct" could mean unusable.
- FR5 offers a genuine either/or: make the weights work, or stop reporting
  them. Both are honest; only the current state is not.
