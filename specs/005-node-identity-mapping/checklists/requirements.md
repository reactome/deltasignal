# Specification Quality Checklist: Every node maps, both ways, including on the diagram

**Purpose**: Validate specification completeness and quality before proceeding to planning
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

Validation notes, iteration 2 (iteration 1 failed two items; both fixed):

- **Iteration 1 failure — implementation detail.** The first draft named
  `nodes.csv` columns, `src/diagram_connectivity.py`, `hasMember`/
  `hasCandidate` and a proposed file layout throughout the requirements.
  Fixed: FR-001..FR-012 are now stated in domain terms (entity, glyph, node,
  resolution) and the file and code references are confined to the "Why this
  feature exists" section, where they are *evidence for the problem* rather
  than a prescribed solution. The plan is the right place for the artifact's
  shape.
- **Iteration 1 failure — untestable requirement.** "The mapping must be
  perfect" was the original phrasing and is not checkable. Fixed by splitting it
  into FR-001/FR-002 (both directions return labelled results), FR-006
  (absence must be declared with a reason) and FR-007 (the check must be
  demonstrated to fail), with SC-001/SC-002 giving the zero-unexplained-
  absences target. "Perfect" is now "no unexplained absence in either
  direction".
- **The spec contradicts its own Input, deliberately.** The input said "there
  really should be an id for each node to uniquely identify them but there
  isn't." There is one, it is measured, and the generator already reads and
  discards it. Stating this in the spec is necessary because it changes the
  design from "invent an identity" to "stop discarding one", and burying it
  would leave the plan solving the wrong problem.
- **No [NEEDS CLARIFICATION] markers.** Three candidates were resolved into
  Assumptions instead: glyph-id uniqueness beyond the two diagrams measured
  (checkable by the same measurement); glyph-id stability across releases
  (assumed NOT stable, which is the conservative direction and drives
  FR-005); and the set-combining rule, which US4 measures rather than
  asserts.
- **Constitution check.** Principle I places this upstream: representing
  what curators authored, including which glyph a node came from, is the
  generator's job, so the mapping is generated rather than reconstructed by
  consumers. Principle III drives SC-007 (losing rules recorded).
  Principle IV drives FR-010 (rules default to current behaviour).
  Principle V's spirit drives FR-006, FR-009 and FR-011 — a partial result
  reported as complete, and a declared-but-empty column, are the same
  dishonesty as a fake convergence report.
- **Denominator warning is load-bearing.** FR-012 and the Assumptions note
  that recovering 204 cases makes post-feature numbers incomparable to
  pre-feature ones. Without it the next report will show an accuracy change
  that is really a coverage change.
