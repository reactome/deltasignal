# Feature Specification: Inverse queries

**Feature Branch**: `specs/158-inverse-queries`

**Created**: 2026-11-15

**Status**: Draft

**Input**: Improvement backlog #255 (the 2026-11-13 approved list, item 5; items 1–4 shipped end to end, item 4 on 2026-11-15).

## The problem

Users ask "what lowers Y?" as often as "what does X do?". The product only answers the second. Everything needed for the first already exists and is validated: the pin protocol enumerates what can be perturbed (specs/140), the forward solve is the measured engine, and the readout ladder values Y under each solve with every honesty field attached (specs/141–157). An inverse query is nothing but the validated forward solves run over candidate perturbations and read at one chosen readout — a product feature, not a new validation axis, and it must be built exactly that way so every served number inherits the record.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Which perturbations move this readout (Priority: P1)

A user names a readout (the same query shape /api/readouts takes: stable ids, optional members, optional abundance mode) and a wanted direction (down, up, or any). The engine enumerates the network's candidate perturbations — exactly the entities the pin protocol can pin (the validated root/root-cycle resolution, specs/140), one candidate per curated entity, each at the protocol's standard doses (knockout 0; over-expression 80) — runs the ORDINARY forward solve per candidate, values the chosen readout through the ORDINARY readout ladder, and returns the candidates whose call at the readout matches the wanted direction, ranked by effect size (|log2 fold|), each row carrying the readout's fold, call, and the full honesty set (railed, confidence, measured ranges, fragility, switch state where present). A candidate that does not move the readout is reported in the summary count, not listed. Nothing in any answer is computed by a shortcut: every number is what the forward product would have served for that perturbation.

**Why this priority**: it is the feature; the inverse question is the second half of the product's reason to exist.

**Independent Test**: on a fixture network where A → B → C and D ⊣ C, an inverse query for "what lowers C" returns the D over-expression and the A/B knockouts with their forward-solve folds exactly equal to individually-run forward solves; "what raises C" returns A/B over-expression and D knockout; an entity with no pinnable root is absent from candidates.

**Acceptance Scenarios**:

1. **Given** an inverse query, **Then** each returned row equals the ordinary forward solve + readout valuation for that candidate, bit for bit (the test runs both paths).
2. **Given** a readout the network cannot resolve, **Then** the endpoint refuses with the readouts endpoint's own honest error, never an empty success.
3. **Given** a direction filter, **Then** only candidates whose CALL matches are listed, and the summary reports how many candidates were solved, how many moved the readout at all, and how many were cut by a cap.

---

### User Story 2 - Nothing else moves (Priority: P1)

The feature is a new endpoint; no existing endpoint, field, or figure changes. The specs/141/142 dumps replay unchanged. The benchmark protocol is untouched.

---

### User Story 3 - The cost is bounded and stated (Priority: P2)

An inverse query is N forward solves. The endpoint takes a candidate cap (registered default resolved in research from measured solve times; the server refuses a cap above its registered maximum) and an optional candidate pre-ranking by the canonical solve's influence scores so the cap keeps the most promising candidates (resolved in research: whether influence at the readout's nodes predicts which candidates move it — measured on a handful of pathways, recorded either way; if it does not, the cap is order-arbitrary and SAYS so in the response). The response always reports solved/total candidates so a capped answer can never masquerade as exhaustive.

---

### Edge Cases

- The chosen readout's nodes include a candidate's own pins (asking "what lowers X" and perturbing X): the candidate is excluded and counted (its answer is the pin, not a propagation).
- Two candidates resolving to identical pin sets (paralogs through one root complex): deduplicated by pin set, both names reported on the row.
- A candidate whose solve does not converge: reported with the standard non-convergence honesty, never silently dropped (FR-006 of the record's conventions).
- An inverse query on a readout the perturbations cannot reach at all: an honest "0 of N candidates move this readout here" — the specs/146 untouched vocabulary, not an error.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001 (forward-path fidelity)**: every served number comes from the ordinary solve + readout code paths; the parity test runs a returned row's candidate through /api/solve + /api/readouts and requires bit equality.
- **FR-002 (the candidate set is the pin protocol's)**: candidates are exactly the distinct pin resolutions of the network's curated entities under the recorded protocol; no other perturbation shapes.
- **FR-003 (additive)**: new endpoint only; the specs/141/142 dumps replay unchanged.
- **FR-004 (honest bounds)**: cap and ordering stated in the response; solved/moved/capped counts always present; refusals reuse the readouts endpoint's honest errors.
- **FR-005 (tests with fixtures that compute)**: the US1 fixture (inverse answers equal forward solves), the direction filter, the dedup and self-exclusion rules, the cap reporting; registered in the asserting list and the CLAUDE.md table.
- **FR-006 (docs)**: docs/API.md gains the endpoint with the cost statement.

## Success Criteria *(mandatory)*

- **SC-001**: the fixture's inverse rows are bit-equal to individually-run forward solves.
- **SC-002**: a capped query can be recognized as capped from the response alone.
- **SC-003**: every honesty field the forward product serves appears unchanged on inverse rows.
- **SC-004**: no existing test, dump replay, or served figure changes.

## Assumptions

- Measured per-solve latency on the served catalog bounds the default cap (research measures it; large pathways may need the cap well below the candidate count).
- WebsiteAngular consumption (a "what moves this?" panel action) is its own spec once served.

## Non-goals

- No new validation axis and no accuracy claim beyond the forward record's.
- No combinatorial perturbations (that is backlog #256), no dose sweeps, no optimization/search beyond ranking the single-pin candidates.
- No async/job machinery in v1: the cap keeps the request synchronous.
- No UI work in this repo.
