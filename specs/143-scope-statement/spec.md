# Feature Specification: The Scope Statement, Served

**Feature Branch**: `specs/143-scope-statement`

**Created**: 2026-11-06

**Status**: Draft

**Input**: Usability programme gap 4 (the programme is listed in specs/140-pin-protocol-api/spec.md; Adam 2026-11-06: "continue trying to figure out how to make our system usable", and approved doing all remaining gaps this session). A user of the product has no way to know what the system is validated for and what it is not. Serve a versioned, provenance-carrying scope statement from the engine, so the product shows the record's answer — not marketing copy that drifts from it.

## The problem

The record is precise about what this system can and cannot do: curator-style directional reasoning inside curated pathway scope is validated on three axes with exact numbers and provenance (docs/RESULTS.md, specs/111/112); drug-response phosphoproteomics, magnitude calibration, switch pathways and proteins the pathway only synthesizes are known limitation classes (specs/138, 108/124, 054, 135). None of this reaches a user. A product that shows a fold to two decimal places but never says "this class of question is outside what was validated" misleads by omission — the same disease gaps 1–3 treated, at the level of the whole tool.

The statement must live beside the record and be served, not hardcoded in a client: numbers are tied to a catalog build and a solver commit, clients multiply, and a statement embedded in UI code drifts the first time the record moves.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A user can read what the numbers mean (Priority: P1)

A user opens the product's scope view and reads: what kind of question the system answers, the measured accuracy on each validated axis with its provenance (build, solver, date, spec), and the named classes of question where the record says the output is not to be trusted.

**Why this priority**: This is gap 4 itself. Gaps 1–3 made the served numbers the measured ones; this one tells the user which questions those measurements cover.

**Independent Test**: Request the scope endpoint and check every validated-axis entry carries a number with complete provenance, and every limitation class names its spec record.

**Acceptance Scenarios**:

1. **Given** a running server, **When** the scope endpoint is requested, **Then** it returns the repo-versioned statement verbatim: a summary of what the system does, validated axes each with value + provenance (catalog build, solver commit, date, spec), limitation classes each with a plain-language statement + spec reference, and the protocol pointers (pins, readouts, rails, convergence — the gap 1–3 services).
2. **Given** the statement file is missing or malformed on the serving host, **When** the endpoint is requested, **Then** it refuses with an allow-listed error — it never fabricates or partially serves a statement.

### User Story 2 - The statement says when it predates what is being served (Priority: P2)

The statement's numbers were measured on a specific frozen build. The response compares that build to the catalog the server is actually serving and says whether they match, so a user on a newer catalog knows the figures are from a predecessor.

**Why this priority**: The 2026 incident record includes an API silently serving a July catalog for two months. A scope statement that silently outlives its build would repeat that mistake at the trust level.

**Independent Test**: Serve a bundle whose catalog id differs from the statement's provenance build; the response flags the mismatch.

**Acceptance Scenarios**:

1. **Given** the served catalog equals the statement's provenance build, **When** requested, **Then** the response carries the served catalog id and `measured_on_served_catalog: true`.
2. **Given** they differ (or the served catalog is unknown), **Then** `measured_on_served_catalog: false` with both ids present (served id `null` when unknown).

### Edge Cases

- Statement file present but invalid JSON, or missing a mandatory section: refuse (allow-listed error), never serve a partial statement.
- A numbers entry missing any provenance field: the repo's own test fails — such a statement cannot ship.
- Server with no catalog information: `served_catalog: null`, `measured_on_served_catalog: false`.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The engine MUST serve the scope statement from a single versioned file in the repository (`src/api/scope.json` — amendment 1, 2026-11-06: the containers mount `src/` but not `docs/`, so a statement under docs/ would be unreadable exactly where it is served; the registration first said `docs/scope.json`), verbatim, via a new GET endpoint; the statement is authored and reviewed as part of the record, not generated.
- **FR-002**: Every validated-axis entry in the statement MUST carry: the measured value, the catalog build, the solver commit, the measurement date, and the spec that recorded it. A repo test enforces this on the file itself, so an unprovenanced number cannot ship.
- **FR-003**: Every limitation entry MUST carry a plain-language statement and its spec reference. The shipped statement MUST cover the recorded limitation classes: drug-response phosphoproteomics / missing cell context (specs/062/138), magnitude calibration and rails (specs/108/124 — and the per-readout `railed` flag of specs/142 as the live signal), switch pathways (specs/054), proteins the pathway only synthesizes (specs/134/135), readouts unreachable in the curated pathway (specs/080), and the curator knowledge ceiling (the 0.733 curator-vs-experiment bound).
- **FR-004**: The response MUST carry the served catalog id beside the statement's provenance build and a `measured_on_served_catalog` flag (US2).
- **FR-005**: A missing or malformed statement file MUST refuse with an allow-listed error constant; no fallback text.
- **FR-006**: Everything is additive: no existing endpoint or response changes.
- **FR-007**: CPPA-derived numbers MUST NOT appear in the statement (the data agreement keeps them private); its verdict may be stated in words.
- **FR-008**: Unit tests: the provenance-completeness test on the shipped file (FR-002/003's enforcement), the happy path, the mismatch flag both ways, and the missing/malformed refusals.

### Key Entities

- **Scope statement** (`docs/scope.json`): summary, validated axes (with provenance), limitations (with spec refs), protocol pointers.
- **Catalog match**: served catalog id vs the statement's provenance build.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A client can render, from one request, what the system is validated for (with numbers and provenance) and the named classes where it is not — no client-side copy required.
- **SC-002**: Every number in the served statement is traceable to build + solver + date + spec, enforced by a test that fails the build otherwise.
- **SC-003**: A deployment serving a catalog other than the statement's build is visible as such in the response.
- **SC-004**: Existing consumers are unaffected; the full API suite passes unchanged apart from new assertions.

## Assumptions

- The statement is hand-written against docs/RESULTS.md and the specs, and reviewed like code; the enforced invariant is provenance completeness, not automatic derivation (deriving prose from RESULTS.md would be a generator nobody reviews).
- Headline figures shipped initially: the frozen record (build `20261026-1349_8b6527e`, solver `05dc867`: curator held-out 89.77% / 19,000 cases; experimental 643/849, equal to MP-BioPath's published figure) and HPN-DREAM SC1A (mean AUROC 0.693, 10th of 74, scored once) — the three axes the record quotes, with the curator ceiling stated beside them.
- The WebsiteAngular rendering of the statement is gap 4's UI half, specced in that repo together with gap 3's honesty surfacing.

## Non-goals

- No UI work in this repo.
- No automatic generation of the statement from the record.
- No new measurements; the statement reports the record.
- No CPPA numbers (FR-007).
- No change to `/api/health` (the statement endpoint reads the same catalog info, additively).
