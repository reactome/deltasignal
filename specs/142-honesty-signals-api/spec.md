# Feature Specification: Honesty Signals Through the API

**Feature Branch**: `specs/142-honesty-signals-api`

**Created**: 2026-11-06

**Status**: Draft

**Input**: Usability programme gap 3 (the programme is listed in specs/140-pin-protocol-api/spec.md; Adam 2026-11-06: "continue trying to figure out how to make our system usable", and approved doing all remaining gaps this session). Surface the solver's honesty signals through the API, so the product can tell a user when a number is not a measured magnitude.

## The problem

The solver knows two things about its own answers that the product hides.

1. **Rails.** A readout at the model's stop — fold at or above 0.99 × ceiling, or at or below 0.01 — is a switch pinned against a bound, not a measured level. specs/108 defines rail-dependence and found 8.0% of correct held-out curator cases and 9.5% of correct experimental cases are rail-dependent; specs/054 traced the RAF switch; specs/082 the PTEN/MKRN1 loop. Every rails% figure in the spec record (specs/107 behaviour checks, specs/124 comparisons) is computed by one function, `bench/analysis/score_124.py railed(fold, ceiling)`: `fold >= 0.99 * ceiling or fold <= 0.01`. The product currently shows such a fold as if it were a magnitude.
2. **Non-convergence.** `/api/solve` already reports `converged`, `residual` and the `nonconv_*` diagnostics (specs/066/073/076), but the consuming UI ignores them, so a user cannot tell a converged steady state from a cap-decided answer.

This spec serves (1) from the engine; (2) needs no engine change. The WebsiteAngular consumption of both is its own spec once this is served.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A railed readout says so (Priority: P1)

A user perturbs a gene and reads an entity's predicted level. When that level is at the model's stop, the response says so, flagged by exactly the statistic every recorded rails% figure used — so the product can present "the model rails this readout at its bound" instead of a magnitude the publication record would not back.

**Why this priority**: This is the gap: a 99.2× fold reads as a precise measurement when it is a switch state. Honest surfacing is the difference between a user trusting a number the record supports and one it does not.

**Independent Test**: Request readouts with activities that drive one entity to the ceiling, one to ≤ 0.01×, and one mid-range; the first two are flagged railed, the third is not, and the flags equal the bench function's answers.

**Acceptance Scenarios**:

1. **Given** a readout whose entity-mean fold is ≥ 0.99 × ceiling, **When** `/api/readouts` values it, **Then** its result carries `railed: true`.
2. **Given** a readout whose fold is ≤ 0.01, **When** valued, **Then** `railed: true`.
3. **Given** a readout whose fold is strictly between the stops, **When** valued, **Then** `railed: false`.
4. **Given** a resolution-only request (no activities), **When** resolved, **Then** no `railed` field is served (there is no fold to test).

### User Story 2 - The client never hardcodes the stops (Priority: P2)

A client renders "at the model's stop (≥ 99×)" from values the response carries, so a deployment running a different `DS_BASELINE` (specs/124) shows its own stops, not the default's.

**Why this priority**: The readout endpoint already resolves baseline and cutoffs env-aware (specs/141); the rail stops must follow the same rule or the first non-default deployment lies.

**Independent Test**: Serve with a non-default `DS_BASELINE` and read `rail_low` / `rail_high` from the response's protocol object; they move with the baseline and the flag agrees with them.

**Acceptance Scenarios**:

1. **Given** the default baseline 0.01, **When** any readouts response is served, **Then** `protocol.rail_low = 0.01` and `protocol.rail_high = 99.0`.
2. **Given** `DS_BASELINE=0.001` (ceiling 1000), **When** served, **Then** `rail_high = 990.0` and a fold of 991 is railed while 989 is not.

### User Story 3 - The flag equals the record (Priority: P1)

The served flag reproduces `score_124.railed` on every scored benchmark case of the frozen build, so the product's "railed" and every recorded rails% figure are one statistic.

**Why this priority**: The programme's method (specs/140/141): never a reimplementation judged by eye — a measured parity gate against the bench's own code.

**Independent Test**: The parity gate (bench/analysis/142/) passes on the frozen build's recorded benchmark results.

**Acceptance Scenarios**:

1. **Given** every scored case's recorded fold (`pred_ui`) in the frozen build's curator and experimental results, **When** the bench's `railed(fold, 100.0)` — imported from `bench/analysis/score_124.py`, never reimplemented — and the engine's served flag are computed for the same folds, **Then** zero cases differ; any difference fails the gate.

### Edge Cases

- A fold exactly at a stop (0.01, or 0.99 × ceiling): railed, as the bench's `>=` / `<=` comparisons say.
- Fold 0 (a knocked-out readout reading straight through to zero): railed low.
- A readout resolved but with every group's uuid missing from activities: it reads baseline (specs/141), fold 1.0, not railed.
- `basis: none` (unresolvable): no value, fold or call today; no `railed` either.
- Abundance-mode readouts: `railed` applies to their fold the same way (`status: trivial` / `no_copy` readouts have no fold and so no flag).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Every `/api/readouts` result valued with activities MUST additionally carry `railed`, true exactly when `fold >= 0.99 * ceiling || fold <= 0.01`, with ceiling = 1 / baseline as specs/141 resolves it (env-aware via `DS_BASELINE`).
- **FR-002**: The response's `protocol` object MUST additionally serve `rail_low` and `rail_high` — the resolved stop values (0.01 and 0.99 × ceiling) the flag was computed with.
- **FR-003**: The parity gate (bench/analysis/142/, the specs/140/141 pattern) MUST compare, for every scored case fold recorded in the frozen build 20261026-1349_8b6527e's benchmark results (curator and experimental axes), the bench's `railed` — imported from `bench/analysis/score_124.py`, never reimplemented — against the engine's served flag for the same fold; any difference fails the gate.
- **FR-004**: Everything is additive: no existing response field changes; `/api/solve` is untouched (it already serves `converged`, `residual` and the `nonconv_*` diagnostics gap 3 needs from it).
- **FR-005**: Unit tests with fixtures that compute (never pinning the intermediate the rule reasons about): a readout railed high through a real aggregation reaching the ceiling, one railed low, one mid-range un-railed, the env-aware threshold moving with `DS_BASELINE`, boundary folds exactly at both stops, and the absence of the flag on resolution-only queries.

### Key Entities

- **Rail flag**: per valued readout, the bench's stop test on its fold.
- **Rail stops**: `rail_low` (0.01, scale-free: a fraction of rest) and `rail_high` (0.99 × ceiling, env-aware), served in `protocol`.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The parity gate passes: zero disagreements between the engine's flag and the bench's `railed` over every scored curator and experimental case fold of the frozen build.
- **SC-002**: A consuming client can distinguish a railed readout from a measured one, and name the stop values, using only response fields.
- **SC-003**: Existing consumers see byte-identical behaviour when they ignore the new fields; the full API test suite passes unchanged apart from the new assertions.

## Assumptions

- The stop test on the fold the user's own solve produced is the honest signal to serve. specs/108's full rail-dependence procedure re-solves the case at a mild 0.5×/2× dose; that is a bench measurement procedure, not a per-request service, and is out of scope (non-goal).
- `rail_low` stays 0.01 under any baseline: it is a fraction of rest (1% of resting level), not of the ceiling — matching `railed`'s constant `0.01` in the bench, which the gate enforces at the default and the unit tests document at non-default baselines.
- No RAF-specific pathway flag: the RAF switch is a known-limitation class for gap 4's scope statement; rails cover its mechanism generically (specs/054's railed p-MAPK readouts hit the stop test).
- The gate runs offline on recorded results (no new catalog build, no re-solve): `railed` is a pure function of the fold, so recorded folds are the complete input.

## Non-goals

- No mild-dose re-solve; no rail-dependence classification of calls (specs/108 remains the bench's procedure).
- No change to bench scoring or to any benchmark number.
- No UI work in this repo (WebsiteAngular consumption is its own spec).
- No new catalog build.
- No change to `/api/solve`.
