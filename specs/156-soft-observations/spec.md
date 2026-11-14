# Feature Specification: Soft observations

**Feature Branch**: `specs/156-soft-observations`

**Created**: 2026-11-14

**Status**: Draft

**Input**: Improvement backlog #253 (the 2026-11-13 approved list, item 3; items 1 and 2 shipped end to end 2026-11-14). The first backlog item that changes the solve itself, so it runs as a registered arm with gates, default OFF.

## The problem

An observation's confidence is documented as a weight and behaves as a binary gate: anything above 1e-6 pins exactly as hard as 1.0, anything at or below is discarded (`OBS_CONFIDENCE_TOL`, pinned by `test_observation_pinning.jl`; CLAUDE.md says so in words). specs/003 FR4 registered the design intent — "a fully confident observation is honoured, a weakly held one can be overridden by surrounding evidence" — and its `Σ ωᵢ(xᵢ − yᵢ)²` term was never implemented (specs/003 research, "Not tested"). The consequence is practical: a user with a NOISY measurement (an RPPA panel, a phospho site, an expression readout) cannot condition a solve on it without asserting it as an absolute truth that overrides the network entirely. Conditioning on one's own measurements is the enabler for user datasets, and specs/144/145 measured (outside the solver) that conditioning on a drug's own measured markers carries real signal.

**A trap this spec exists to avoid**: every recorded benchmark pin was sent with confidence 0.9 as a convention (the bench's `[activity, 0.9]`), and the product sends 1.0. Reinterpreting the EXISTING confidence field as a soft weight would silently soften every recorded figure. Soft observations therefore arrive in their OWN request field; the existing `observations` stay hard forever, whatever their confidence says.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A measurement conditions the solve without overriding it (Priority: P1)

A user perturbs a drug's targets (hard pins, today's protocol) and ALSO supplies measured levels for a handful of readouts they assayed, each with a weight. The solve honours the hard pins exactly as today, and each soft-observed node settles at a weighted reconciliation of its measurement and what the network computes for it — a weight of 1 behaves as a hard pin, a weight of 0 as no observation, and in between the network's own evidence can pull the value away from a noisy measurement. Downstream nodes read the reconciled value.

**Why this priority**: it is the feature — the solver-level form of the conditioning specs/144/145 could only do by post-hoc shrinkage.

**Independent Test**: on a fixture chain A → B → C with A pinned, soft-observing B at a value the network disagrees with yields, at weight 1, exactly the hard-pin answer; at weight 0 exactly the unobserved answer; at intermediate weights a value strictly and monotonically between the two, and C reads it.

**Acceptance Scenarios**:

1. **Given** a solve with only hard `observations`, **When** the feature is enabled, **Then** every number is byte-identical to today (soft machinery engages only when soft observations are present).
2. **Given** a soft observation at weight 1, **Then** the node serves the observed value exactly as a hard pin would; at weight 0 (or at/below the existing tolerance) it is discarded as today.
3. **Given** a soft observation at weight w in between, **Then** the node's served value lies between the unobserved solve's value and the observed value, monotonically in w, and the solve still reports convergence honestly.

---

### User Story 2 - The recorded figures cannot move (Priority: P1)

Every recorded benchmark, parity dump and served figure is untouched: the feature is OFF by default (it changes the solve — constitution IV), the existing `observations` field keeps hard semantics under every setting, and with the flag on but no soft observations in the request the solve is bit-identical. The benchmark protocol never uses soft pins.

**Independent Test**: the specs/141/142 dumps replay unchanged with the flag on; the fixture equality test runs the same request with the flag off and on.

---

### User Story 3 - The registered measurement: does solver-level conditioning help? (Priority: P2)

Before anything is served, a registered arm answers the measurable question: conditioning the solve on a drug's own measured readouts, does prediction of a HELD-OUT readout improve? On korkut (FIT), per (drug set, readout): pin the targets hard exactly as the recorded protocol does, soft-observe the drug set's OTHER measured readouts at their measured folds, predict the held-out readout; score against the recorded measured call. Sweep the weight over {0.25, 0.5, 0.75} (registered grid, no fitting). The gate: balanced accuracy on held-out readouts must beat the unconditioned recorded baseline with a 95% bootstrap interval (resampling drug sets) excluding zero, for some registered weight, with the specs/051 discipline (korkut is FIT; nothing outside FIT is touched). The curator and experimental axes are a no-harm guard only (they have no measurements to condition on; they run with the flag off, unchanged by construction).

**Why this priority**: the serving decision rests on it; a solver feature that doesn't beat not-conditioning serves only API surface, and the honest outcome is recording that.

**Acceptance Scenarios**:

1. **Given** the registered arm's run, **Then** per-weight held-out balanced accuracy with bootstrap intervals is recorded in research.md whatever it shows, and the serving decision follows the gate.
2. **Given** a gate failure at every weight, **Then** the feature may still ship as OFF-by-default API surface ONLY if Adam decides the user-dataset enabler is worth serving unmeasured — recorded as his call, never a silent default.

---

### Edge Cases

- A soft observation on a node that is also hard-pinned: the hard pin wins; the soft entry is reported as overridden, never silently averaged.
- A soft observation on a pool-managed or cofactor-held node follows the same precedence the hard pin has today (the existing rules decide; no new special case).
- Several soft observations on one node with different values: an error, as hard observations already treat contradictions.
- Weights outside [0, 1]: an error, as activity ranges already are.
- Convergence: a soft-observed solve reports residuals over the same test as today; the reconciliation must not mask non-convergence.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001 (the solve unchanged by default and by absence)**: flag default OFF; with the flag on and no soft observations, bit-identical solves; the specs/141/142 dump replays unchanged with the flag on. The existing `observations` field keeps hard semantics under every setting.
- **FR-002 (the reconciliation rule)**: a soft observation's node settles, at the solve's fixed point, at the weight-blend of its observed value and its network-computed value; weight 1 equals the hard pin, weight at/below the existing tolerance is discarded; monotone in the weight on the fixture. The exact iteration form is resolved in research from the solver's own structure (the damped sweep and the pools) and pinned by tests that compute.
- **FR-003 (own request field)**: soft observations arrive separately from `observations` (API and CLI), each as (node, value, weight); the response reports which nodes were soft-observed, their weights, and any overridden by a hard pin. Malformed entries are 400s, the house pattern.
- **FR-004 (the registered arm)**: US3 exactly as stated — korkut leave-one-readout-out, weight grid {0.25, 0.5, 0.75}, balanced accuracy vs the unconditioned record, 95% drug-set bootstrap excluding zero; role-guarded to FIT; numbers recorded whatever they are.
- **FR-005 (tests with fixtures that compute)**: the US1 fixture (equality at weight 1, absence at 0, monotone between, downstream propagation), the precedence and error cases, and FR-001 equality; registered in the asserting list and the documentation table.
- **FR-006 (honest reporting)**: the solve's diagnostics name the rule, the count of soft-observed nodes, and the overridden ones; a non-converged soft solve is reported exactly as a non-converged hard solve is.

## Success Criteria *(mandatory)*

- **SC-001**: with the feature off, and with it on but unused, every existing test, dump replay and served figure is unchanged.
- **SC-002**: the US1 fixture's three-point behaviour (hard at 1, absent at 0, monotone between) holds on every registered weight.
- **SC-003**: the US3 measurement is on the record with intervals, and the serving decision cites it.
- **SC-004**: a user request carrying hard pins plus soft measurements returns a converged (or honestly non-converged) solve whose soft-observed nodes and weights are named in the response.

## Assumptions

- korkut's measured folds (the specs/155 pairs machinery) supply the soft-observation values for the arm without re-running any recorded benchmark.
- The arm runs on the served catalog with the recorded protocols; nothing is fitted (the weight grid is fixed up front).
- WebsiteAngular consumption (a UI for entering measured levels) is its own spec, only if the gate passes and the feature serves.

## Non-goals

- No reinterpretation of the existing `observations` confidence — ever.
- No fitted weights, no per-node learned confidences (specs/128–131 closed fitting at current data); the grid is fixed and registered.
- No change to any benchmark protocol or recorded figure; the bench never sends soft pins.
- No UI work in this repo.
