# Feature Specification: Serve bistable components as switches

**Feature Branch**: `154-switch-states-api`

**Created**: 2026-11-13

**Status**: Draft

**Input**: Improvement backlog #251 (the 2026-11-13 improvement list, item 1; Adam approved starting it: "I agree with starting with item 1", and confirmed the intent — "if we know something is either on or off like a receptor then we would treat it that way mathematically" — with the refinement that the equations already have two stable states, so the feature changes what is *served*, never how the solve propagates).

## The problem

Some curated loops are genuinely bistable dynamical systems: positive feedback with enough gain creates two stable resting states and nothing stable between. The solver already lives in that reality — every RAS/RAF perturbation rails p-MAPK at about 50× because 50× is the coordinates of the ON state, not a magnitude (specs/054 traced RAF, flip threshold about ±3.5% on HRAS; specs/082 traced the PTEN/MKRN1 loop; specs/108 defines rail-dependence). Today the product serves "fold 50×, railed: true" — honest, but it still presents a state's coordinates as if they were a level. The scale: 16% of scored benchmark folds rail (specs/142's 3,922 of 24,356), and 8.0% of correct held-out curator / 9.5% of correct experimental answers are rail-dependent (specs/108). The served scope statement's `switch_pathways` limitation describes this class in words; the API should describe it in data, per answer.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A railed readout is served as a switch state, not a magnitude (Priority: P1)

A user perturbs HRAS in the RAF/MAPK cascade. Today every MAPK readout shows "~50×, at the model's bound". With this feature the answer carries, in data: this readout sits on a switch; its two resting states read ~50× (high) and near zero (low); your perturbation selects the high state. The UI (its own spec, later) can then say the biologically correct thing instead of a pseudo-magnitude.

**Why this priority**: it is the whole feature — the largest honest error class after pathway scope, turned from a flagged number into the right kind of answer.

**Independent Test**: solve HRAS over-expression in the RAF pathway on the served catalog; the railed MAPK readouts carry switch fields whose two states match the specs/054 record; a readout outside any bistable component carries `switch: false` and no state fields.

**Acceptance Scenarios**:

1. **Given** the RAF/MAPK network and an HRAS over-expression pin, **When** readouts are requested for p-MAPK entities, **Then** each railed MAPK readout is served `switch: true`, `state: "high"`, and an `other_state_fold` near the low resting state — and its existing `fold`, `call` and `railed` fields are unchanged from today byte for byte.
2. **Given** the same network and a perturbation that selects the low state, **When** readouts are requested, **Then** the same readouts serve `state: "low"` with `other_state_fold` near the high state's value.
3. **Given** a readout in a pathway with no bistable component on its route, **When** readouts are requested, **Then** it serves `switch: false` and no `state`/`other_state_fold` fields.

---

### User Story 2 - Honest absence: no second point claimed where none was found (Priority: P1)

A user solves a pathway whose cyclic component does not converge, or where the search finds no second stable point within budget. The product must not guess: a component with one found point serves `switch: false`; a non-converged component serves no switch fields at all, because its single point is not trusted (specs/066).

**Why this priority**: the feature's value is honesty; a guessed switch is worse than a railed number.

**Independent Test**: a monostable fixture serves `switch: false` everywhere; a non-converged solve serves no switch fields on the affected readouts.

**Acceptance Scenarios**:

1. **Given** a monostable cyclic network, **When** readouts are requested, **Then** every readout serves `switch: false` and the protocol object states the detection rule ran.
2. **Given** a component still above tolerance at the cap (non-converged), **When** readouts are requested, **Then** its readouts serve no switch fields, and the existing non-convergence reporting is unchanged.
3. **Given** a component where the second search hits its iteration budget without converging, **When** readouts are requested, **Then** `switch: false` is served (absence of evidence stated as absence), never a half-converged second state.

---

### User Story 3 - The detected-switch set explains the recorded rails (Priority: P2)

Before anything ships, a registered coverage report measures, on the benchmark record: what share of the 3,922 railed scored folds lie on readouts the detector flags as switches; that every railed RAF MAPK readout under RAS/RAF perturbations is flagged with two states matching the specs/054 record; and that the PTEN/MKRN1 case (specs/082) is flagged. This is a coverage measurement, not an accuracy gate — no call or score may change by construction.

**Why this priority**: it is the evidence that the detector finds the real switches rather than an arbitrary subset; it gates serving, not scoring.

**Independent Test**: the report script runs on the recorded case dumps plus a detection pass over the served catalog, and prints the coverage shares and the two named traced cases' verdicts.

**Acceptance Scenarios**:

1. **Given** the benchmark record and a detection pass, **When** the coverage report runs, **Then** it reports the share of railed folds on flagged readouts, and both named traced cases (RAF, PTEN/MKRN1) are flagged with states matching their specs.
2. **Given** two identical-content catalog builds with different uuids, **When** detection runs on both, **Then** the same readouts are flagged with the same state values to served precision (the specs/013 label-invariance lesson).

---

### Edge Cases

- A component with more than two stable fixed points: serve the two the search found (canonical + the most distant second), and the protocol object says the count is a lower bound.
- Two found points closer than the regeneration noise floor: they are one point; `switch: false`.
- A readout resolving across several components, some bistable and some not: switch fields describe the readout's own value difference between the solve's two states; if the readout's value does not differ beyond the served threshold, `switch: false` even when some upstream component is bistable.
- An observation (pin) inside the bistable component: pins win as always; the second solve honours the same pins, so both states are states *under the user's constraints*.
- Resolution-only queries (no activities): no switch fields, as with `railed` (specs/142).
- The latency budget: the second solve is env-capped; hitting the cap serves honest absence, never a stall.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001 (additive only)**: with the feature on, canonical activities, folds, calls, railed flags and every existing response field are byte-identical to today. A parity check replays the specs/142 rails dump and the specs/141 readout dump against the feature branch: zero differences.
- **FR-002 (detection)**: after the canonical solve, each converged cyclic component is re-solved from at least a high-basin and a low-basin initialization under the same tolerance; a second stable fixed point exists only if the re-solve converges and differs from the canonical point beyond the regeneration noise floor. The canonical solve's own point is always one of the served states; detection never changes it.
- **FR-003 (served fields)**: readout entries valued with activities additionally serve `switch` (bool); when true, `state` ("high" | "low", by the readout's own value at the two points) and `other_state_fold` (the readout's fold at the unselected point, same scale and rounding as `fold`). The response protocol object names the detection rule, the state-difference threshold, and the iteration budget, so clients never hardcode them.
- **FR-004 (honest absence)**: no second point found within budget → `switch: false`; non-converged component → no switch fields on its readouts; resolution-only queries → no switch fields.
- **FR-005 (determinism and label-invariance)**: detection's outcome (flag and state values to served precision) must not depend on uuid order, and must reproduce on a rebuilt identical-content catalog (specs/013; verified as US3 scenario 2).
- **FR-006 (coverage report, registered)**: before serving, a report measures on the benchmark record: the share of the 3,922 railed scored folds on flagged readouts; that every railed RAF MAPK readout under RAS/RAF perturbations is flagged with states matching specs/054 (high ~50×, low near baseline/zero); that the PTEN/MKRN1 PIP3 case is flagged (specs/082). Reported, whatever the shares are; the two named traced cases are the pass/fail bar for serving.
- **FR-007 (unit tests with fixtures that compute)**: a minimal positive-feedback network with known bistability (both points found from the two basins; selection correct under a pin to each side); a monostable control (no second point claimed); noise-floor rejection (two points within the floor are one); API field presence/absence per FR-003/FR-004. No fixture pins the intermediate the rule reasons about.
- **FR-008 (cost bound)**: detection adds at most one bounded re-solve per basin per cyclic component, env-capped iterations; the cap's value is served in the protocol object, and hitting it is FR-004 absence.

### Key Entities

- **Bistable component**: a cyclic component of the solved network with two found stable fixed points under the user's pins, distinct beyond the noise floor.
- **Switch state**: which of the two points the canonical solve selected, named by the readout's own relative value ("high"/"low").
- **Coverage report**: the registered measurement tying the detector to the recorded rails (specs/142's folds, specs/108's rail-dependence, the specs/054 and specs/082 traced cases).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: on the served catalog, the HRAS→RAF scenario serves switch states whose two values match the specs/054 record, while every existing field stays byte-identical (the FR-001 parity replay: zero differences on 24,356 rail folds and 983+272 readouts).
- **SC-002**: the coverage report states what share of the 3,922 railed benchmark folds lie on flagged readouts, and both named traced cases are flagged with matching states — the serving bar.
- **SC-003**: detection reproduces exactly (flags and state values to served precision) across a rebuild of identical content and across uuid relabelling.
- **SC-004**: a full benchmark pathway's solve with detection stays within the served latency budget (the env cap), and hitting the cap is observable as honest absence, never an error or a stall.

## Assumptions

- The canonical solve's fixed point is trusted as-is; detection never replaces or tunes it (Constitution I and IV: default behaviour unchanged; the feature is additive honesty metadata like `railed` and `confidence`).
- The regeneration noise floor (the recorded identical-content rebuild noise) is the right distinctness bar for "two points"; the exact numeric floor is resolved at planning time from the recorded figure.
- "High"/"low" is defined per readout by its own value at the two points, not globally per component — a readout inhibited by the switch is "low" when the component is on.
- The flip threshold ("how hard to flip", e.g. HRAS's ±3.5%) is a follow-up user story once two-state serving stands — out of scope here.
- No change to scoring, benchmarks or any recorded figure; no Boolean semantics anywhere (the continuous propagator is untouched); no UI work in this repo; no new catalog build.
