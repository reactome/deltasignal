# Feature Specification: Solve the objective the model specifies

**Feature Branch**: `003-solver-objective`

**Created**: 2026-09-10

**Status**: Draft

**Input**: The solver runs a fixed-point iteration in place of the specified
minimisation. Implement the objective. Goal is to exceed MP-BioPath.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Every perturbation returns an answer (Priority: P1)

A researcher perturbs a gene in a pathway containing a feedback loop and asks
for the downstream result. Today, for one pathway in ten, the answer is
whatever the solver happened to be computing when it gave up — and running it
longer produces a *different* answer, not a better one.

**Why this priority**: 122 of 564 scored cases are in this state. Until it is
fixed, 22% of every measurement is unfalsifiable, and no other improvement can
be honestly evaluated in that region.

**Acceptance Scenarios**

1. **Given** any pathway including one with cyclic components, **When** a
   perturbation is solved, **Then** a result is returned with a reported
   quality figure, and there is no outcome equivalent to "gave up".
2. **Given** the same network and perturbation, **When** the solve is run with
   a larger computational budget, **Then** the answer does not change
   materially. Today it changes 40 predictions.
3. **Given** a cyclic component with no exact consistent state, **When** it is
   solved, **Then** the best available state is returned together with an
   honest measure of how far from consistent it is — rather than a failure.

### User Story 2 - Measurements are weighted, not absolute (Priority: P2)

An experimentalist supplies a measurement they are only partly confident in.
Today confidence is ignored for the solve: an observation is nailed to its
stated value and the rest of the network must accommodate it.

**Why this priority**: it is the same missing formulation as Story 1 — the
weights exist in the design and in the code's own parameter record, and are
read by nothing. Independent of Story 1 in value, but delivered by the same
change.

**Acceptance Scenarios**

1. **Given** an observation with low confidence that conflicts with strong
   surrounding evidence, **When** the network is solved, **Then** the returned
   value for that node may differ from the stated observation.
2. **Given** an observation with full confidence, **When** the network is
   solved, **Then** that node's value matches the observation.
3. **Given** a node with no observation and no upstream influence, **When** the
   network is solved, **Then** it stays at baseline rather than drifting.

### User Story 3 - Reported configuration is real (Priority: P3)

Anyone reading a result file can trust that the parameters recorded in it
affected the computation.

**Why this priority**: lower user impact, but it is a truthfulness problem —
results and provenance records currently name parameters that had no effect.

**Acceptance Scenarios**

1. **Given** a results file recording the model-consistency and baseline-prior
   weights, **When** either is changed, **Then** the results change.

### Edge Cases

- A cyclic component with several distinct consistent states: the returned one
  must be reproducible for the same inputs, not dependent on incidental
  ordering.
- A component where no consistent state exists: must return the best available
  state and report the shortfall, not fail. This may be biologically
  meaningful — a genuine oscillator has no steady state.
- Acyclic pathways already solve exactly in one pass and must not become
  slower or less accurate.
- Very large components must remain tractable; a correct answer that takes
  hours is not usable in an interactive tool.
- Conflicting observations that cannot all hold simultaneously must produce a
  weighted compromise, not an arbitrary winner.

## Requirements *(mandatory)*

### Functional Requirements

- **FR1** — Every solve MUST return a result. No input may produce an outcome
  meaning "did not converge".
- **FR2** — Every solve MUST report how far the returned state is from
  self-consistency, on a scale comparable across solves.
- **FR3** — Repeating a solve with a larger computational budget MUST NOT
  materially change the answer.
- **FR4** — Observation confidence MUST influence the result: a fully
  confident observation is honoured, a weakly held one can be overridden by
  surrounding evidence.
- **FR5** — The model-consistency and baseline-prior weights MUST affect the
  computation, or MUST NOT appear in results and provenance records.
- **FR6** — Acyclic portions MUST retain their current exact single-pass
  treatment. Only cyclic components change.
- **FR7** — Aggregation MUST stay n-ary. Inputs to one reaction are combined
  once, not split into a chain of pairwise steps, so saturation is applied a
  single time.
- **FR8** — Every change MUST be measured against MP-BioPath and the
  structural baseline on identical scored cases, reporting per-pathway
  attribution.

### Key Entities

- **Node state** — the activity assigned to each network node; what the solve
  determines.
- **Consistency shortfall** — how far the returned state is from satisfying
  the model at every node. Currently a pass/fail verdict; becomes a reported
  quantity.
- **Observation weight** — a measurement's confidence, currently used only to
  decide whether to pin a node at all.
- **Cyclic component** — a group of mutually dependent nodes. The only part of
  the network whose treatment changes.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC1** — DeltaSignal **exceeds** MP-BioPath's 407 of 564 on identical
  scored cases. Matching is the floor, not the goal.
- **SC2** — DeltaSignal exceeds the structural baseline's 393 of 564.
- **SC3** — Zero cases report a failed solve, down from 122.
- **SC4** — Increasing the computational budget forty-fold changes at most a
  handful of predictions, down from 40.
- **SC5** — Accuracy on the currently-stable subset does not fall below its
  present 237 of 442. The unstable cases must be fixed without breaking the
  stable ones.
- **SC6** — Macro-F1 improves, so that "predict no change" is not rewarded.
- **SC7** — A typical pathway solve stays within a few seconds, keeping the
  tool usable interactively.

## Assumptions

- The ten MP-BioPath experimental pathways remain the evaluation set, scored
  against experimental ground truth at the harness's default cutoffs, with
  perturbations at root inputs.
- Comparisons share one catalog build (see 002 research R1).
- MP-BioPath's per-case predictions are fixed published values.
- The existing decomposition into cyclic and acyclic parts is correct; it was
  independently verified against a second algorithm over 400 random graphs.
- The aggregation rules themselves (multiply for AND, divide for negative,
  average for OR) are correct and out of scope here. The assembly clamp is
  handled in feature 002.

## Evidence *(measured, carried into this spec)*

Same code, same networks, only the stopping point differs:

| iteration budget | correct | failed solves | max iterations used |
|---|---|---|---|
| 500 | 334/564 | 122 | 500 |
| 20,000 | **295/564** | **122** | 20,000 |

| subset | budget 500 | budget 20,000 |
|---|---|---|
| solved in both (442) | 237 | **237** — identical |
| failed in both (122) | 97 | **58** |

All 40 changed predictions are in the failed subset. Forty times the budget
removes not one failure, and weaker damping makes failures worse (122 → 152),
so this is neither a budget shortfall nor oscillation — the iteration does not
reach a consistent state at all, and what it reports is an artefact of when it
was stopped.

The same networks solve perfectly when unperturbed (shortfall 5.2e-18), so the
difficulty is introduced by pinning observations, not by the network alone.

MP-BioPath solves the objective this spec describes and scores 407 of 564
against DeltaSignal's 334 on identical cases.

## Out of Scope

- The assembly aggregation rule (feature 002).
- Network coverage: ~21% of cases have no directed path from perturbation to
  readout. Unaffected by this change.
- The evaluation harness itself.
