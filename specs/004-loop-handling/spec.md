# Feature Specification: Handle loops properly, and compare honestly

**Feature Branch**: `004-loop-handling`

**Created**: 2026-09-10

**Status**: Draft

**Input**: the loops published with MP-BioPath were removed manually, while the LNG pathways retain both positive and negative loops. The MP-BioPath optimisation model did not handle loops properly, so DeltaSignal has to. Open question raised with it: are our accuracy problems caused by loops?

## Why this feature exists

Two claims, one measured and one open.

**Measured, and it changes how every previous comparison must be read.** The
MP-BioPath networks that were published with the paper are effectively
acyclic. Across the ten benchmark pathways there is not one self-loop, six of
the ten contain no strongly connected component larger than a single node,
and the largest component anywhere is eleven nodes. The corresponding
DeltaSignal networks contain components of 836, 643, 465 and 443 nodes, with
between 9% and 52% of every pathway's nodes living inside a cycle. MP-BioPath's
407 of 564 was earned on hand-acyclicised networks; DeltaSignal's 365 was
earned on cyclic ones. **These two numbers have never described the same
task**, and every report that has placed them side by side — including the
one written yesterday — overstated the like-for-like deficit by an unknown
amount.

**Open.** Whether the cycles are also the cause of the accuracy gap. The
obvious mechanism is ruled out: the deficit is not non-convergence. Of the 81
cases DeltaSignal gets wrong and MP-BioPath gets right, 63 are in solves that
converged, and on the 179 non-converged cases the two models score exactly
equally (143 each). If loops are hurting, they hurt by making the converged
fixed point a *different and worse* answer, not by preventing one being found.

The governing constraint is constitution principle I. The generator represents
pathways as curators intended; when a faithful representation is hard to
solve, **the fix belongs in DeltaSignal**. Removing edges upstream to make the
graph easy is available as a diagnostic and is not, by itself, an acceptable
outcome.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Know whether it is the propagator or the networks (Priority: P1)

The lead needs to know which half of the system is responsible for the 42-case
deficit before investing in either. Today it is unattributable: DeltaSignal's
propagator and DeltaSignal's networks are only ever measured together.

Running the DeltaSignal propagator over MP-BioPath's *own* published acyclic
networks, scoring the same cases, separates them. If DeltaSignal scores near
407 there, the propagator is sound and the entire gap is network structure.
If it still scores near 365, the networks are exonerated and the propagator
is the problem. Either answer redirects the next three features.

**Why this priority**: it is the cheapest experiment with the largest
decision value, it has never been run, and every other arm in this feature is
easier to interpret once it has.

**Independent Test**: DeltaSignal produces a score on MP-BioPath's networks
for the same scored cases, reported next to both existing numbers.

**Acceptance Scenarios**:

1. **Given** MP-BioPath's published network for a pathway, **When** it is
   loaded and solved by DeltaSignal, **Then** a prediction is produced for
   every case that MP-BioPath itself scored on that pathway, or the case is
   reported as unmappable with a reason.
2. **Given** the control has been run on all ten pathways, **When** results
   are reported, **Then** DeltaSignal-on-MPB-networks, DeltaSignal-on-LNG-
   networks and MP-BioPath's own published predictions appear in one table on
   one case set.
3. **Given** any case is unmappable between the two identifier schemes,
   **When** the comparison is reported, **Then** that case is excluded from
   all three arms rather than counted against one of them.

---

### User Story 2 — State the comparison honestly whatever else happens (Priority: P1)

Whether or not any intervention works, the manuscript claim must stop
comparing a cyclic-network score against an acyclic-network score without
saying so.

**Why this priority**: it is a correctness obligation on work already
published internally, it does not depend on any experiment succeeding, and it
is the one deliverable that cannot be dropped if the interventions all fail.

**Independent Test**: the reported comparison splits cases by whether the
DeltaSignal network they traverse is cyclic, and states the acyclicity
difference in the same breath as the headline number.

**Acceptance Scenarios**:

1. **Given** the scored case set, **When** results are reported, **Then**
   DeltaSignal and MP-BioPath are compared separately on cases whose
   DeltaSignal readout is cycle-resident and on those whose readout is not.
2. **Given** the acyclicity asymmetry, **When** any headline
   DeltaSignal-versus-MP-BioPath figure is stated, **Then** it is accompanied
   by the largest-SCC figures for both network sets.

---

### User Story 3 — Find out whether loops are actually the cause (Priority: P2)

Test the hypothesis by intervention, in increasing order of how much curated
structure each intervention destroys, stopping at the first that both works
and is defensible.

The three cycle sources are already identified and are pathway-specific.
`diagram_bridge` edges are 99.9% cycle-resident and are a connectivity
heuristic rather than curated causality; removing them collapses Cell Cycle
Checkpoints from 777 cycle-resident nodes to 16 and does nothing to TP53.
TP53's 836-node component is catalyst-driven and collapses to 52 without
catalyst edges. WNT's 387 responds to neither.

**Why this priority**: it is the actual question, but it is worth less before
User Story 1 says whether network structure is where the deficit lives.

**Independent Test**: each arm is A/B'd against the current default on the
shared catalog and reported with its per-pathway breakdown, including the
arms that lose.

**Acceptance Scenarios**:

1. **Given** an intervention arm, **When** it is benchmarked, **Then** it is
   reported with macro-F1, per-pathway net change, and the count of changed
   predictions where both arms converged.
2. **Given** an arm improves the score by removing curated causal structure,
   **When** it is evaluated, **Then** it is recorded as a diagnostic result
   and NOT adopted as a default.
3. **Given** an arm reduces cycle-resident nodes but does not improve
   macro-F1, **Then** that negative result is recorded permanently in
   `research.md` with its numbers.

---

### User Story 4 — Establish the ceiling acyclicity could buy (Priority: P3)

A generic feedback-arc-set DAGification is not shippable — it deletes
whatever edges are cheapest to delete, with no biological justification — but
it bounds the prize. If perfect acyclicity is worth four cases, the
hypothesis is dead and no principled loop-breaking is worth building. If it
is worth forty, it justifies real investment.

**Why this priority**: it only informs how much to spend, and only matters if
User Story 3's principled arms are ambiguous.

**Independent Test**: a fully acyclic variant of each LNG network is scored
and reported explicitly as an unshippable upper bound.

**Acceptance Scenarios**:

1. **Given** a DAGified catalog, **When** it is benchmarked, **Then** the
   result is labelled an upper bound and the number of deleted edges is
   reported alongside it.

---

### Edge Cases

- **A pathway has no cycles to break.** Four of the ten LNG pathways have
  comparatively small components. An intervention must be a no-op there, and
  a score change in a pathway an intervention did not touch is evidence of
  measurement noise, not of the intervention working.
- **An intervention disconnects a readout.** Breaking a cycle can remove the
  only path from perturbation to readout, converting a wrong answer into an
  unscoreable one. Cases that become unreachable must be reported separately
  and must not be silently dropped, since dropping them inflates accuracy.
- **The cycle is real biology.** Genuine negative feedback exists in these
  pathways and is part of what the model should capture. An intervention that
  cannot distinguish a curated feedback edge from a recycling artifact is a
  diagnostic, not a candidate default.
- **Identifier mismatch in the control.** MP-BioPath networks use Reactome
  database identifiers; DeltaSignal networks use generated uuids carrying
  stable identifiers. Cases that cannot be mapped must be excluded from every
  arm equally.
- **The control is unexpectedly worse.** If DeltaSignal on MP-BioPath's own
  networks scores below MP-BioPath's published predictions on the same cases,
  that is a propagator finding and must be reported as such rather than
  explained away.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST be able to solve a network supplied in
  MP-BioPath's published four-column form (parent identifier, child
  identifier, polarity, conjunction flag) and produce predictions for the
  same cases the existing benchmark scores.
- **FR-002**: The comparison MUST report DeltaSignal-on-MP-BioPath-networks,
  DeltaSignal-on-DeltaSignal-networks, and MP-BioPath's published predictions
  on one identical case set, with any case unmappable in any arm excluded
  from all arms.
- **FR-003**: The system MUST classify every scored case by whether its
  readout node is cycle-resident in the DeltaSignal network it was solved on,
  and report accuracy and macro-F1 separately for each class.
- **FR-004**: Any reported DeltaSignal-versus-MP-BioPath headline figure MUST
  be accompanied by the largest-strongly-connected-component size for both
  network sets on the same pathways.
- **FR-005**: Each loop intervention MUST be selectable independently and
  MUST default to current behaviour, so that no intervention changes
  predictions until it has been measured.
- **FR-006**: Each intervention arm MUST be evaluated with macro-F1, a
  per-pathway net change, and the number of changed predictions where both
  arms converged.
- **FR-007**: The evaluation MUST report, per arm, how many cases became
  unreachable relative to the baseline, and MUST NOT count a case that became
  unscoreable as anything other than a loss of coverage.
- **FR-008**: Negative and neutral arm results MUST be recorded permanently
  with their numbers, including arms that reduce cyclicity without improving
  accuracy.
- **FR-009**: An intervention that improves the score by deleting curated
  causal structure MUST NOT be adopted as a default; it MUST be recorded as a
  diagnostic bound.
- **FR-010**: Cycle structure measurement MUST be reproducible from the
  shipped catalog by a documented command, reporting per pathway the number
  of components larger than one node, the total cycle-resident nodes, and the
  largest component.

### Key Entities

- **Cycle-resident node**: a node belonging to a strongly connected component
  of more than one node in the directed network as solved.
- **Cycle source**: the classification of an edge by whether removing its
  whole edge type collapses a given component; the currently identified
  sources are diagram bridges, catalyst edges, and depletion edges.
- **Recycling artifact**: a component arising from one or two reaction stable
  identifiers instantiated as many virtual reactions, where a catalyst is
  consumed and regenerated — a representational cycle rather than biological
  feedback.
- **Control arm**: DeltaSignal's propagator evaluated on an externally
  supplied acyclic network, used to attribute error between propagator and
  network.
- **Case mappability**: whether a benchmark case's perturbation and readout
  both resolve to nodes in a given network; determines inclusion in the
  common case set.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The attribution question is answered — DeltaSignal's score on
  MP-BioPath's own networks is reported on the common case set, so the 42-case
  deficit is divided into a propagator share and a network share with numbers.
- **SC-002**: No comparison between DeltaSignal and MP-BioPath is stated
  without the acyclicity difference stated alongside it.
- **SC-003**: Accuracy and macro-F1 are reported separately for cycle-resident
  and acyclic readouts, on both models, so the size of the loop effect is a
  measured quantity rather than a hypothesis.
- **SC-004**: Every loop intervention tried is recorded with its macro-F1,
  per-pathway breakdown, both-arms-converged count and coverage change —
  including the ones that fail.
- **SC-005**: The upper bound on what acyclicity could buy is known, expressed
  as a case count and macro-F1 delta against the current default.
- **SC-006**: Any intervention adopted as a default improves macro-F1 on the
  shared catalog against BOTH ground truths, or is not adopted.
- **SC-007**: The cycle-structure measurement is reproducible from the shipped
  catalog by a documented command.

## Assumptions

- The ten-pathway shared catalog build (LNG main 4ff0408, Reactome Release97)
  is the comparison basis; all arms read the same catalog directory.
- MP-BioPath's published networks in `mp-biopath-pathways/pathways/` are the
  networks its published predictions were computed on. If a spot check shows
  otherwise, User Story 1 is invalid and must be re-scoped.
- The four-column network files encode parent, child, polarity as ±1, and a
  conjunction flag; the exact conjunction convention is to be confirmed
  against MP-BioPath's own reader before the control is trusted.
- the account that loops were removed manually is taken as given; this
  feature measures the consequence rather than re-deriving the provenance.
- Non-convergence is out of scope here. It is being addressed in
  `specs/003-solver-objective`, and 179 of 564 cases currently do not
  converge. Where the two features interact, this one reports the interaction
  and does not attempt to fix it.
- The uuid4 sweep-order artifact means non-converged solves differ run to run;
  the established A/B protocol (shared catalog, per-pathway breakdown,
  both-arms-converged count) applies to every arm.
