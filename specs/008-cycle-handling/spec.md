# Feature Specification: Handle cycles correctly in the propagator

**Feature Branch**: `010-loop-handling`

**Created**: 2026-09-16

**Status**: Draft — measurement complete, intervention not yet chosen

**Input**: User description: "we need deltasignal to do the right thing with
these sorts of things where two reactions are connected by inputs and outputs
making it a loop of two. or other cases."

## Context

MP-BioPath removed loops from their networks **by hand**, a curator-reviewed
process we cannot replicate. So DeltaSignal has to handle what they never faced,
and this is now the largest measured deficit with a mechanism attached to it.

On 21,450 scored curator cases at Release97:

| readout sits | cases | accuracy |
|---|---|---|
| acyclic | 16,955 | **0.8641** |
| in a small cycle (2–9) | 115 | 0.7391 |
| in a large SCC (≥10) | **1,922** | **0.6904** |

A 17-point gap over ~9% of the scored set, worth roughly **+1.6pp overall** if
closed to the acyclic rate.

It is the **only case-level discriminator that has survived scrutiny**. Every
other candidate failed:

| candidate | false-change | correct change | verdict |
|---|---|---|---|
| saturated at a 0/100 rail | 50.5% | 51.8% | no |
| shortest path length | median 10 hops | 9 | barely |
| readout in-degree, raw | mean 45.1 | 11.9 | **confounded** |
| readout in-degree, within pathway | median difference **+0.0** | | no |

The in-degree case is the cautionary one: raw it looks like a 4× effect and
would have justified a branch-dilution fix. Controlled per pathway it is higher
in 3 pathways, **lower in 6**, equal in 26. The signal was "false positives
concentrate in dense networks, and dense networks have high-in-degree readouts".

Corroborated from the other direction by the control already on file:
DeltaSignal scored **403/562** on MP-BioPath's own acyclic networks against
their 405, versus **364** on ours — with **0/845 non-converged there against
179/564 here**.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A loop does not collapse to a root that its inputs forbid (Priority: P1)

Someone solves a pathway containing a cycle whose external inputs sit at or
above baseline. Every node in the cycle reports a value consistent with those
inputs, not zero.

**Why this priority**: this is the defect in its own terms. All-zero is a
self-consistent fixed point — every node computes zero from its zero
loop-inputs, the residual is zero, the math is satisfied — but it is the wrong
root when the inputs feeding the loop are ≥ baseline. The iteration can simply
land in the wrong basin. It is a numerical problem with a numerical answer, not
a modelling judgement about mass balance.

**Independent Test**: a synthetic cycle with external inputs pinned at baseline
or above must not report any cycle member at or near zero, from any starting
state.

**Acceptance Scenarios**:

1. **Given** a two-reaction cycle with two external supplies at baseline,
   **When** one supply is knocked out, **Then** the cycle settles at a reduced
   but non-zero level held up by the survivor.
2. **Given** the same cycle, **When** the surviving supply is raised, **Then**
   the cycle scales with it rather than staying pinned.
3. **Given** a cycle whose external inputs are *genuinely* all zero, **When**
   it is solved, **Then** it is still allowed to reach zero — the guard must
   exclude the spurious root without forbidding the legitimate one.

---

### User Story 2 - A stable solve is not reported as non-converged (Priority: P1)

Someone reads the convergence flag and it means what it says.

**Why this priority**: equal-first, because the reported non-convergence rate is
the evidence anyone would use to size the loop problem, and part of it is a
measurement artifact. Reported figures currently overstate the dynamics problem
and understate how much is reporting.

**Independent Test**: a solve whose node values are bit-identical across
iteration budgets must report converged.

**Acceptance Scenarios**:

1. **Given** a minimal two-reaction loop driven 80× above baseline, **When** it
   is solved at iteration caps of 50 through 5000, **Then** the values are
   identical **and** the convergence verdict agrees with that stability.
2. **Given** any solve, **When** the sweep's stopping rule fires, **Then** the
   final verdict is computed on the same quantity the stopping rule used.

---

### User Story 3 - Cyclic readouts stop costing accuracy (Priority: P2)

Someone benchmarking sees cyclic-readout cases score closer to acyclic ones.

**Why this priority**: it is the outcome the work is for, but it is downstream
of the two defects above and cannot be pursued directly without tuning.

**Independent Test**: on the wide curator set, the accuracy gap between
cyclic-readout and acyclic-readout cases narrows, with no loss on acyclic cases.

---

### Edge Cases

- **A cycle whose inputs really are zero.** Zero is then correct, and a guard
  that forbids it would be worse than the defect. US1 scenario 3.
- **A cycle with no external input at all.** Undefined by the inputs; must not
  oscillate or produce non-finite values.
- **Self-loops.** Two exist in the catalog and both are real autocatalysis —
  active caspase-8 in "Caspase-8 processing in the DISC" and PAK1 in
  "Interaction of PAK1 with Rac1-GTP" — not construction artifacts. They must
  keep working.
- **Nested or overlapping cycles.** Cycles span 2–24 reactions, median 4.
- **The giant components.** 251 entity-level SCCs covering ~18,000 nodes come
  from just 30 reaction-level cycles via variant multiplicity — in Class I MHC a
  single reaction carries 488 nodes. Any per-component cost is paid at entity
  scale even though the biology is reaction scale.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: A cycle whose external inputs are at or above baseline MUST NOT
  report cycle members at or near zero.
- **FR-002**: A cycle whose external inputs are genuinely zero MUST still be
  able to reach zero. The guard excludes a spurious root; it does not impose a
  floor on legitimate ones.
- **FR-003**: The convergence verdict and the sweep's stopping rule MUST be
  computed on the same quantity.
- **FR-004**: A solve whose free-node values are bit-identical across iteration
  budgets MUST report converged.
- **FR-005**: Any classification of cycles MUST be order-independent, and its
  stability MUST be asserted — the same component must receive the same
  classification across pathways and across runs.
- **FR-006**: Classification, if used, MUST operate at the reaction level. The
  `precedingEvent` annotation is reaction-to-reaction while the network is
  bipartite entity–reaction, so an entity-level test cannot carry it.
- **FR-007**: No edge may be removed from a logic network, and no cycle may be
  broken structurally in the generator.
- **FR-008**: Every arm MUST be measured on the wide curator set. The
  nine-pathway experimental set MUST NOT be used to accept or reject an arm.
- **FR-009**: The synthetic cycle fixtures MUST live in `test/` with assertions,
  not in a scratch directory.

### Key Entities

- **Cycle (reaction level)**: a strongly connected component of the
  reaction graph, where reaction A precedes B when A outputs an entity B
  consumes. 30 exist across the catalog.
- **Spurious root**: the all-zero assignment, self-consistent but excluded by
  the cycle's external inputs.
- **Convergence verdict**: the reported claim that the solve reached a fixed
  point, which must agree with the stopping rule that ended it.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: No synthetic cycle with inputs at or above baseline reports a
  member at or near zero.
- **SC-002**: A cycle with genuinely zero inputs still reaches zero.
- **SC-003**: The two-reaction fixture driven 80× reports converged, with values
  unchanged from today.
- **SC-004**: On the wide curator set, the cyclic-versus-acyclic accuracy gap
  narrows from its current 17 points.
- **SC-005**: Accuracy on acyclic-readout cases does not fall.
- **SC-006**: Re-running classification on the same catalog yields identical
  assignments.

## Assumptions

- The wide curator set at Release97 on one shared catalog build is the decision
  basis, paired and conditioned on the experiment being unchanged.
- "Near zero" means at or below the baseline of 0.01 on the internal scale.
- Cycles are a solver concern; the networks stay as curators recorded them.

## Non-Goals

- Deleting edges from the logic networks. ~98.5% of edges trace to curator
  assertions — diagram bridges come from curator-drawn diagrams, assembly and
  dissociation from curated complex composition — and only the 4,930 depletion
  edges (1.5% of 334,844) are our own modelling inference. The cofactor feature
  established the pattern: the generator stays faithful, the solver decides how
  to process.
- Breaking cycles structurally in the generator.
- Tuning the anti-collapse magnitude against the evaluation set.
- Re-running the entity-level classification. Its failure is recorded below so
  nobody repeats it.

## Prior Findings This Feature Must Not Repeat

### Entity-level classification does not work

Three successive formulations, all against the same catalog:

| formulation | SCCs judged "real" | why it is wrong |
|---|---|---|
| any internal edge annotated | 49% | a large SCC always contains some annotated link |
| closing/back-edge annotated | 18% | **DFS-order dependent** — gave the same 1,127-node component two different types in two pathways |
| survives on annotated edges only | 0.4% | most edges connect a node to its own reaction and can never be annotated |

At the reaction level the test is well-posed and order-independent: delete every
link not marked "Has Preceding Event" and ask whether the component is still
strongly connected. That yields Type I 2, Type II 16, Type III 4, Type IV 8 —
60% artifact, 40% real.

### Classification by provenance is probably the wrong axis

Since ~98.5% of edges are curator-derived, there is very little to sort on. The
alternative — make the cyclic dynamics correct regardless of provenance — is
simpler, carries no risk of misclassifying real feedback as artifact, and needs
no new generator artifact.

### The anti-collapse machinery already exists and is off

`DS_INHIBITOR_FLOOR_SCOPE` already defaults to `loops`; `DS_INHIBITOR_FLOOR`
defaults to `0.0`, which makes it inert. It was tested once, judged not to help,
and left off — on the nine-pathway experimental set, which is precisely where
cyclic readouts barely exist. That is the same trap that has now reversed five
separate results in this project, and is why FR-008 exists.

### The convergence defect, stated exactly

On a minimal two-reaction loop driven 80× above baseline the solver reports
`converged=false` with residual `4.045534594765421e-6` against a `1e-6`
tolerance. That residual is **identical** at tolerances 1e-6, 1e-8, 1e-10 and
1e-12, and at iteration caps of 50, 100, 200, 400, 1000 and 5000, with node
values bit-identical throughout. The sweep has stopped moving and the final
global check disagrees with it by 4×.

Ruled out while diagnosing: multiple reactions per target (each node has exactly
one), `DS_SCC_BREAK_CATALYST` supply freezing (defaults off), iteration budget,
and damping (λ=0.1 is strictly worse, matching the existing code comment).

**The 4e-6 is not yet explained.** This spec does not claim otherwise, and the
first task of the feature is to explain it.
