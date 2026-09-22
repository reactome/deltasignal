# Feature Specification: Handle cycles correctly in the propagator

**Feature Branch**: `010-loop-handling`

**Created**: 2026-09-16

**Status**: **CLOSED — NOT PROCEEDING**, 2026-09-16. All three premises failed
re-measurement on the same day the spec was written. The cyclic-readout accuracy
gap is a between-pathway confound, not a loop effect. Kept in full as a negative
result per constitution III. Do not restart this without new evidence that
survives a within-pathway control.

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

**Why this priority**: this is the defect in its own terms — a cycle whose
inputs are at or above baseline should not read zero.

> **REVISED 2026-09-16.** The stated *cause* was wrong. Iterating the same map
> from seven starting states between 0.0 and 1.0 converges to one identical
> root, so there is no second basin to land in and nothing spurious to exclude.
> The outcome this story wants is still right, and OR-joined cycles already
> deliver it; AND-joined cycles do not. See "What the real defect is".

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

**Why this priority**: WITHDRAWN 2026-09-16 — the defect does not reproduce.
At stock config the residual tracks the tolerance, values differ between
tolerances, and every non-converged case clears with more iterations. Kept as a
regression property (FR-003, FR-004), not as work to do.

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
- **Spurious root**: ~~the all-zero assignment, self-consistent but excluded by
  the cycle's external inputs~~ — RETIRED. Measurement found a unique fixed
  point, so this entity does not exist.
- **Recycled co-input**: a cycle member consumed by the reaction that
  regenerates it. Under AND aggregation it becomes a required co-input, so no
  external input can compensate for a knockout elsewhere in the cycle. This is
  the mechanism that replaces the spurious root.
- **Convergence verdict**: the reported claim that the solve reached a fixed
  point, which must agree with the stopping rule that ended it.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: No synthetic cycle with inputs at or above baseline reports a
  member at or near zero.
- **SC-002**: A cycle with genuinely zero inputs still reaches zero.
- **SC-003**: ~~The two-reaction fixture driven 80× reports converged, with
  values unchanged from today.~~ WITHDRAWN — it already reports converged.
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

**Measured at scale, 2026-09-16, and NOT adopted.** `DS_INHIBITOR_FLOOR=0.1` on
the wide curator set, paired and conditioned on one catalog build:

| arm | scored | correct | accuracy | macro-F1 |
|---|---|---|---|---|
| floor 0.0 (default) | 21,450 | 18,083 | 0.8430 | 0.8079 |
| floor 0.1 | 21,450 | 18,100 | 0.8438 | 0.8090 |

93 predictions changed, net **+17** (+29 / −12); per pathway −9
Transcriptional_regulation_by_RUNX2, −1 Signaling_by_Insulin_receptor, +2
Signaling_by_ERBB2. So the nine-pathway rejection does not survive at scale, as
FR-008 anticipated.

**But the mechanism check fails.** If the floor works by excluding the spurious
all-zero root, the gain must land on cyclic readouts. Splitting the 93 changes:

| readout | gained | lost | net |
|---|---|---|---|
| cyclic | 2 | 1 | **+1** |
| acyclic | 19 | 10 | **+9** |

The gain is on acyclic readouts, where anti-collapse should not act. An acyclic
readout can sit downstream of a loop, so this is not proof against the
mechanism — but +1 where the mechanism is supposed to operate does not
demonstrate it either. Per US1, the requirement is to show the spurious root is
excluded, not to bank a +17 of unknown origin.

**This is why FR-001 is stated as a root-selection property and not as an
accuracy target**, and why "tuning the anti-collapse magnitude against the
evaluation set" is a non-goal. One magnitude has been measured and recorded;
sweeping for a larger number would be fitting a constant to the benchmark with
no mechanism behind it. The +17 stays on file as an unexplained general
inhibition effect to be revisited once the dynamics are understood.

### The convergence defect — DOES NOT REPRODUCE (2026-09-16)

Recorded earlier: a minimal two-reaction loop driven 80x reports
`converged=false` with residual `4.045534594765421e-6` against a `1e-6`
tolerance, that residual **identical** across tolerances 1e-6 to 1e-12 and
iteration caps 50 to 5000, with node values bit-identical throughout.

**Re-run at stock config on clean `src`, that is not what happens.** The same
fixture script, same parameters, DS_* env verified inside the container:

| cap | tol 1e-6 | tol 1e-8 | tol 1e-12 |
|---|---|---|---|
| 50 | false, 1.15e-5 | false, 1.15e-5 | false, 1.15e-5 |
| 100 | **true, 7.00e-7** | **true, 7.01e-9** | false, 5.19e-10 |
| 200+ | **true, 7.00e-7** | **true, 7.01e-9** | **true, 7.01e-13** |

The residual **tracks the tolerance** (~0.7x it, as expected from a rule that
stops on the first sweep below tol), node values **differ** between tolerances
(0.6019277… / 0.6019291… / 0.60192915…), and every non-converged cell is an
honest iteration-budget shortfall that clears with more iterations. Per-node
attribution puts the whole residual inside the 4-node SCC, correctly detected
and correctly ordered, with the acyclic nodes exactly consistent.

The recorded numbers also disagree with a verbatim re-run in the values
themselves (A=79.98 recorded, A=60.19 now), so **that run carried an
unrecorded configuration**. The measurement cannot be reproduced and nothing
should be built on it.

**US2 and SC-003 are therefore withdrawn** pending a reproducible case. FR-003
and FR-004 remain as regression properties worth asserting — they are just not
currently violated.

### The spurious root — the framing is WRONG (2026-09-16)

The spec was built on "all-zero is self-consistent and the iteration lands in
the wrong basin". Tested directly on a two-supply cycle with one supply knocked
out and the other raised 50x, iterating the identical map from starting states
0.0, 0.01, 0.1, 0.3, 0.5, 0.8 and 1.0:

    every start converges to the SAME state (A=0.00026, B=0.07756)

There is **one fixed point, not two**. The solver is not selecting the wrong
root — the low value *is* the unique root. A floor that pushes the iterate away
from zero is therefore not "excluding a spurious root"; it is displacing the
only root the model has. That is consistent with the floor's +17 landing on
acyclic readouts: it was never acting on a basin problem, because there is no
basin problem.

### What the real defect is: AND semantics on a recycled input

The OR-joined cycle already satisfies all three US1 acceptance scenarios today
(baseline = 1.0 throughout):

| case | A | B |
|---|---|---|
| both supplies | 1.0 | 1.0 |
| one supply out | 0.334 | 0.668 |
| both out | 0.0022 | 0.0022 |
| one out, other 50x | 13.8 | 31.5 |

Held up by the survivor, scales with it, still reaches zero when the inputs
genuinely go. Nothing to fix.

The **AND-joined** cycle, same topology, is where it breaks:

| case | A | B |
|---|---|---|
| both supplies | 1.0 | 1.0 |
| one supply out | 0.00013 | 0.00278 |
| one out, other **50x** | 0.00026 | 0.078 |

One external input 50x above baseline, and the cycle still sits ~13x *below*
baseline. This is the behaviour described, and it is a property of the AND
semantics, not of root selection: the cycle's own recycled product is treated
as a required co-input, so a knockout anywhere in the cycle cannot be
compensated by any other input.

This matters because cycle-internal positive edges **are** AND in the generated
networks (`and_or` is a function of edge sign), so the artificial-looking
fixture is the common real case, not a corner.

**Not yet established**: that this is what drives the 17-point cyclic gap. The
whole-set error structure does not obviously support it — errors skew to false
UP (1,165) over false DOWN (883) — and the cyclic subset has not been split out
that way. The mechanism is reproduced; its contribution to the gap is not
measured. That measurement is the next task, replacing the 4e-6 as the feature'"'"'s
starting point.


## Why this feature is closed

Six independent checks, all on the wide curator set at Release97 or on
fixtures, all pointing the same way:

1. **The floor does not act on cycles.** `DS_INHIBITOR_FLOOR=0.1` gains +17
   overall but lands on acyclic readouts (net +9) not cyclic ones (net +1).
2. **The convergence defect does not reproduce.** Residual tracks tolerance;
   values differ between tolerances; every failure clears with more iterations.
3. **There is no spurious root.** Seven starting states, one fixed point.
4. **Cycles do not amplify.** Log-log gain through a 2-cycle is 0.995 versus
   0.997 through an acyclic chain of the same length — if anything slightly
   *less* responsive.
5. **OR-joined cycles already behave correctly** on every US1 scenario.
6. **The cyclic penalty is a between-pathway confound.** This is the one that
   closes it.

### The confound, stated exactly

Raw, the gap replicates: cyclic readouts score 0.7192 over 2,678 cases against
0.8580 acyclic over 19,080. On truly-NORMAL cases cyclic readouts call a change
42% of the time versus 8.4% acyclic.

Reach explains most of it. False-change rate on NORMAL cases by the fraction of
the pathway that can reach the readout:

| reach | cyclic n | cyclic FC | acyclic n | acyclic FC |
|---|---|---|---|---|
| 0–20% | 160 | 0.194 | 8,484 | **0.023** |
| 20–40% | 992 | 0.461 | 1,893 | **0.116** |
| 40–60% | 126 | 0.294 | 1,552 | **0.285** |
| 60–80% | 76 | 0.355 | 344 | **0.483** |

The acyclic column alone rises monotonically from 0.023 to 0.483. Over-coupling
is the axis.

Controlled **within pathway** — 11 pathways with at least 20 cyclic and 20
acyclic NORMAL cases — the cyclic effect vanishes:

    median difference -0.013; cyclic worse in 5, better in 6

That is the same signature as readout in-degree, which looked like a 4x effect
raw and gave a within-pathway median of +0.0. Cyclic readouts concentrate in
hard pathways; they are not themselves harder.

Two pathways do show a large positive difference — DNA_Double-Strand_Break_
Repair (+0.346) and Transcriptional_regulation_of_pluripotent_stem_cells
(+0.195). If loop work ever restarts, it starts there, on those two, not on a
catalog-wide mechanism.

### What survives

The **AND-on-a-recycled-input collapse is real as behaviour** — the fixture
proves an external input 50x above baseline cannot hold a cycle up — but it is
not measurably costing accuracy. It stays pinned by `@test_broken` in
`test/test_cycle_handling.jl` so that if it is ever fixed for another reason,
the change is visible.

The tool is `bench/analysis/cyclic_error_direction.py`. It takes an existing
case dump and a catalog; it does not re-solve anything.
