# Feature Specification: Derived-edge recycling closures do not carry loop signal

**Feature Branch**: `feat/018-derived-edge-loops`

**Created**: 2026-09-19

**Status**: Draft

**Input**: steer — understand how loops are broken without necessarily breaking them ourselves, and make rules that improve accuracy without foreclosing a later, biologically faithful representation. Findings in [research.md](research.md) §1–3.

## Why

Our giant strongly-connected components are not Reactome's loops. Removing
our own derived edge classes — assembly (member → complex) and depletion
(complex ⊣ free subunit) — takes TP53's component from 836 nodes to 34 and
DSB's from 1,127 to 126; only 2,077 assembly and 847 depletion edges lie
inside cycles, but they weld separate reaction-level cycles into one. MP-
BioPath's hand-curated networks keep the reaction backbone inside our cycles
(input 48%, catalyst 64%, output 49% present) and essentially never contain
our derived classes (dissociation 0%, depletion 0%, assembly 21%); they have
26 cyclic components where we have 258. Every rule applied uniformly to an
SCC this month failed because the SCC was several small loops welded by our
edges, not one loop.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A welded component falls apart into its reaction-level cycles (Priority: P1)

A modeller perturbs a gene in TP53. Today the solve iterates one 836-node
component and lands on a label-dependent side of a knife-edge. With the rule,
the assembly and depletion edges whose source and target share a component
are recognised as recycling closures: they still carry their feed-forward
meaning (a complex cannot exceed its scarcest member; an abundant complex
drains its subunit) but do not feed a value back around the cycle. Component
detection is re-run without them, and the solver iterates the small
reaction-level cycles that remain.

**Why this priority**: it is the one loop intervention grounded in what the
curators actually did, and it removes the mechanism behind this month's
coin-flips without touching the network.

**Independent Test**: a fixture in which one assembly edge welds two
two-node reaction cycles into one component splits into two components under
the rule; the assembly edge still limits the complex feed-forward.

**Acceptance Scenarios**:

1. **Given** two reaction cycles joined only by an assembly edge from a member released by cycle 1 to a complex in cycle 2, **When** solved with the rule, **Then** the solve reports two cyclic components (one before), and a knockout of the member outside both cycles still lowers the complex.
2. **Given** the same fixture, **When** the rule is off, **Then** the solve is bit-identical to the current solver.
3. **Given** a depletion edge from a complex to a subunit inside the same component, **When** solved with the rule, **Then** the depletion still suppresses the subunit when the complex is raised from *outside* the component, and no longer feeds back inside it.
4. **Given** an isomorphic relabelling, **When** solved with the rule, **Then** activities are bit-identical.

---

### User Story 2 - Each derived role is measured on its own (Priority: P1)

Freezing catalyst closures alone was measured in June (TP53 +104, held-out
−65). The roles must be selectable independently so the A/B can attribute:
`assembly,depletion` (our constructs), `catalyst` (Reactome's recycling
closures, re-measured with component recomputation), and all three.

**Acceptance Scenarios**:

1. **Given** the role list `assembly,depletion`, **When** a component's only closures are catalyst edges, **Then** it is left intact and iterated.
2. **Given** a misspelt role, **When** the solver starts, **Then** it raises a configuration error naming the allowed roles.

---

### User Story 3 - The default is untouched and the effect is reported (Priority: P1)

**Acceptance Scenarios**:

1. **Given** the role list unset, **When** any network is solved, **Then** activities equal the current solver's exactly, and the existing `DS_SCC_BREAK_CATALYST` knob behaves exactly as before.
2. **Given** a role list, **When** a network is solved, **Then** the solve reports how many edges were treated as closures, per role, and the cyclic component count before and after.

---

### User Story 4 - Decided on the wide curator set, both axes, prediction first (Priority: P2)

Arms on the production catalog against the current-tree `fixed_point`
control; held-out split with concentration; experimental axis; McNemar p;
per-pathway net; false change; the TP53 AKT1/AKT2-KO cases; relabel churn.
Predictions committed before the arms.

### Edge Cases

- A closure edge whose source is pinned: the pin is read (pins are never
  overwritten), which is feed-forward anyway.
- A closure edge whose ends fall into different components after
  recomputation: it is now an ordinary feed-forward edge and reads the
  upstream component's final value.
- A closure edge whose ends are still in one component after recomputation
  (another cycle through them survives): it reads the component-entry value
  — no feedback — while the surviving cycle is iterated as before.
- Components made only of derived closures (no reaction cycle): become
  acyclic and are evaluated once, exactly.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: A new setting MUST accept a list of roles from {`catalyst`,
  `assembly`, `depletion`}; unset MUST be bit-identical to the current
  solver; the existing `DS_SCC_BREAK_CATALYST` MUST be unchanged.
- **FR-002**: An edge of a listed role whose source and target lie in the
  same strongly connected component (first pass) MUST be marked a closure.
- **FR-003**: Component detection and topological order MUST be recomputed
  with the closure edges excluded (second pass); the solver MUST iterate the
  second-pass components.
- **FR-004**: A closure edge MUST be read at its component-entry value
  (upstream final value when its source is upstream; the initial value when
  it is not) — the existing `supply` mechanism — for activator roles and for
  depletion alike.
- **FR-005**: The solve result MUST report closures per role and cyclic
  component counts before and after recomputation; the API `scc` object
  carries them additively.
- **FR-006**: Output MUST be invariant to node labels and edge order.
- **FR-007**: A misspelt role MUST be rejected at configuration time.
- **FR-008**: The measurement of US4 MUST be pre-registered and committed
  before any arm runs, with the three role lists as separate arms.

## Success Criteria *(mandatory)*

- **SC-001**: On the fixtures, welded components split as specified and
  feed-forward limiting is preserved.
- **SC-002**: On TP53 the largest component after recomputation under
  `assembly,depletion` is ≤ 60 nodes (836 today) and the AKT1-KO solve
  converges.
- **SC-003**: Relabelling the catalog moves 0 predictions under the rule
  (14 today).
- **SC-004**: Held-out net vs the current-tree control ≥ 0 with p < 0.05 for
  `assembly,depletion`, distributed over ≥ 5 pathways and ≥ 10 genes, false
  change not rising, experimental axis not negative; the 102 TP53 AKT cases
  lose ≤ 10 of their 100. Otherwise a recorded negative result.

## Assumptions

- Assembly, depletion and catalyst are the only derived-or-recycling roles
  that close cycles in numbers that matter (census: assembly 2,077,
  depletion 847, catalyst 10,238 inside cycles; regulator 231 is left alone).
- Dissociation edges never close a cycle (their targets are sinks) and are
  not listed.
- Non-goals: no generator change; no edge deletion; no claim about genuine
  negative feedback, which stays on the backbone and is iterated as before.
