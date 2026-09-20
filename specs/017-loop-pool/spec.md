# Feature Specification: Loop as a conserved pool

**Feature Branch**: `feat/017-loop-pool`

**Created**: 2026-09-19

**Status**: Draft

**Input**: User description: Adam's proposal — "if you have a loop you could treat it like a node somehow. where say there are four things going in A,B,C,D with values 1,1,3,2 then the cycle should have the value of 1*1*3*2=6 around. unfortunately you lose the granularity of different parts of the loop having different values and I don't know how you handle negative interactions." Treat a strongly-connected component as one pool whose level is set by what flows into it, instead of iterating around it.

## Why

Every internal edge of a positive cycle has gain exactly 1 under
fold-multiplication AND, so a cycle's baseline is a knife-edge: any leak rails
it to 100x or collapses it to the all-zero root, and which happens depends on
sweep order and node labels (specs/013, specs/014). On 2026-09-19 that coin
decided the sign of three otherwise-correct structural fixes (specs/016:
sharing flipped TP53 by −92; composition edges flipped DSB and TP53 in both
directions). The loop taxonomy (2026-06-12) showed the giant components —
TP53 1,250 nodes, WNT 622, DSB — are catalytic-recycling artifacts dominated
by one or two Reactome reactions, where iterating around the cycle is
meaningless: the recycled species is one pool whose level is set by external
supply. `DS_SCC_BREAK_CATALYST` (freeze the recycled catalyst at its entry
value) was a half-step: TP53 +104, held-out −65, because it reads
signal-carrying catalysts at baseline.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - A loop reads its external supply, not its own echo (Priority: P1)

A pathway modeller perturbs a gene upstream of a recycling loop. Today the
loop's members land at 100x, at 0, or anywhere between depending on node
labels. With the pool rule every member reads the combined external drive:
one value, deterministic, label-independent, and the loop can no longer
amplify or swallow a signal on its own.

**Why this priority**: it removes the mechanism that has set the sign of every
structural fix this month; nothing downstream of a loop is measurable until it
is gone.

**Independent Test**: a positive two-node loop driven from outside at 2x reads
2x at both members (not 100x, not 0); driven by a knockout (0) it reads 0;
relabelling every node of an isomorphic network gives identical activities.

**Acceptance Scenarios**:

1. **Given** a component whose only external input is at fold 2, **When** solved under the pool rule, **Then** every non-pinned member is at fold 2 (±1e-9) and the solve reports the component as pooled.
2. **Given** the same component with the external input at fold 0, **When** solved, **Then** every non-pinned member is at 0.
3. **Given** a component with entries at folds 1, 1, 3, 2 that are co-required (AND) by distinct member reactions, **When** solved, **Then** the pool reads fold 6.
4. **Given** two entries that are alternative producers (OR) of one member, **When** solved, **Then** the pool reads their OR combination, not their product.
5. **Given** an isomorphic relabelling of any of the above, **When** solved, **Then** activities are bit-identical under the relabelling.

---

### User Story 2 - The default is untouched and a fallback is never silent (Priority: P1)

A benchmark operator must be able to run the current solver and the pooled
solver side by side and know, from the output, which components were pooled
and which fell back to the damped fixed point.

**Independent Test**: with the mode unset or at its default, every existing
assertion test file passes and a differential solve matches the previous code
bit for bit; with the mode set, the solve result carries counts of pooled and
iterated components; a misspelt mode is a startup error.

**Acceptance Scenarios**:

1. **Given** the mode variable unset, **When** any network is solved, **Then** the activities equal the pre-feature solver's exactly.
2. **Given** the mode set to a pooling variant, **When** a network with three cyclic components is solved of which one falls back, **Then** the result reports pooled = 2, iterated = 1.
3. **Given** a misspelt mode value, **When** the solver starts, **Then** it raises a configuration error naming the allowed values.

---

### User Story 3 - Negative edges inside the loop are measured, not assumed (Priority: P2)

Adam is not sure how internal negative interactions should be treated. Two
rules are shipped and measured against each other:

- **`pool`** (conservative): a component containing any internal inhibitor or
  depletion edge is not pooled; it keeps the existing damped fixed point.
- **`pool_parity`**: each member reads each entry's fold to the power +1 or −1
  according to the parity of negative edges on its internal path from that
  entry, so a member downstream of an odd number of internal negatives moves
  opposite to the entry. A component whose internal signs are inconsistent (an
  odd negative cycle — genuine negative feedback, e.g. p27 ↔ Cdk2 in Mitotic
  G1) is not pooled under either rule.

**Independent Test**: the traced TP53 case (specs/016): X → M enters a
component in which M depletes T; with X at 0.5, `pool_parity` reads T at 2x and
M at 0.5x, while `pool` falls back and reports the component as iterated.

**Acceptance Scenarios**:

1. **Given** the X → M, M ⊣ T fixture with X = 0.5, **When** solved under `pool_parity`, **Then** M reads 0.5x and T reads 2x, and the component is reported pooled.
2. **Given** the same fixture, **When** solved under `pool`, **Then** the component is reported iterated and activities equal the fixed-point solver's.
3. **Given** a component with an odd negative cycle, **When** solved under either rule, **Then** it is reported iterated.

---

### User Story 4 - Decided on the wide curator set, both axes, prediction first (Priority: P2)

Both variants are measured on the production catalog against the current
production dump, on the held-out split with concentration columns, on the
experimental axis, with McNemar p and per-pathway net, and with the TP53
AKT1/AKT2-knockout cases scored directly. The predictions are committed to git
before the arms run.

**Acceptance Scenarios**:

1. **Given** the pre-registration commit, **When** the arms finish, **Then** each prediction is scored held / failed in the same document, and the verdict follows the pre-stated decision rule.

### Edge Cases

- A component in which every member is pinned by an observation: nothing to
  pool; members keep their pinned values.
- A pinned member inside a component: it is an entry (its fold enters the pool)
  and it is never overwritten.
- A component with no external input at all (a closed loop): pool fold is 1;
  every member reads baseline.
- Entry folds that multiply past the cap: the pool is capped exactly as
  hill_sat caps a reaction (100x), and 0 × anything is 0.
- Members whose baseline differs from the common 0.01: each member reads *its
  own* baseline × pool fold.
- A member reachable from an entry by two internal paths of different parity
  under `pool_parity`: the component's signs are inconsistent and it falls
  back.
- Downstream acyclic nodes and later components read the pooled values in the
  existing topological order; nothing outside the component changes rule.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST offer two new component-resolution modes, `pool`
  and `pool_parity`, selected the same way as the existing modes; the existing
  default MUST remain the default and be bit-identical in output.
- **FR-002**: For a pooled component, each member reaction's *external fold*
  MUST be computed with the normal AND/OR/inhibition rules, with every
  in-component input held at that input's baseline and every pinned member at
  its pinned value.
- **FR-003**: The pool fold MUST be the product of the external folds of the
  entry reactions (those whose external fold ≠ 1), with 0 absorbing and the
  product capped at the same ceiling the AND rule uses.
- **FR-004**: Every non-pinned member MUST read its own baseline × the pool
  fold (`pool`) or × the parity-signed pool fold (`pool_parity`); no iteration
  is performed inside a pooled component.
- **FR-005**: Under `pool`, a component with any internal negative edge
  (inhibitor or depletion) MUST fall back to the existing damped fixed point.
- **FR-006**: Under `pool_parity`, each member's exponent for each entry MUST be
  (−1)^(number of negative internal edges on the path from that entry); a
  component whose sign assignment is inconsistent MUST fall back under both
  modes.
- **FR-007**: The solve result MUST report the number of components pooled and
  the number iterated, so a benchmark cannot silently measure the old solver.
- **FR-008**: Output MUST be invariant to node labels and to edge order for
  pooled components (asserted with the existing relabelling harness).
- **FR-009**: A misspelt mode MUST be rejected at configuration time through
  the existing allow-list mechanism.
- **FR-010**: Both variants MUST be measured on the production catalog against
  the current production dump: held-out and tuning splits with concentration
  columns, McNemar p, per-pathway net, false-change counts, the experimental
  axis, and the TP53 AKT1/AKT2-KO cases scored directly; the pre-registered
  predictions MUST be committed before the arms run.

### Key Entities

- **Component**: a maximal set of mutually reachable nodes; trivial (one node,
  no self-loop) or cyclic.
- **Entry reaction**: a member reaction whose external fold differs from 1 —
  the place a signal enters the loop.
- **External fold**: a member reaction's output ÷ its target's baseline when
  all in-component inputs sit at baseline.
- **Pool fold**: the product of entry folds, capped; the one value the loop
  carries.
- **Parity**: for `pool_parity`, the sign a member reads an entry with,
  determined by the count of negative internal edges between them.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: On the synthetic fixtures, a 2x external drive yields 2x at every
  member (no rail, no collapse) and a knockout yields 0, under both modes.
- **SC-002**: Relabelling an isomorphic catalog moves 0 predictions under a
  pooling mode (the current solver moves 14–22, specs/013).
- **SC-003**: The AKT1-KO → TIGAR case reads UP under `pool_parity`; the 102
  AKT1/AKT2-KO TP53 cases score ≥ 60 correct (currently 2).
- **SC-004**: Held-out net vs the production dump is ≥ 0 for at least one
  variant with p < 0.05, distributed over ≥ 5 pathways and ≥ 10 genes, and
  false change does not rise; the experimental axis is not negative
  (conditioned). Otherwise the feature is a recorded negative result.
- **SC-005**: Every cyclic pathway in the catalog reports its pooled / iterated
  component counts in the arm log.

## Assumptions

- The existing SCC machinery (Tarjan components, topological order, pinned-set
  handling) is reused; the pool replaces only the per-component iteration.
- "Negative internal edge" means an inhibitor or depletion edge whose source and
  target are both in the component; activator edges are positive regardless of
  edge type.
- The product rule follows Adam's statement literally across *entries*; within
  a reaction the existing AND/OR rules apply, so alternative producers do not
  multiply.
- Genuine negative-feedback components are rare (loop taxonomy: Mitotic G1's
  92-node component is the named example); falling back to iteration for them
  is acceptable for this measurement.
- Non-goals: no generator change; no edge deletion; no tuning of a pool
  exponent against the evaluation set; no claim that pooling is biology for
  genuine negative-feedback loops.
