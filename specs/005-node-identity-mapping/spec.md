# Feature Specification: Every node maps, both ways, including on the diagram

**Feature Branch**: `005-node-identity-mapping`

**Created**: 2026-09-10

**Status**: Draft

**Input**: User description: "we need to be able to map every node in the database to the LNG node perfectly. and same with the other way around" — and: "two nodes in the reactome pathway diagram that are in the same compartment could be the same thing but in two different places. We need to be able to know which one the uuid was from. this is especially important when interacting with deltasignal through the pathwaydiagram. I know this is confusing. there really should be an id for each node to uniquely identify them but there isn't."

## Why this feature exists

### There IS a per-node id, and we are throwing it away

Adam's premise is that Reactome has no unique id for a diagram node. It does.
Every glyph in a diagram layout carries its own `id`, distinct from the
`reactomeId` naming the entity it draws — measured, R-HSA-1257604 has 114
glyphs with 114 unique ids, and R-HSA-69620 has 154 with 154. That id is
what the pathway browser selects and highlights, so it is exactly the handle
this feature needs.

More pointedly: **the generator already reads these ids and discards them.**
`src/diagram_connectivity.py` collects glyph ids when it pairs producer and
consumer reactions for diagram bridges, then drops them. Nothing records
which glyph a node came from.

So this is not a missing identifier in Reactome. It is a missing column in
our own export, and the fix is to stop discarding what we already handle.

### The duplicate-glyph problem is real and measured

Adam's "same thing in two different places" happens: R-HSA-1257604 draws 4
entities more than once, accounting for 24 of its 114 glyphs (ADP 9, ATP 9,
H2O 3); R-HSA-69620 draws 4 entities across 43 of its 154 (ADP 19, ATP 19,
Ub 3). These are mostly cofactors, which Reactome deliberately draws once
per reaction rather than as one shared hub.

### The correspondence is many-to-many in both directions

This is why a single extra column will not do it. Measured on
R-HSA-1257604: 94 entities have glyphs, 572 generated nodes carry an entity
identifier, and only 56 identifiers appear on both sides. Of that overlap
only 45 have matching counts.

| case | example | glyphs | generated nodes |
|---|---|---|---|
| many glyphs collapse to one node | ADP | 9 | 1 |
| one glyph explodes into many nodes | R-HSA-6811515 | 1 | 32 |
| drawn but not generated | 38 entities | ≥1 | 0 |
| generated but not drawn | 516 identifiers | 0 | ≥1 |

The last row is expected — a diagram covers only its own top-level
reactions, while generation descends further and adds reaction nodes and set
variants — but it has never been *stated*, so nobody can tell an intended
absence from a bug.

### The export promises more than it delivers

`nodes.csv` declares `source_sets`, `chosen_members` and `compartment`.
`compartment` is hardcoded empty, and all three are empty in every row of
every pathway measured; the first two populate only for one node kind.
`diagram_entity_id` holds an entity identifier, not a glyph id, despite its
name. A reader trusting these columns is misled.

### It is already costing 204 benchmark cases

204 of 847 cases are discarded before the solver is consulted because the
readout is an EntitySet. Generation splits sets into their members, so the
set itself has no node, and nothing links the set's identifier to its
members' nodes. All 20 blocked readouts are sets, and they are the canonical
set-shaped readouts of the best-known pathways — phospho-AKT, p-S9/21-GSK3,
phospho-FOXO1/3/4/6, phospho-MAPK dimers — with 116 of the 204 in PIP3 and
33 in RAF.

The members are all there. 15 of the 20 resolve completely by one hop of set
membership. **The remaining 5 are nested sets whose members are themselves
sets**, so resolution must recurse to leaves; with recursion all 20 resolve
and 177 of the 204 cases become scoreable.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — A set-valued readout can be scored (Priority: P1)

Ask for "phospho-AKT" and get an answer, rather than silence. The entity has
no node of its own because it was split into members, so resolving it means
finding those members and combining their values into one figure for the
set.

**Why this priority**: it is 177 recoverable cases, it unblocks the largest
single distortion in the benchmark denominator, and it is the smallest slice
that delivers value on its own.

**Independent Test**: a case whose readout is a set produces a prediction,
and the value is reproducible from the member values by a stated rule.

**Acceptance Scenarios**:

1. **Given** a readout that is a set with no node of its own, **When** it is
   resolved, **Then** the members it was split into are returned, with the
   relationship recorded, and the case is scored.
2. **Given** a set whose members are themselves sets, **When** it is
   resolved, **Then** resolution recurses to leaves and every leaf is
   returned, with its depth.
3. **Given** a member that appears at several positions, **When** the set
   value is formed, **Then** those positions are combined into one member
   value first, so a member split into three positions does not outweigh one
   split into two.
4. **Given** a set that cannot be fully resolved, **Then** it is reported as
   partially resolved with the missing members named — never silently
   averaged over whatever was found.

---

### User Story 2 — Every entity resolves, and the gaps are declared (Priority: P1)

For any entity in a pathway, state which nodes represent it and how — or
state explicitly that it is not represented and why. The same in reverse:
for any node, which entity or entities it stands for.

**Why this priority**: this is Adam's actual requirement, it is what makes
User Story 1 a consequence rather than a special case, and without the
reverse direction a solved result cannot be explained back to a user.

**Independent Test**: a completeness report over a pathway lists zero
unexplained absences in either direction.

**Acceptance Scenarios**:

1. **Given** any entity participating in a reaction, **When** it is looked
   up, **Then** either its nodes are returned with the relationship, or it
   appears in the declared exclusion list with a stated reason.
2. **Given** any generated node, **When** it is looked up, **Then** at least
   one entity is returned.
3. **Given** an entity represented several ways at once, **Then** all are
   returned, each labelled, rather than one being chosen silently.
4. **Given** a deliberately corrupted mapping, **When** the completeness
   check runs, **Then** it FAILS. A completeness check that cannot fail is
   the specific defect this project has shipped before.

---

### User Story 3 — Clicking a glyph selects the right node (Priority: P2)

In the pathway browser, clicking one of the two ATP glyphs must resolve to
what that glyph means, and a solved node must resolve back to the glyph to
highlight — not to every glyph that happens to draw the same entity.

**Why this priority**: it is the interaction Adam named, but it depends on
User Story 2's mapping existing first, and the benchmark does not need it.

**Independent Test**: a glyph identifier resolves to nodes, and a node
resolves back to glyph identifiers, for every glyph in a diagram.

**Acceptance Scenarios**:

1. **Given** a glyph, **When** it is resolved, **Then** the nodes derived
   from that glyph are returned — not the nodes of every glyph sharing its
   entity.
2. **Given** an entity drawn several times, **When** each of its glyphs is
   resolved, **Then** the results are distinguishable where the generated
   nodes differ, and where they are genuinely the same node that collapse is
   reported rather than hidden.
3. **Given** a solved node, **When** it is resolved for display, **Then**
   the glyphs to highlight are returned, empty if it derives from nothing
   drawn.
4. **Given** an entity drawn in a diagram but not represented at all,
   **Then** asking for it returns an explicit "not represented" with a
   reason, so the interface can grey it out rather than appear broken.

---

### User Story 4 — Choose how a set's members combine (Priority: P2)

A set's value has to come from its members somehow, and the choice changes
predictions. Measure the candidates rather than asserting one.

**Why this priority**: User Story 1 needs *a* rule to function and a
defensible default exists, so this is refinement, not a blocker.

**Independent Test**: each rule is benchmarked on the same cases and
reported with macro-F1 and a per-pathway breakdown.

**Acceptance Scenarios**:

1. **Given** the candidate rules, **When** each is benchmarked, **Then**
   each is reported with macro-F1, per-pathway net change, and the count of
   changed predictions where both arms converged.
2. **Given** a rule is adopted, **Then** the arithmetic justifying it is
   recorded, not only its score.

---

### Edge Cases

- **Cycles in set membership.** Reactome's set graph can contain them;
  recursion without a visited-set does not terminate.
- **A set with one member.** Must behave identically to the member itself.
- **A member unreachable from the perturbation.** It sits at baseline and
  drags an averaged set value toward "no change" — the recorded TP53
  signal-dilution failure. This is the main risk in User Story 4 and the
  reason a reachability-restricted variant is measured.
- **A glyph whose entity was never generated.** 38 per diagram measured.
  Must resolve to an explicit empty answer with a reason.
- **Glyph ids across Reactome releases.** Glyph identifiers are diagram
  layout artifacts and there is no evidence they are stable across releases.
  Anything persisted must record the release it came from. Version skew has
  already produced one false finding on this project.
- **Duplicate cofactor glyphs collapsing.** Nine ATP glyphs mapping to one
  node is correct behaviour, not a bug; the mapping must express "these nine
  are the same node" rather than appear lossy.
- **The count mismatch is not always benign.** One glyph to 32 nodes is
  positional decomposition, expected. One glyph to zero nodes may be
  expected or a defect; only a declared exclusion reason distinguishes them.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: For any Reactome entity in a generated pathway, the system
  MUST return every node representing it, each labelled with the
  relationship (the entity itself, a member of a set it was split into, a
  component of a complex, a variant, a collapsed duplicate) and the depth at
  which it was reached.
- **FR-002**: For any generated node, the system MUST return every entity it
  represents, with the same labelling.
- **FR-003**: Set resolution MUST recurse through nested sets to leaves, and
  MUST terminate on a cyclic membership graph.
- **FR-004**: Where a node derives from something drawn in a diagram, the
  system MUST record which glyph, so that a glyph resolves to its nodes and
  a node resolves back to its glyphs.
- **FR-005**: Any persisted mapping MUST record the Reactome release it was
  derived from, and MUST NOT be used against artifacts from another release
  without an explicit override.
- **FR-006**: An entity with no representation MUST appear in a declared
  exclusion list with a stated reason. Silent absence is not permitted.
- **FR-007**: A completeness check MUST run in both directions and MUST be
  demonstrated to fail on a deliberately corrupted mapping.
- **FR-008**: Resolving a set to a single value MUST combine positions of
  the same member before combining across members, so that how finely a
  member was split does not change its weight.
- **FR-009**: A partially resolved set MUST be reported as partial, naming
  what is missing, rather than being combined over the members that happened
  to resolve.
- **FR-010**: Each candidate combining rule MUST be selectable and MUST
  default to current behaviour until measured.
- **FR-011**: Columns that are declared MUST be populated or removed. A
  column that is empty in every row misleads its reader.
- **FR-012**: The number of cases the benchmark discards for want of a
  resolvable readout MUST be reported alongside any accuracy figure, since
  it sets the denominator.

### Key Entities

- **Entity**: a Reactome biological entity, identified by its stable
  identifier.
- **Glyph**: one drawing of an entity in one diagram, with its own
  identifier. An entity may have many; this is the identity Adam needs and
  the one currently discarded.
- **Node**: one node in a generated network, identified by its uuid.
- **Resolution**: one entity-to-node correspondence, carrying the
  relationship kind, the depth, and the glyph where one applies. Grouped one
  way it answers "what represents this entity"; grouped the other, "what
  does this node stand for".
- **Exclusion**: an entity deliberately not represented, with a reason —
  what makes completeness checkable rather than aspirational.
- **Set value rule**: how member values combine into one value for a set.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Every entity participating in a reaction of a generated
  pathway either resolves to at least one node or appears in the exclusion
  list with a reason — zero unexplained absences across all ten benchmark
  pathways.
- **SC-002**: Every generated node resolves to at least one entity — zero
  unexplained absences, both directions verified.
- **SC-003**: All 20 currently blocked set readouts resolve to their leaf
  members, and the number of benchmark cases discarded for an unresolvable
  readout falls from 204 to no more than 27.
- **SC-004**: Every glyph in a diagram resolves to its nodes or to an
  explicit empty answer with a reason, and every node deriving from a glyph
  resolves back to it.
- **SC-005**: An entity drawn more than once in the same compartment is
  distinguishable by glyph wherever the generated nodes differ, and where
  they do not, the mapping says so.
- **SC-006**: The completeness check fails when the mapping is deliberately
  corrupted.
- **SC-007**: Each candidate set-combining rule is reported with macro-F1
  and a per-pathway breakdown, including the rules that lose.
- **SC-008**: No declared column is empty in every row.

## Assumptions

- Glyph identifiers are unique within a diagram — verified on two diagrams
  (114/114 and 154/154) and assumed for the rest, checkable by the same
  measurement.
- Glyph identifiers are NOT assumed stable across Reactome releases; this
  drives FR-005. If they turn out to be stable, that is a simplification to
  take later, not one to rely on now.
- Reactome Release97 and the ten-pathway catalog are the basis; the mapping
  is generated per pathway alongside the network.
- Recovering these 204 cases changes the denominator of every published
  benchmark figure, so no number produced after this feature is comparable
  to one produced before it without restating both.
- Building the set-combining rule also supplies the aggregator the project
  constitution names as its open worked example — set-valued catalysts
  marked OR, costing 51 of 223 cases. That is not claimed as delivered here,
  only as unblocked.
- Reaction proxies stay out of scope as a substitute for set membership: a
  producing reaction's activity is a different quantity from its product's
  abundance, and substituting one for the other would inflate coverage
  without measuring what the case asks about.
