# Data Model: Every node maps, both ways

Phase 1. The central artifact is one table per pathway. Grouped by
`stable_id` it answers Reactome → LNG; grouped by `uuid`, LNG → Reactome.

## Resolution (the table)

Supersedes `node_reaction_context.csv`, which is the same idea with three
of its four roles' rows broken (research.md R3) and no glyph or set
information.

| field | type | rule |
|---|---|---|
| `stable_id` | string | Reactome stable id. **Stable id is primary**; a database id may be carried in a separate column but never as the key. |
| `uuid` | string | the generated node |
| `relation` | enum | how this node represents this entity — see below |
| `depth` | int | 0 for `self`; the number of containment hops otherwise. Observed maximum for sets is 5. |
| `role` | enum \| empty | `input`, `output`, `catalyst`, `regulator` where the node was minted at a reaction position; empty otherwise |
| `reaction_stid` | string \| empty | the reaction the node was minted against; empty for boundary/collapsed nodes |
| `glyph_id` | int \| empty | the diagram glyph this node derives from; empty when nothing draws it |
| `diagram_stid` | string \| empty | which diagram `glyph_id` belongs to — a glyph id is only unique **within** a diagram |
| `release` | string | the Reactome release, on every row (FR-005) |

`(stable_id, uuid, relation, role, reaction_stid)` is the natural key. The
same pair can legitimately appear more than once with different relations —
an entity that is both an input here and a catalyst there — and FR-003
requires all of them, not one chosen silently.

### `relation` values

Derived from the code paths that actually mint nodes, not invented:

| value | meaning |
|---|---|
| `self` | the node *is* this entity |
| `set_member` | the node is a member the entity was split into; `depth` ≥ 1 |
| `complex_component` | the node is a component of this complex |
| `variant` | a `set_variant` node — one chosen combination of a set-valued participant |
| `reaction` | a virtual-reaction node standing for a reaction |
| `boundary_collapsed` | one node shared across positions via the boundary cache |
| `dissociation_sink` | a terminal sink node |

Any node whose relation cannot be determined is a defect, not an `other`
bucket. There is deliberately no catch-all value.

## Exclusion

| field | rule |
|---|---|
| `stable_id` | the entity with no node |
| `reason` | free text, required, non-empty |
| `release` | as above |

**Expected to be empty.** R5 measured every currently-unresolved entity as
an EntitySet (29 of 29 in PIP3, 9 of 9 in Cell Cycle Checkpoints), so once
set membership is mapped there is no residual unmappable class. A non-empty
exclusion list is a finding to investigate, not a category to accept — this
is what makes SC-001 a gate rather than a formality.

## Glyph

Read from the diagram layout, not invented.

| field | rule |
|---|---|
| `glyph_id` | unique **within** its diagram — verified 114/114 and 154/154 |
| `diagram_stid` | required alongside `glyph_id`; without it the id is meaningless |
| `entity_stable_id` | what the glyph draws |
| `reaction_stid`, `role` | the position it occupies |

**`(reaction_stid, entity_stable_id, role)` uniquely identifies a glyph** —
0 of 156 diagram triples resolved to more than one glyph. This is what lets
an entity drawn nine times be disambiguated: it is drawn once per reaction.

## Set value

How a set readout becomes one number, two levels so that decomposition
granularity cannot skew it (FR-008).

| level | over | rule |
|---|---|---|
| 1 | positions of the *same* member | combine first — `R-HSA-202074`'s members map to 3, 2 and 2 uuids, and without this AKT1 outweighs AKT2 by decomposition alone |
| 2 | across members | the rule under test |

Level-2 candidates, all measured (US4), default unchanged until then:

| rule | meaning | note |
|---|---|---|
| `mean` | pool semantics | with uniform baseline x₀, `Σxᵢ/(n·x₀) = mean(xᵢ)/x₀`, so mean **is** the sum-of-abundances reading; matches Adam's "OR should be an average" |
| `max` | up if any member is up | today's multi-uuid default; expected to over-call UP |
| `mean_reachable` | mean over members reachable from the perturbation | guards the dilution failure where inert members drag the set toward no-change |

A set that resolves only partially is reported partial and **not**
combined over whatever resolved (FR-009).
