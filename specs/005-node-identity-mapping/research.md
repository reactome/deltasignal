# Research: Every node maps, both ways, including on the diagram

Phase 0. Five questions, all resolved by measurement. **One of them
retracts a claim I put in the spec**, which is recorded first because the
spec is wrong until amended.

## R1 — RETRACTION: Reactome's set-membership graph has NO cycles

**The spec asserts cycle protection is required "because Reactome's set
graph contains cycles". That is unverified and false.**

Measured on Release97:

| check | result |
|---|---|
| EntitySets whose member is an EntitySet | 5,440 |
| direct self-membership | **0** |
| membership cycles of length 2 / 3 / 4 | **0 / 0 / 0** |
| max nesting depth, set → leaf | **5** |

Nesting is real and deep enough that single-hop resolution is wrong (which
is the substantive point, and it stands). Cycles are not.

**Decision**: keep a visited-set guard anyway, but justify it as cheap
insurance against a future curation error, **not** as a response to
observed cycles. Recursion depth is bounded at 5 today, so a depth cap is a
reasonable secondary guard.

**Why this matters beyond the fix**: the claim went into a spec as fact on
the strength of it sounding likely. The constitution's second principle
exists because of exactly this. Amend `spec.md`'s edge-case entry.

## R2 — The provenance record already exists, keyed positionally

**Decision**: extend the existing export rather than invent a new artifact.

`logic_network_generator.py:472` mints node ids through a registry keyed
`(entity_dbId, reaction_uuid, role)`. **A node's identity is already
positional — entity × reaction × role.** That triple is precisely the
provenance this feature needs, and it is fully materialised in memory at
generation time; nothing needs re-deriving.

It is even already exported, as `node_reaction_context.csv`
(`context_node, reaction_id, role`). What is missing from it is the glyph
id, the entity's stable id, and any relation other than the four roles.

**Alternatives considered**: a fresh `node_identity.csv` built from
scratch — rejected, it would duplicate a partially-correct export and leave
two sources of truth. Reconstructing the mapping in DeltaSignal from
`nodes.csv` — rejected under constitution principle I and because the
information is lossy by the time it reaches the consumer.

## R3 — BLOCKER: that export is 31.3% orphaned, and it is the causal 31%

`node_reaction_context.csv` rows whose `context_node` does not appear in
`logic_network.csv` at all, across the ten-pathway catalog:

| role | rows | orphaned | |
|---|---|---|---|
| input | 3,641 | 0 | 0.0% |
| output | 3,491 | 0 | 0.0% |
| **catalyst** | **2,228** | **2,228** | **100.0%** |
| **regulator** | **1,018** | **1,018** | **100.0%** |
| TOTAL | 10,378 | 3,246 | 31.3% |

This is LNG issue #67, previously filed and now quantified. Every catalyst
and every regulator row points at a node that does not exist in the
network — the export writes the fetch-row uuid rather than the endpoint
`append_regulators` actually used.

**Decision**: fixing #67 is a **blocking prerequisite**, not a related
cleanup. Catalysts and regulators are where the causality lives, and R4
shows they are also most of the glyph-join shortfall.

## R4 — Glyph identity: `(reaction, entity, role)` is a unique key

**Decision**: join glyphs to nodes on `(reaction, entity, role)`. No new
identifier scheme is needed.

Measured on R-HSA-1257604:

- 156 `(reaction, entity, role)` triples in the diagram, and **0 of them
  resolve to more than one glyph.**

That is the whole answer to Adam's question. An entity drawn nine times is
drawn *once per reaction*, so the triple that already keys our node
registry also uniquely identifies the glyph. "Which one did the uuid come
from" is answerable exactly, with no ambiguity to resolve.

Current join coverage, before any fix:

| | triples | |
|---|---|---|
| joinable today | 114 | 73.1% of diagram, 21.7% of LNG |
| diagram-only | 42 | **25 catalyst**, 9 input, 8 output |
| LNG-only | 411 | generation descends below the diagram |

The 25 catalyst misses are R3's orphaning, so fixing #67 should take
diagram coverage from 73% to roughly 89% without any new joining logic.
The 411 LNG-only triples are expected — a diagram draws only its own
top-level reactions — and the point of recording them is that "expected"
becomes checkable rather than assumed.

## R5 — The exclusion list is empty. Every gap is a set.

**Decision**: target zero exclusions, not "few".

Entities participating as input or output of any reaction in the pathway,
against what currently resolves to a node:

| pathway | participating | resolvable | unresolved |
|---|---|---|---|
| PIP3 (R-HSA-1257604) | 157 | 128 | 29 (18.5%) |
| Cell Cycle Checkpoints (R-HSA-69620) | 101 | 92 | 9 (8.9%) |

**All 29 and all 9 are EntitySets** — 15 CandidateSet + 14 DefinedSet, and
6 DefinedSet + 3 CandidateSet respectively. There is no residual class of
genuinely unmappable entity.

So the exclusion list required by FR-006 should come out **empty** once set
membership is mapped, and any non-empty entry is a defect to investigate
rather than a category to accept. That makes SC-001 a real gate.

## Open questions carried into implementation

1. Are glyph ids stable across Reactome releases? Unresolved — only the
   Release97 diagrams are on disk. The design assumes not (FR-005).
2. Does the 89% projected glyph coverage materialise once #67 is fixed, and
   what are the residual 17 input/output misses?
3. Which set-combining rule wins? Unmeasured by design; US4.
