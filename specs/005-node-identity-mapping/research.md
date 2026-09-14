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

---

# Phase 3 findings (T008)

## R6 — Recursion resolves all 20 blocked readouts

Validated against the real blocked readouts rather than fixtures:

| | one hop | recursive |
|---|---|---|
| fully resolved | 15 | **18** |
| partially resolved | 2 | 2 |
| unresolvable | **3** | **0** |

Depth distribution 15×1, 4×2, 1×3 — recursion is load-bearing for 5 of 20,
which is exactly the nested-set count predicted in R1.

Of the 204 blocked cases, **194 are fully resolvable**; the remaining 10 sit
in the two partial sets and, under FR-009, are not scored by combining over
the members that happened to resolve.

## R7 — "227 dropped set members" is 80% by design. Nearly reported as a bug.

Across the catalog, **227 of 841 set-member leaves (27%) have no node**, in
every pathway, 0% to 66.7%. That reads as a large generation defect, and the
granularity hypothesis — that a phospho-form leaf is represented by its base
protein — was tested and **failed completely: 0 of 227 have a same-gene node
present.** With both of those in hand the obvious conclusion was a serious
bug.

**The tell was in the counts.** Almost every pathway showed exactly 14 or 28,
which is not what scattered losses look like. They are the same entities
everywhere: `R-HSA-68524` "Ub" and `R-HSA-113595` "Ub [cytosol]", whose
members are the individual UBB/UBC/UBA52/RPS27A repeat units —
`UBB(1-76)`, `UBC(153-228)` and so on. Those are the **deliberately atomic
modifier sets** from the modifier-collapse work, and the generator publishes
the authoritative list as `get_modifier_isoform_entity_set_ids()` (46 sets).

Splitting on it:

| pathway | missing | intentionally atomic | genuine |
|---|---|---|---|
| Cell_Cycle_Checkpoints | 28 | 28 | **0** |
| PIP3_activates_AKT_signaling | 50 | 28 | **22** |
| Signaling_by_ERBB2 | 26 | 14 | **12** |
| Mitotic_G1_phase_and_G1_S_transition | 18 | 14 | 4 |
| Transcriptional_Regulation_by_TP53 | 17 | 14 | 3 |
| Signaling_by_WNT | 30 | 28 | 2 |
| HDR / RAF | 15 each | 14 each | 1 each |
| S_Phase | 28 | 28 | **0** |
| Mitotic_Prophase | 0 | 0 | **0** |
| **TOTAL** | **227** | **182 (80%)** | **45** |

**The real gap is 45, not 227**, concentrated in PIP3 (22) and ERBB2 (12).
Those are Complexes that reach their reactions only *via an EntitySet
participant* — never as a direct input or output — and structurally
identical siblings differ in whether they got a node (`R-HSA-1963593` and
`R-HSA-1248703` have nodes; `R-HSA-1963583` and `R-HSA-1250316` do not).
That inconsistency is a genuine defect and it is what makes the two partial
readouts partial.

**Consequences for the plan:**

1. `node_exclusions.csv` will **not** be empty, contradicting R5's
   expectation. R5 measured *direct* input and output participants and was
   right about those; set members are a layer it never looked at. Roughly
   182 entries are the intentional modifier sets and need the reason
   `atomic_modifier_set`; about 45 need investigating.
2. An exclusion reason is doing real work here rather than being a
   formality — it is the only thing separating a design decision from a bug
   in the same list.
3. The 45 are a separate defect from this feature. File, do not absorb.

**Method note.** Two hypotheses were tested and rejected (granularity, 0 of
227) before the right one was found, and the right one was found by looking
at the *shape* of the numbers rather than their size. A count alone would
have shipped "227 members are being dropped", which is true and useless.
