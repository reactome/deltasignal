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

## R8 — Set readouts recovered: a validity win, not an accuracy win (T013–T017)

Coverage on the ten-pathway catalog, same networks (cat7 is structurally
identical to the original catalog, 10/10 on edge counts, node counts and the
edge_type/sign/and_or/stoichiometry signature):

| | scored | of 847 |
|---|---|---|
| baseline | 564 | 66.6% |
| **with set resolution** | **742** | **87.6%** |

Remaining unscored: 10 `set_partially_resolved`, 3 `absent_from_network` —
13 against the SC-003 target of ≤27. 194 cases resolve through `set_members`.

**Scores, and the denominators are different so the rates are not
comparable:**

| arm | correct | accuracy | macro-F1 |
|---|---|---|---|
| baseline | 365/**564** | 0.6472 | **0.5781** |
| with set resolution | 491/**742** | 0.6617 | **0.5771** |

**Macro-F1 is flat — 0.5781 to 0.5771.** On the metric this project
optimises, recovering 178 cases changes nothing. The accuracy rise is a
coverage effect, not better prediction, and quoting it as an improvement
would be the denominator mistake the plan warned about.

The recovered cases are not harder than the ones we were already scoring —
**123 of 178, 69.1%**, against the baseline's 64.7% — and their class
balance is even (84 UP, 85 DOWN, 9 NO_CHANGE). 116 of the 178 are PIP3.

**What this actually fixes is validity, and that matters more than the
score.** Every number this project has published was computed on 564 of 847
cases, and the discarded third was not random: it was systematically the
set-shaped readouts — phospho-AKT, phospho-ERK, phospho-FOXO — of the
best-characterised pathways in the benchmark. That is a selection effect in
every published figure, and it is now closed.

**Noise check.** 20 previously-scored cases changed prediction. All 20 are
`exact`-mode cases that set resolution cannot touch, **0 of 20 converged in
both arms**, and they sit in TP53 (17) and ERBB2 (3). That is the known
uuid4 sweep-order signature, seen twice before with the same shape, arising
from regeneration rather than from this change. Net +3, i.e. noise.

### Restated comparison on the corrected case set (742 cases)

| model | correct | accuracy | macro-F1 |
|---|---|---|---|
| DeltaSignal on our networks | 491/742 | 0.6617 | 0.5771 |
| shortest signed path | 518/742 | 0.6981 | 0.6060 |
| **MP-BioPath published** | **546/742** | **0.7358** | **0.6682** |

And on the 740 the control also scores:

| model | correct | accuracy | macro-F1 |
|---|---|---|---|
| DeltaSignal on our networks | 490/740 | 0.6622 | 0.5776 |
| **DeltaSignal on MP-BioPath's networks** | **538/740** | **0.7270** | 0.6448 |
| MP-BioPath published | 544/740 | 0.7351 | 0.6678 |

**Feature 004's control survives the bigger case set intact**: given the same
networks DeltaSignal is 6 cases behind MP-BioPath (538 vs 544) and 48 ahead
of itself on our networks (538 vs 490). The propagator is still exonerated
and the gap is still network structure — now measured on 740 cases instead
of 562, so it is no longer resting on the narrower set.

DeltaSignal still loses to a model-free signed traversal (491 vs 518), which
remains the sharpest statement of what is unfixed.

## R9 — NEGATIVE: the set-combining rule does not matter, and mean is slightly worse

Three arms, identical catalog, identical 742 scored cases:

| rule | correct | accuracy | macro-F1 | F1_DOWN | F1_NC | F1_UP |
|---|---|---|---|---|---|---|
| **`max`** (default, current behaviour) | **491/742** | 0.6617 | **0.5771** | 0.748 | 0.265 | 0.718 |
| `mean` | 488/742 | 0.6577 | 0.5740 | 0.748 | 0.259 | 0.716 |
| `mean_reachable` | 488/742 | 0.6577 | 0.5740 | 0.748 | 0.259 | 0.716 |

**I predicted mean would win and it does not.** The arithmetic argument in
data-model.md is still sound — with uniform baseline x₀, a set's pool fold is
`Σxᵢ/(n·x₀) = mean(xᵢ)/x₀`, so mean *is* the sum-of-abundances reading, and it
matches Adam's stated rule that OR configurations average. The data simply
does not support acting on it.

**Why the rule barely matters — the mechanism, not just the null.** Of the
178 scored set cases, **131 (73.6%) have `max` exactly equal to `mean`**,
because every member carries an identical value. Set members are typically
paralogs sitting at structurally symmetric positions receiving the same
upstream signal, so they land on the same number and there is nothing for a
combining rule to choose between. 44 cases differ by more than 10% in value,
but only **6 cross a classification cutoff** — all in RAF_MAP_kinase_cascade,
all converged in both arms, net **−3**.

`mean_reachable` is **byte-identical to `mean`**. It restricted the member
set on 13 cases (S_Phase 8, TP53 4, WNT 1) and changed **no** prediction, so
the dilution failure it was designed to guard against — inert members
dragging an average toward no-change — **does not occur on this case set**.
That hypothesis is answered, not merely untested.

**Decision: keep `max`.** Not because it is more principled — it is less so —
but because changing a default requires evidence and there is none. Six
discriminating cases cannot separate two rules, and what little signal exists
points the other way. `mean` and `mean_reachable` stay available as flags,
and this is recorded so the next person does not re-run it blind.

**What would change this.** A case set where set members genuinely diverge.
If OR-cluster loss is ever modelled (the constitution's open worked example,
where a set-valued catalyst marked OR costs 51 of 223 cases because
`max(and, or)` discards the loss entirely), members would stop moving in
lockstep and the rule would start to matter. Re-run then, not before.
