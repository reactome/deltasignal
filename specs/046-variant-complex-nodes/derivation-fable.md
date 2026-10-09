# specs/046 — blind derivation B (Fable), 2026-10-09

Independent derivation of the variant key and the end-to-end variant wiring,
from the requirement in `spec.md` (I1–I5) and the generator code at
`logic-network-generator` main `6990015`. I did not read any other
`derivation-*.md`. Every point I am not sure of is marked **UNSURE**.

Sources: `spec.md`; generator `CLAUDE.md`, `docs/UUID_DESIGN.md`;
`src/reaction_generator.py` (`break_apart_entity`, `get_broken_apart_ids`,
`get_uids_for_iterproduct_components`, `decompose_by_reactions`, `CAPPED_IDS`);
`src/best_reaction_match.py`; `src/logic_network_generator.py`
(`_register_phase1`, `_get_or_create_entity_uuid`, `_complex_variant_leafsets`,
`_matching_leaves`, `_map_annotated_entity_to_nodes`, `_resolve_vr_entities`,
`_decompose_regulator_entity`, `_expand_complex_variants`, `append_regulators`,
`_emit_vr_inputs`, `create_pathway_logic_network`,
`_emit_boundary_decomposition_edges_inner`, `_node_leaves`, `export_nodes`,
`_derive_sets_and_chosen`, `find_pools`); `src/neo4j_connector.py`;
`~/deltasignal-catalogs/analysis/046/excess_classify3.py` and `size_census.py`.
Neo4j Release97, read-only, bounded queries only (one level at a time, or
`hasEvent*1..6` which is not a set/member path).

---

## 0. What the code does today that the spec does not say

These are the facts the design rests on. Each was read from the code or
queried; none is inferred.

**0.1 The matching-layer combination uid is a content hash over the chosen
components, not scoped by the complex.** `get_uids_for_iterproduct_components`
computes `sha256(str(sorted(component_to_input_or_output.values())))`, where
the values are the chosen leaf stIds and the nested combination uids. The
complex stId is stored in the row (`reactome_id`) but is not in the hash, so
two complexes with identical chosen content get the same uid (the docstring of
`_resolve_vr_entities` admits this). The hash is sorted, so it is seed-stable.

**0.2 A set at stoichiometry n > 1 inside a complex is ONE slot.**
`get_complex_components` returns `{member: stoichiometry}`, so a set that
appears twice is one entry with stoichiometry 2; `break_apart_entity` calls
`break_apart_entity(set)` once and the product has one factor for it. Both
copies therefore choose the SAME member, and the stoichiometry is attached to
the chosen member (`_direct_component_stoichiometry`). Heteromers of a set
are not representable. Neo4j examples (a bounded query another census listed;
I verified the FGFR2 one myself): CaMKII [endocytic vesicle membrane]
R-HSA-417004 = 4 × set "CaMKII subunits delta/gamma" (2 members); PKA tetramer
R-HSA-111922 = 2 × "PKA regulatory subunit" (4 members); AKAP9:KCNQ1:KCNE dimer
R-HSA-5577100 = 2 × KCNEs (5); 3×CLTC:3(CLTA,B) R-HSA-5138447; and in FGFR2
every "Activated FGFR2b homodimer bound to FGF" member complex (e.g.
R-HSA-190411) holds "FGFR2b-binding FGFs" R-HSA-189967 (5 members) at
stoichiometry 2.

**0.3 Two incompatible `::variant::` conventions already exist, both dead
under the default flags.** `_map_annotated_entity_to_nodes` (under
`LNG_COMPLEX_AS_NODE=0`) names a variant by the complex stId plus the FLAT
SORTED SET OF ALL TERMINAL LEAVES (`_complex_variant_leafsets`, built with
`get_terminal_components`, so a chosen member that is itself a complex is
dissolved to its leaves). `_expand_complex_variants` (negative regulators,
same flag) names it by the complex stId plus the sorted CHOSEN MEMBER IDS,
with nested complexes contributing their own `::variant::` ids. Neither
records which set a chosen member came from. The flat-leaf form collides for a
set that holds a monomer and a dimer of the same protein — "p-T,Y MAPK
monomers and dimers" R-HSA-5674340, a catalyst of three RAF-pathway reactions
(R-HSA-5674496, 5675194, 5675198), is literally that.

**0.4 F9 is both a determinism problem and a semantic one.**
`decompose_by_reactions` hands `find_best_reaction_match` the input and output
combinations as Python sets; `find_best_match_both_decomposed_reactions` does
`list(...)` on them, so the row/column order — and therefore which of several
equal-cost assignments the Hungarian solver returns — follows the string hash
seed (`bin/create-pathways.py` pins `PYTHONHASHSEED`, so it is reproducible,
not meaningful). The overlap count is over
`component_id_or_reference_entity_id`, which `get_reference_entity_id`
resolves to an HGNC-linked `ReferenceEntity` for an EWAS and to the entity's
own stId otherwise. A set member that is a simple complex therefore carries no
reference entity: in "RAF phosphorylates MAP2K dimer" R-HSA-5672978 the input
slot "MAP2K homo/heterodimers" R-HSA-5672716 = {MAP2K1 dimer, MAP2K2 dimer,
MAP2K1:MAP2K2} and the output slot "p-2S MAP2K homo/heterodimers"
R-HSA-5672721 = {p-MAP2K2 dimer, p-MAP2K1 dimer, p-hetero} overlap by ZERO on
every pairing, so the three pairings are an exact tie. (That reaction is
capped today, so the tie is moot there; it is not moot in general.)

**0.5 Phase 2 unifies only across `reaction_connections` pairs.**
`_get_or_create_entity_uuid` is called for `(p_vr, f_vr)` taken from
`reactome_to_vr` over the rows of `reaction_connections` (precedingEvent,
unioned with diagram connectivity by the caller). Two reactions that share a
produced entity but have no such pair keep positional copies by design
(`docs/UUID_DESIGN.md`; `_register_phase1` docstring: "bridging those was
measured harmful four times"). So I2 as worded ("the same node wherever it
appears") is STRONGER than what the generator does for a plain entity today.
See §3.2.

**0.6 Catalyst == input is common, and the catalyst's sets are the input's
sets.** R-HSA-5672978: catalyst = R-HSA-5672720 = its own substrate complex.
FGFR2 (R-HSA-5654738): "Activated FGFR2 phosphorylates FRS2" R-HSA-5654397
(catalyst = input "Activated FGFR2:FRS2" R-HSA-5654178), likewise SHC1, PLCG1,
FRS3, PTPN11, both autophosphorylations (catalyst = input SET "FGFR2b homodimer
bound to FGF" R-HSA-192615), and PP2A/SPRY2 — eight of the twelve
set-bearing-catalyst reactions I listed in FGFR2 are of this shape. A design
that multiplies catalyst variants by input variants would make, e.g., the
FGFR2b-copy of FRS2 phosphorylation require the FGFR2c-variant catalyst.

**0.7 Output sets are mostly the same set stIds as the input sets.** RAF
R-HSA-5672972 "MAP2Ks and MAPKs bind to the activated RAF complex": inputs
MAPKs R-HSA-169291 (2), MAP2K dimers R-HSA-5672716 (3), scaffolds
R-HSA-5672717 (CandidateSet, 8), and complex R-HSA-5672718 (4 × 4 × 5 = 80
variants over "activated RAF/KSR1", "mature p21 RAS", "'activator' RAFs");
output complex R-HSA-5672720 has EXACTLY those four as its components. Product
2 × 3 × 8 × 80 = 3,840 > 512. In R-HSA-5672978 the output R-HSA-5672723 keeps
three of the four slots by stId and replaces R-HSA-5672716 with
R-HSA-5672721; the members pair one-to-one by reference-entity multiset
(MAP2K1 dimer {MAP2K1×2} ↔ p-MAP2K1 dimer {MAP2K1×2}, etc.).

**0.8 The F10 anatomy, on FGFR2.** "Activated FGFR2" R-HSA-5654152 is a
DefinedSet of two DefinedSets (R-HSA-192606, R-HSA-192616), each of two
complexes, each complex = HS + 2 × FGFs-set + 2 × p-8Y-FGFR2 (e.g.
R-HSA-190411). R-HSA-190408 "Autocatalytic phosphorylation of FGFR2b" OUTPUTS
the bare set R-HSA-192606 and precedes R-HSA-5654399 "Activated FGFR2 binds
FRS2", which INPUTS the bare set R-HSA-5654152. Today both bare sets go
through the set branch of `_map_annotated_entity_to_nodes`, whose
`_matching_leaves` dissolves each chosen member complex into {HS, FGFx,
p-8Y-FGFR2}. The consumer copy therefore requires HS and FGFx as separate root
inputs, and the produced complex (whole) is never a node. Note the nuance for
I3's counting: in this case the signal still conducts through the shared
p-8Y-FGFR2 leaf (both sides dissolve to it and Phase 2 merges it), so the
"hand-off cut" is a cut of the complex's identity and of the FGF/HS members,
not always a total disconnection. `excess_classify3.py` counts the root
subunit; I suggest it also report whether the copy has ANY produced input
from the same holder (§5.4). **UNSURE** how many of the 2,504 are total cuts.

**0.9 Consumers of the id format.** `split("::variant::")[0]` is the parent
everywhere (also the solver, `reaction_model.jl:206`, which reads nothing
else). The tail is parsed as `_`-joined stIds by `_parse_variant_members`,
`_node_leaves`, `export_nodes` (then `_derive_sets_and_chosen` re-derives sets
and chosen members from Neo4j by leaf intersection), `find_pools.parts`,
`export_node_resolution`, `export_drugs`. The boundary layer's `_is_complex`
and `_is_set` return False for any id containing `::variant::`, so a variant
ROOT would today be left undecomposed, silently. The benchmark resolves a gene
to nodes through `nodes.csv` `member_leaves` (`benchmark_mpbiopath_cases.py`
registers every node under its own stable id and every leaf), so a wrong
`member_leaves` on a variant node is a coverage loss that would read as an
accuracy loss.

**0.10 Catalyst/regulator node lookup is "first positional copy".**
`append_regulators` builds `stid_to_existing_uuid` from the registry in
insertion order and takes the first uuid per stId; a catalyst produced by the
preceding reaction may thus attach to a different positional copy than the
one that reaction produced. Pre-existing, independent of variants.

**0.11 Cap bookkeeping mixes levels.** `get_broken_apart_ids` adds the
`reactome_id` it was called with to `CAPPED_IDS`: the COMPLEX stId when called
from `break_apart_entity`, the REACTION stId when called from
`decompose_by_reactions`. `_resolve_io` tests `reaction_id in CAPPED_IDS`;
`_complex_variant_leafsets` keeps its own `_variant_capped`. Three caps, two
registers; they agree today only because each is a pure function of the
entity and `MAX_VARIANTS`.

---

## 1. The variant key

### 1.1 Definition

A **choice** σ is a map from *slots* to chosen members. A slot is a set
occurrence inside a participant, identified by the set's stId (see 1.3 for
the nested-path qualification). The chosen member of a slot is always a
non-set entity: nested sets are flattened first (a set whose member is a set
contributes that inner set's members), exactly as `break_apart_entity` does.

```
VK(e, σ):
  e is a leaf (EWAS, SimpleEntity, GenomeEncodedEntity, Drug, Polymer,
      OtherEntity, Cell, or unknown label)                    -> e
  e is a modifier-isoform set (modifier_isoform_set_ids())   -> e
  e is a set S                                               -> VK(σ[S], σ)
  e is a complex with no set anywhere inside it              -> e
  e is a complex C with |variants(C)| > MAX_VARIANTS         -> C
  e is a complex C, otherwise                                ->
      C + "::variant::" + "_".join(sorted(
          f"{slot}={VK(σ[slot], σ)}"  for slot in slots(C) ))
```

where `slots(C)` is every set occurrence reachable from C through
`hasComponent` and through the chosen members' own `hasComponent` (a chosen
member that is a complex with sets contributes its slots too), and `variants(C)`
is the product over `slots(C)` of member counts (the same number
`_complex_variant_leafsets` computes, so "capped" is decided identically
wherever C appears).

Worked examples (Release97):

- Activated FGFR2:FRS2 R-HSA-5654178, choosing FGFR2b long dimer with FGF7:
  `R-HSA-5654178::variant::R-HSA-189967=R-HSA-<FGF7>_R-HSA-5654152=R-HSA-190411`
  (slot 5654152 is the outer set-of-sets; its chosen member is the complex
  R-HSA-190411, which has no further variant key of its own because its only
  set, R-HSA-189967, is listed as a slot of the OUTER key — see 1.3).
- A bare-set input "MAPKs" R-HSA-169291 choosing MAPK1: `R-HSA-59282`
  (the set dissolves; the node is the species, as today).
- R-HSA-5672720 (3,840 variants): `R-HSA-5672720` (capped).

### 1.2 Why chosen MEMBERS, labelled by SLOT, and not leaves

- Leaves collide: a set holding the monomer and the dimer of one protein
  (R-HSA-5674340, §0.3) gives both variants the same leaf set. Members do
  not collide: distinct members have distinct stIds.
- Unlabelled members collide when two slots overlap in membership: C with
  S1 = {A, B} and S2 = {A, B, D} (two different sets of the same proteins —
  e.g. the inner RAF sets "activated RAF/KSR1" R-HSA-5672713 and "'activator'
  RAFs" R-HSA-5672709 share ARAF/RAF1/BRAF/KSR1 by gene but not by stId, so
  they do not collide today; a case that overlaps by stId would). Labelling
  by slot makes (S1→A, S2→B) and (S1→B, S2→A) distinct.
- A chosen member that is a complex contributes its stId, not its leaves,
  so the member complex keeps its identity. This is what fixes F10: the
  consumer of "Activated FGFR2" and the producer of "Activated FGFR2b
  homodimer bound to FGF" both name the chosen member R-HSA-190411 and the
  chosen FGF.

### 1.3 Nested sets, sets of complexes containing sets, and slot naming

- **Set of sets** (R-HSA-5654152 ⊃ R-HSA-192606 ⊃ complexes): the slot is the
  OUTERMOST set on that path from the participant (5654152); its chosen
  member is the flattened leaf of the set tree (the complex R-HSA-190411).
  The inner set id does not appear in the key. Reason: a producer that outputs
  the INNER set (R-HSA-190408 outputs R-HSA-192606) and a consumer that inputs
  the OUTER set must still name the same variant, and they do, because for a
  BARE set participant the set dissolves (VK(S, σ) = VK(member)) — the key is
  the member's key `R-HSA-190411::variant::R-HSA-189967=<FGF>` on both sides.
  The slot label matters only INSIDE a complex key, where it is the outer
  complex's own set occurrence, which both sides see identically because they
  see the same complex stId.
- **Set member that is a complex with sets**: its slots are hoisted into
  the enclosing complex key as `slot=member` pairs (1.1), so the key is flat
  (one `::variant::`, `_`-separated `slot=member` tokens, no nesting). A bare
  set participant whose chosen member is such a complex gets that complex's
  key. This keeps every consumer's parse (`split("::variant::")`, split on
  `_`) a one-level parse; only the token format changes (`=` inside tokens).
- **Path qualification.** If one participant contains the same set stId at
  two different `hasComponent` paths (C = D:E with D ⊃ S and E ⊃ S, D ≠ E),
  the two occurrences are different slots and must be told apart: the slot
  label is then `D/S` and `E/S` (the path of complex stIds below the
  participant, `/`-joined). When a set occurs once, the label is the bare
  set stId. **UNSURE** whether this case exists at all; the census (§5.5)
  should count it. The simple rule "label = set stId" is what the shared-slot
  logic in §2 needs, so path qualification is the exception, not the rule.
- **Stoichiometry > 1 of one set (0.2)**: one slot, one chosen member, the
  stoichiometry carried on the edge, exactly as today. This is the only
  choice that keeps the output complex's key derivable from the input's
  (Reactome's own output complex has the same single set at stoichiometry n,
  so it cannot express a heteromer either). Extension, NOT proposed now:
  n sub-slots `S#1..S#n` with an unordered (sorted) assignment, giving
  C(m+n−1, n) variants; the census should report how many copies that would
  add before anyone tries it.
- **CandidateSet vs DefinedSet**: identical (both are "one of"; the
  connector already reads `hasMember|hasCandidate` as one relation). The RAF
  scaffolds set R-HSA-5672717 has 7 candidates and 1 member; nothing in the
  key distinguishes them. A later arm could weight candidates differently in
  the solver, but that is a solver matter and the key should not encode it.
- **Modifier-isoform sets**: atomic (`VK = set stId`), decided by the same
  `modifier_isoform_set_ids()` on both sides. Unchanged.
- **OpenSet**: `get_set_members` treats it like any set; `_is_set` in the
  boundary layer lists it explicitly. Same here.

### 1.4 Determinism

The key is a pure function of (participant stId, σ) built from sorted
strings; no uuid4, no dict order, no hash seed. The σ enumeration (§2) sorts
member lists (`sorted(get_set_members(...))`) before `itertools.product`, so
copy ORDER is deterministic as well, which matters for `_register_phase1`'s
first-wins cache and `append_regulators`' first-existing lookup. Test:
regenerate under two `PYTHONHASHSEED`s and compare the multiset of edge
(source VK, target VK or reaction stId, type) — must be identical (§5.2 T6).

### 1.5 The cap in the key

"Capped" is decided per COMPLEX stId from its own `variants(C)` against
`MAX_VARIANTS`, so VK(C) = C on every side at once (producer, consumer,
catalyst, regulator). The REACTION-level cap (product over all slots of the
reaction > `MAX_VARIANTS`) does NOT change any participant's key; it changes
how many copies the reaction has and how the single copy is wired (§2.5).

---

## 2. Virtual reactions when inputs, catalysts and regulators are all expanded

### 2.1 One variable per set per reaction, not per participant

Let `V(R)` be the set of slot labels over ALL of R's inputs, outputs,
catalysts and positive regulators (negative regulators are handled in 2.4).
Two occurrences of the same set stId anywhere in R — bare input and inside
the output complex (R-HSA-5672972, all four slots), catalyst and input
(R-HSA-5672978, R-HSA-5654397), inside two different inputs — are ONE
variable. A copy is one assignment σ ∈ Π_{v ∈ V(R)} members(v).

Justification: §0.6 and §0.7. Independent choice would (a) multiply copies by
the catalyst's variants even when the catalyst IS the substrate complex, and
(b) require a matcher to re-pair the output's slots with the input's — which
is F9. Sharing by set identity removes both. The cost is one modelling
assumption: when a set appears twice in one reaction, the curator meant the
same member in both places. For dimerisations Reactome already encodes that
assumption (stoichiometry 2 of one set, §0.2); for input/output it is how the
output complex is built; for catalyst = input it is literal. **UNSURE** only
for the rare "same set in two DIFFERENT input complexes" case, which the
census should count (§5.5).

### 2.2 Output slots with no input counterpart: binding by member correspondence

When an output participant has a slot `w` whose set stId is not in any input,
catalyst or positive regulator (R-HSA-5672978: R-HSA-5672721 replaces
R-HSA-5672716), it is bound to an unmatched input-side slot `u` if the two
member lists pair one-to-one. The pairing rule, applied in order and recorded:

1. **Reference-entity multiset** of each member's terminal leaves (MAP2K1
   dimer → {MAP2K1:2}); a member of `u` and a member of `w` match when the
   multisets are equal. If every member of `w` has exactly one match in `u`,
   `w` is bound to `u` through that map and is no longer a variable.
2. If several members of `w` have the same multiset (two phospho-forms of one
   protein in the same set), pair them with the matching members of `u` by
   **canonical rank**: sort both tied groups by stId and pair by position.
   Count it (`slot_binding_ties`) and log the reaction. This replaces F9's
   hash-order tie with a stated, stable rule; it is still a guess, and the
   log is how anyone finds out whether it mattered.
3. If `|members(w)| ≠ |members(u)|` or some member has no match, the slots
   are **partially bound**: matched members bind, unmatched members of `w`
   are produced by no copy, and copies whose `u`-member has no `w`-match
   produce the output participant WITHOUT that slot's contribution — i.e.
   they are one-sided for that output. **UNSURE**; today's matcher pairs the
   surplus with the best-overlap output instead (`argmax`), which can make a
   copy produce a variant made from a member it did not consume. I prefer
   emitting nothing over emitting a wrong variant, and counting it
   (`slot_binding_partial`). Measurable as a two-arm question.
4. If `w` matches no input-side slot at all (a set that exists only in the
   output), it stays a free variable: copies then differ only in the output
   variant they make. Phase 1 sharing makes their input nodes one node, so
   this is "one reaction, OR over output variants", which is right.

With 2.1 and 2.2 the Hungarian step and the content-hash pairing are not
needed for variant reactions: the output keys are DERIVED from σ, so the
input→output correspondence is exact by construction. The content hash can
stay as the uid of `decomposed_uid_mapping` rows for provenance, but it must
not decide pairing.

### 2.3 Catalysts and positive regulators

For each copy σ, each catalyst c gets the node VK(c, σ) and each positive
regulator likewise (AND, as today). Slots of c that are not shared with any
input are variables in V(R) (2.1), so a catalyst set S = {E1, E2} on a
reaction with input set T = {A, B} gives four copies (A,E1), (A,E2), (B,E1),
(B,E2); their outputs for the same T-member coincide (VK depends on T only),
so the output node receives four `output` edges (`and_or = "or"`,
`entity_producer_count > 1`), which is the OR over copies the spec asks for.

This is the semantic change the spec pre-registers: today a set catalyst is
one pool combined by `product` (knock out E1 ⇒ every copy stops); under
expansion knocking out E1 stops the E1 copies and the OR keeps the E2 copies.

### 2.4 Negative regulators

Not variables. A negative regulator n is enumerated into its own variants
(§1), and the question is which variants inhibit which copies. Rule: a slot
of n that is also in V(R) is bound to the copy's choice (the FGFR2b-variant
of the feedback complex "Activated FGFR2:p-8T-FRS2" inhibits the FGFR2b copy
of R-HSA-5654397); slots of n not in V(R) are enumerated and every such
variant inhibits every copy (OR, `and_or = "or"`, as today). Rationale:
"any one inhibitor blocks" is reaction-level logic, not a copy-defining
choice; but an inhibitor that shares the substrate's identity with a copy
should not be made to inhibit the OTHER isoform's copy.

**UNSURE** on two counts. (a) A negative-regulator complex with k unshared
variants puts k inhibitor edges on every copy where today there is one node;
under `DS_INHIBITION_MODE=divide` each inhibitor divides, so k inhibitors at
baseline are harmless but k perturbed ones compound. (b) The solver's
self-inhibitor rule (`DS_SELF_INHIBITOR_WEIGHT`, specs/022/040) reads
`containment.csv`; variant inhibitor nodes need containment rows (§4.7) or
the rule goes inert on them, which would move predictions for a reason that
has nothing to do with this design. Both are testable before the arm.

### 2.5 The cap, and what "one OR pool per set" means under expansion

`copies(R) = Π_{v ∈ V(R)} |members(v)|` after binding (shared slots counted
once). If `copies(R) > MAX_VARIANTS`, R has ONE copy, wired as follows:

- **Bare set input / catalyst / positive regulator** → one pool node per set,
  fed by `set_member` (OR) edges from the member nodes (specs/045 for
  inputs, specs/033 for catalysts/regulators; unchanged). This is I4.
- **Complex participant that is itself capped** (VK = stId) → one node,
  as today. R-HSA-5672720 and R-HSA-5672723 are this.
- **Complex participant that is NOT itself capped but whose reaction is**
  (C1 with 30 variants × C2 with 30 = 900 > 512) — two sides:
  - as an INPUT: the single copy reads a pool node `C1 (any)` fed by
    `set_member` edges from every C1 variant node the pathway produces, plus
    the plain C1 root if none is produced. The pool stands for "any variant".
  - as an OUTPUT: the single copy writes `output` edges to every C1 variant
    node that some downstream copy consumes (lazily, by key), and to the
    plain C1 node if nothing consumes it. **UNSURE**; the alternative is to
    write only the plain node and let downstream consumers' variant inputs
    become ROOTS, which is exactly the I3 cut re-created by the cap. Writing
    every consumed variant is the "bundle = all variants" reading of the cap
    and keeps I3 at zero; its cost is that the one copy drives every variant
    equally, which is the resolution the cap gives up anyway.
  I read the spec's I4 as the bare-set case; this bullet is the extension it
  needs for complexes, and it should be stated in the spec one way or the
  other.

The three cap registers (§0.11) collapse to one: `capped_complexes` (by
complex stId, from `variants(C)`) and `capped_reactions` (by reaction stId,
from `copies(R)`), both pure functions, both exported in `BUILD`-level logs.

### 2.6 Copy identity and provenance

Copy (VR) node ids stay uuid4 (the solver is label-free since specs/013). Each
copy records σ as `slot=member` pairs in `node_reaction_context.csv` (or a new
`copies.csv`), so a trace can say which isoforms a copy is for. The
`decomposed_uid_mapping` content hashes are not needed for emission; keep the
table for the cache/provenance contract or retire it — **UNSURE** which;
retiring it is cleaner but touches the cache fingerprint (§5.1).

---

## 3. Emission wiring

### 3.1 Registry and phases, keyed by VK

The registry key stays `(node_id, vr_uid, role)` with `node_id = VK`.
`_register_phase1` shares per `(VK, reaction stId, role)`: copies of one
reaction that chose the same member for a participant share its node; copies
that differ do not. Phase 2 merges `outputs(p_vr) ∩ inputs(f_vr)` by VK over
`reaction_connections` pairs — unchanged code, new ids. Root and terminal
caches share per VK across the pathway.

### 3.2 I2 — what "the same node" can mean

Within what the generator does today, I2 is: *for every `(p, f)` in
`reaction_connections`, every VK in `outputs(p) ∩ inputs(f)` resolves to one
uuid*. That is a consequence of 3.1 once producer and consumer compute the
same key, which §1 guarantees. The design does NOT add global unification
(every producer of VK with every consumer of VK regardless of linkage),
because that is the positional-uuid decision of `docs/UUID_DESIGN.md`, which
variants did not create. **UNSURE** whether Adam intends global; if he does,
it is a separate, measurable switch (it is the "silo bridge" family that
failed three times, specs/019). The I2 test should report both numbers: linked
pairs with mismatched uuids (must be 0) and unlinked same-VK producer/consumer
pairs (a census figure, not a failure).

### 3.3 A consumer needs a variant no reaction produces

The VK is a root: one root node per VK (root cache), decomposed by the
boundary layer by its STRUCTURE (§4.1): fixed components as today (existing
upstream node, else nested build, else leaves), and each set slot resolved to
its CHOSEN member — an existing upstream copy of that member if one is not
downstream of the root, else a leaf/nested node. No F7 set node is built for a
resolved slot (the set was decided by the key). Example: a consumer copy
needing C[S=B] while the pathway produces only C[S=A] gets C[S=B] as a root
fed by P (fixed) and B — which is the truth: B-containing C enters from
outside.

### 3.4 A produced variant with no consumer

Terminal node; dissociation sinks to its leaves (fixed leaves ∪ chosen
members' leaves), as today for a complex. A member complex chosen in a slot
is dissociated to ITS leaves (not left as a sink node), matching
`get_terminal_components` behaviour for plain complexes.

### 3.5 Catalysts and regulators attach to the copy's own node first

`append_regulators` resolves each (copy, catalyst) to VK(c, σ_copy) and
takes, in order: the copy's own input node for that VK
(`registry[(VK, vr_uid, "input")]`, which is the catalyst = input case,
§0.6), else the first existing node for VK (today's rule, §0.10), else a
fresh node (a root, decomposed by 3.3). Depletion edges
(`_emit_substrate_depletion_edges`) must compare PARENT stIds, not VK strings,
in their "a catalyst does not deplete itself" guards, and
`_pool_contains` must see the new pool nodes (2.5).

---

## 4. Interaction with existing layers

### 4.1 `boundary_edges.csv` (specs/044, F7) — I5

Kinds unchanged: `assembly`, `dissociation`, `set_member`. Required change:
`_is_complex`/`_is_set` must recognise a VK id (parse parent and slots), and
`_decompose_hier` must, for a variant root, join the chosen member for each
slot instead of building an F7 set node (3.3). F7 set nodes remain for plain
roots (a root complex with sets that is NOT a variant because its consumer
is capped) and for the capped-complex case. `_existing_upstream` must look up
nested set-bearing components by their VK (the produced VARIANT copy of a
nested complex), because the plain stId node of such a complex no longer
exists. `containment_structure.csv` is unchanged (Neo4j structure).

### 4.2 Set pools (specs/033) — replaced where expansion applies, kept past the cap

For uncapped reactions, a set catalyst/positive regulator becomes copies
(2.3), so `LNG_SET_POOL` has no effect there; past the cap the pool is the
wiring (2.5), and for negative-regulator sets a pool and direct OR edges are
equivalent (keep the pool for node economy). The flag stays meaningful as the
comparison arm (`LNG_VARIANT_NODES=0` must reproduce today's catalog
byte-for-byte modulo uuids). `DS_SET_POOL_MODE` applies only to pools that
remain.

### 4.3 Variant node sharing (specs/020)

Mechanism unchanged (`_register_phase1`), now keyed by VK: it shares the SAME
variant across copies and keeps different variants apart, which is what its
docstring wanted and could not have under plain-stId ids. The "no pathway
gains cyclic nodes" property should be re-measured: splitting a node into
variants can only remove cycles that went through the merged node, but the
new catalyst copies add catalyst→reaction edges per copy, and the specs/018
weld check (assembly leaves reusing downstream nodes) runs on more roots.

### 4.4 The July `LNG_COMPLEX_AS_NODE` silo

Dissolved by construction: the produced node and the catalyst node are the
same VK. The flag's remaining meaning is the capped-complex rule (VK = stId),
which is now part of the key function, so the flag can expire after the arm
under the specs/020 policy. The dead `LNG_COMPLEX_AS_NODE=0` paths in
`_map_annotated_entity_to_nodes` and `_expand_complex_variants` (§0.3) are
replaced by the one key function and deleted.

### 4.5 Matching layer (`reaction_generator.py`, `best_reaction_match.py`)

`break_apart_entity`/`get_broken_apart_ids` enumerate per-participant
products and the matcher pairs them (§0.4). Under 2.1–2.2 the enumeration is
per-REACTION over V(R) and outputs are derived, so the matcher is not called
for variant reactions. Minimal transitional fix if the matcher is kept for a
while: sort `inputs`/`outputs` by uid before `linear_sum_assignment` (removes
the seed dependence, F9's determinism half) and add the same-set-stId overlap
as a first-tier cost (F9's semantic half). Cleavage reactions (one input,
several outputs, zero overlap) are not variant reactions and keep today's
"every output is produced" behaviour.

### 4.6 pools.csv (specs/039), drugs.csv (specs/032), node_resolution.csv, nodes.csv

Every tail parser (§0.9) is replaced by one `parse_variant_key(id) ->
(parent, [(slot, member)])` and one `variant_leaves(id)` = terminal leaves of
the parent's fixed components ∪ terminal leaves of each chosen member.
`export_nodes` writes `source_sets` = the slots and `chosen_members` = the
members straight from the key (no Neo4j re-derivation by leaf intersection,
which is what `_derive_sets_and_chosen` does and which is ambiguous for the
monomer/dimer set). `find_pools.proteins` uses `variant_leaves`. The
`member_leaves` column is what the benchmark resolves genes through, so T12
(§5.2) checks it on a fixture.

### 4.7 The solver (deltasignal)

Reads the parent via `split("::variant::")[0]` only; the new token format
(`slot=member`) is inside the tail and invisible to it. `containment.csv` must
have a row per VK id (keyed by the VK string, as `export_containment` keys by
`reactome_id_to_uuid` values today) or the self-inhibitor rules
(specs/022/040) and the `DS_CONSERVED_MODE`/`DS_DRUG_MODE` lookups go inert
on variant nodes. Check before the arm: the solve reports
`self_inhibitor_leaf_pairs` and `drugs_held`; they must not fall to zero in
pathways that had them.

---

## 5. Failure modes, tests, measurement

### 5.1 Failure modes

1. **Key collision** (two σ, one key): guarded by slot labels and member ids
   (§1.2); T7.
2. **Key mismatch across sides** (the F10 class): one function, one place;
   any second encoder (e.g. a leaf-based fallback in a boundary helper) is a
   regression. T1/T2 with nested-set fixtures.
3. **Cap asymmetry**: a complex capped on one side only (impossible if
   `variants(C)` is per-stId); a reaction-level cap re-creating I3 cuts
   downstream (2.5, design α); T4 includes a capped producer with an
   uncapped consumer.
4. **Wrong output binding** (2.2 ties/partials): counted, logged, and the
   RAF/FGFR2 traces check specific reactions by name.
5. **Explosion**: catalyst and regulator slots multiply copies; shared slots
   (2.1) cut the growth, the cap bounds it. The census (5.5) must be run
   BEFORE the code, as the spec orders.
6. **Silent root variants**: `_is_complex` returning False for variant ids
   (§0.9) means I5 passes while every variant root loses its assembly edges;
   T5 fails on that.
7. **Benchmark coverage drop masquerading as accuracy**: `member_leaves`
   wrong for the new tail format; T12 and a catalog-level check that every
   gene in `bench/catalog_pathways.tsv`'s case set resolves to at least as
   many nodes as before.
8. **Cache reuse**: `_cache_is_reusable` adopts a pre-fingerprint cache
   (generator CLAUDE.md); the new enumeration must be in the fingerprint, and
   the first regeneration must be done with `rm -rf <out>/*/cache`.
9. **Depletion self-guard by string** (3.5): a variant catalyst that is the
   substrate's parent would deplete it unless guards compare parents.
10. **Edge dedup** (`create_pathway_logic_network` tail): the dedup key omits
    stoichiometry; copies that differ only by a slot that maps to the same
    node (stoichiometry-2 slots) are intended duplicates and must collapse —
    they do, since the keys are equal.
11. **Phase-1 "unmapped" copies**: a copy without a reaction stId is
    registered per-variant; with more copies per reaction the inconsistency
    grows. Make it an error, not a warning.

### 5.2 Test fixtures (pure-Python, monkeypatched connector, like
`tests/test_logic_network_generator.py`)

- **T1 (I1)** C = {P, S={A,B}}, R: C + X → CX with CX = {P, S, X}. Expect two
  copies; copy A has inputs {C[S=A], X} and output CX[S=A]; no copy's input
  or output set contains both A and B. Bare-set form: S + X → S:X. Nested
  (F10) form: input set T = {D1, D2}, D1 = {Q, S}: copies (D1,A), (D1,B),
  (D2); the input node is `D1::variant::S=A`, never {Q, A}.
- **T2 (I2)** R1: Y → C, R2: C + X → CX, linked. For each member m:
  uuid(R1 output C[S=m]) == uuid(R2 input C[S=m]). Unlinked control: distinct
  uuids and the "unlinked same-VK" counter = 1.
- **T3 (I3)** R1 produces C[S=A] only (its set S' = {A}); R2 consumes C with
  S = {A, B}. Expect C[S=B] a root with assembly edges from P and B; expect NO
  `input`/`and` edge from a node for P into any copy of R2. Catalog-level:
  `excess_classify3.py` hand-off-cut = 0 and same-set = 0 (I1's 156).
- **T4 (I4)** `MAX_VARIANTS = 2`: (a) C with 3 variants → VK = C, one node on
  both sides; (b) reaction with S (3 members) × T (3 members) = 9 > 2 → one
  copy, one pool per bare set fed by 3 `set_member` edges, reaction reads each
  pool once; (c) capped producer R1 of C1 (uncapped, 3 variants) feeding an
  uncapped consumer R2: R1's copy writes output edges to the C1 variants R2
  consumes (design α) — no copy of R2 requires a root P.
- **T5 (I5)** `boundary_edges.csv` edge types ⊆ {assembly, dissociation,
  set_member} before and after; a variant ROOT gets assembly edges (fails
  today's `_is_complex`); a plain root with a set still gets an F7 set node;
  a terminal variant gets dissociation sinks for fixed ∪ chosen leaves.
- **T6 (determinism)** Build a fixture pathway under `PYTHONHASHSEED=0` and
  `=1` (subprocess) and with `get_set_members` returning members in reversed
  order; the multiset of (source id, target id or reaction stId, edge_type,
  and_or, stoichiometry) must be identical, with VK strings as ids.
- **T7 (key canonicity)** VK(C, σ) identical for σ given in any order;
  S1={A,B}, S2={A,B}: (A,B) ≠ (B,A); monomer/dimer set: different keys;
  capped C: VK = C regardless of σ; modifier set: atomic; set-of-sets: the
  outer slot label, the flattened member.
- **T8 (catalyst == input)** R: C → C', catalyst C: each copy's catalyst node
  uuid == that copy's input node uuid; copies = |S|, not |S|².
- **T9 (output binding)** u = {A2 {A:2}, B2 {B:2}, AB {A:1,B:1}}, w = {pA2,
  pB2, pAB} with equal multisets → bound 1:1; tie fixture (two w-members with
  multiset {A:1}) → paired by sorted stId and `slot_binding_ties == 1`;
  size-mismatch fixture → `slot_binding_partial == 1` and the unmatched copy
  produces no w-variant.
- **T10 (stoichiometry > 1)** C = {2 × S}: two copies (A,A) and (B,B), input
  edge stoichiometry 2, no (A,B).
- **T11 (negative regulator)** inhibitor N = {S, Z={z1,z2}} on R with input
  C ⊃ S: copy A is inhibited by N[S=A,Z=z1] and N[S=A,Z=z2] only.
- **T12 (consumers)** `export_nodes` on a variant id gives `member_leaves` =
  fixed leaves ∪ chosen members' leaves, `source_sets` = slots,
  `chosen_members` = members; `find_pools.parts`, `_node_leaves`,
  `export_drugs`, `export_node_resolution` agree with `parse_variant_key`.
- **T13 (flag off)** `LNG_VARIANT_NODES=0` reproduces today's edge multiset
  (by stId strings) on the fixture suite.

### 5.3 Pre-registered directions for the arm (both axes, held-out split)

- Pathways whose reactions have set CATALYSTS and are under the cap (RAF's
  R-HSA-5674496/5675194/5675198 via "p-T,Y MAPK monomers and dimers"; FGFR2
  autophosphorylations via the FGFR2 dimer sets): a single-member knockout
  (one MAPK isoform, one FGF) now stops only its copies. Expect FEWER
  predicted changes for those cases; held-out effect sign unknown — specs/033
  measured `mean`-like pooling as neutral and `max` as harmful, and this is
  between them. Pre-register: curator cases perturbing MAPK1 or MAPK3 alone
  in RAF move toward "no change"; cases perturbing a gene present in every
  member do not move.
- FGFR2, FGFR1/3/4, R-HSA-194315, Mitotic G1, DDX58, EPH-ephrin (the I3
  worst list): expect MORE predicted change reaching readouts downstream of
  the formerly cut hand-offs. Pre-register: the 2,504 cut cases' readouts
  gain non-baseline predictions; direction correctness is the measurement.
- RAF itself: mostly unchanged, because its big complexes are capped either
  way (§0.7) — if RAF moves a lot, something other than this design moved it.
- Node count: up where sets are inside produced complexes, down where the
  dissolved leaves disappear; edge count: catalyst copies up. Report both
  and the cyclic-node count per pathway.

### 5.4 Measurement plan

1. Size census (spec step 1), with the shared-slot rule applied, reporting
   per pathway: copies today, copies under 2.1, copies with independent
   catalysts (to show what 2.1 saves), reactions over cap, participants with
   one set at two paths (1.3), partial bindings (2.2.3), and stoichiometry>1
   slots. (`size_census.py` exists and computes the independent-product
   figure; it should be extended rather than duplicated.)
2. `excess_classify3.py` extended to split "root subunit" into "no produced
   input from the holder" (total cut) and "another leaf conducts" (§0.8),
   so the before/after is honest about what was disconnected.
3. Tests T1–T13 green; catalog built with `scripts/catalog.sh build` against
   a generator worktree at the flag-on commit; `BUILD.json` records the flag.
4. Arm against the canonical build (`0b1d535` results): held-out and
   experimental, McNemar, per-pathway net, distinct readouts moved, cyclic
   nodes per pathway, convergence, and the label-invariance re-check
   (specs/013) since variant ids change the tie structure in `flow` order.
5. Traces: RAF (R-HSA-5672972/5672978 — confirm capped behaviour unchanged),
   FGFR2 (R-HSA-190408 → 5654399 → 5654397 — confirm one node per variant
   across the three), PIP3, against MP-BioPath's networks to the first
   divergence.

### 5.5 Census numbers from this derivation

Method: the same one-hop bulk load `size_census.py` uses (every
PhysicalEntity with its direct children; no variable-length set paths), then
per reaction: `today` = product over input participants; `independent` =
that × catalysts × positive regulators; `shared` = product over DISTINCT set
stIds across inputs, catalysts and positive regulators (rule 2.1). Each
figure is `min(product, 512)` summed over reactions, as `size_census.py`
counts. Caveat: for a set whose member complexes carry DIFFERENT inner sets
the script multiplies where it should add, so `shared` and `independent` are
upper bounds there (RAF's and FGFR2's member complexes share their inner
sets, so those two are exact or nearly so).

RAF R-HSA-5673001, FGFR2 R-HSA-5654738, PIP3 R-HSA-1257604 together
(210 reactions):

| figure | value |
|---|---|
| reactions with any set slot | 164 |
| reactions with catalyst / positive-regulator slots | 67 |
| … of which the catalyst slot is also an input slot (§0.6) | 18 |
| reactions needing output-slot binding (2.2) | 34 (35 of 345 output slots) |
| capped-copies today (inputs) | 12,225 |
| capped-copies, independent product | 19,105 |
| capped-copies, shared-slot rule | 16,679 |
| reactions over the cap: today / independent / shared | 8 / 22 / 14 |

Per pathway (shared ÷ today; copies today → shared; over-cap today → shared):
PIP3 1.75× (855 → 1,493; 0 → 1), FGFR2 1.37× (6,358 → 8,712; 0 → 3),
RAF 1.29× (5,012 → 6,474; 8 → 10). So in these three the expansion costs a
third more copies and pushes six more reactions over the cap; the shared-slot
rule saves 2,426 copies and eight cap crossings against independent
expansion.

Two of the §1.3/§2.1 **UNSURE** structural cases were counted directly
(bounded `hasComponent*1..5` from `DISTINCT` participant complexes): a set
reachable by two different paths inside one participant complex — **0** in
the three pathways; the same set inside two different INPUT participants of
one reaction — **0**. The simple slot label (set stId) therefore suffices
there; the catalog-wide count is in the appended run below if it finished.
A first query that omitted `DISTINCT c` reported 10 for the first case; that
was the number of reactions using the complex, not paths — recorded so no
one repeats it.

Output slots absent from the input/catalyst side (the binding case, by a
separate bounded query): RAF 19 of 53 reactions with output sets, FGFR2 5 of
38, PIP3 13 of 46.

**Catalog-wide** (`bench/catalog_pathways.tsv`, same script, same caveat;
the two zero-counts above were NOT re-run catalog-wide):

| figure | value |
|---|---|
| reactions | 4,975 |
| reactions with any set slot | 2,799 |
| reactions with catalyst / positive-regulator slots | 1,050 |
| … of which the catalyst slot is also an input slot | 449 (43%) |
| reactions needing output-slot binding | 564 (634 of 4,542 output slots, 14%) |
| capped-copies today / independent / shared-slot | 113,145 / 204,811 / 144,213 |
| reactions over the cap: today / independent / shared | 111 / 259 / 155 |

So full expansion under the shared-slot rule is +27% capped-copies and +44
reactions over the cap catalog-wide, against +81% and +148 for independent
expansion; the rule matters most where catalysts are the substrate complex
(449 reactions). Largest growth factors (shared ÷ today): R-HSA-2871796
16.6× (66 → 1,094; two reactions newly capped), R-HSA-2404192 9.4× (62 → 585),
R-HSA-74752 3.1× (438 → 1,341), R-HSA-453279 Mitotic G1 2.3× (787 → 1,832),
R-HSA-5654743 FGFR4 1.9× (683 → 1,318, but 4,332 independent). Those are
where to look first when the arm moves. The `size_census.py` figure
(independent product) is the one this replaces; both scripts should be
reconciled before the spec's step 1 is called done.

---

## 6. Open questions (for Adam)

1. **I2 scope**: linked pairs only (positional uuids kept) or global by key?
   I derived linked-only (3.2) and flagged global as a separate switch.
2. **Cap on a produced, uncapped complex inside a capped reaction** (2.5):
   write every consumed variant (my choice) or the plain node only?
3. **Partial output binding** (2.2.3): produce nothing for the unmatched
   copy (my choice) or best-overlap as today?
4. **Negative-regulator variants** (2.4): shared slots bound, unshared
   enumerated — and does the solver's inhibitor stacking under `divide` need
   a cap on k?
5. **Same set in two different input complexes of one reaction**: one
   variable (my rule) — the census should say whether the case exists.
6. **Stoichiometry > 1 heteromers**: not now; the census says what it costs.
7. **Retire `decomposed_uid_mapping`** or keep it for provenance/cache?
8. **Token format** `slot=member` inside the existing `::variant::` tail:
   acceptable to every consumer (one parser, §4.6)? The solver is untouched.
