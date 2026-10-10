# specs/046 derivation A (Opus): variant key, copies, emission

A blind derivation against `spec.md` (Adam's design, invariants I1-I5).
- Read-only on code. Generator: `logic-network-generator` main `6990015`.
- Neo4j: Release97. Build inspected: `20261004-1959_6990015`.
- Scratch scripts: `bind.py`, `stoich.py`, `tree.py` in the session
  scratchpad. They are not committed; the numbers below say what each one
  measured.

Every point I am unsure of is marked **[UNSURE]** or **[DECIDE]** (a choice
for Adam rather than a fact).

---

## 0. What the code does today

Node identity is chosen in three places, and each uses a different naming
rule.

1. **Matching layer** (`reaction_generator.py`).
   - `break_apart_entity` returns three things: a set's members (flattened
     through nested sets), the combinations of a complex that contains sets,
     and the stId of anything else.
   - A combination uid is `sha256(sorted(component ids))`. It is
     **content-addressed and not keyed by the complex stId**, so two
     complexes with the same chosen leaves collide. That is why emission
     stopped trusting it.
   - `complex_components` is a dict `member -> stoichiometry`. A set at
     stoichiometry 2 is therefore ONE choice: both copies get the same
     member (a homodimer).
   - Input combinations are paired with output combinations by
     `find_best_reaction_match`. That is Hungarian assignment on
     `|components ∩|`, where a component is mapped to its **HGNC gene**
     reference (`get_component_id_or_reference_entity_id`). Ties go to
     whichever assignment `linear_sum_assignment` returns, and surplus rows
     go to `argmax`, which takes the first index (F9).
2. **Emission** (`_resolve_vr_entities` → `_map_annotated_entity_to_nodes`).
   It works from the reaction's annotated entities:
   - a complex becomes its plain stId (`LNG_COMPLEX_AS_NODE=1`);
   - a bare set becomes `member_set ∩ _matching_leaves(set)`.

   `_matching_leaves` descends into any complex that contains a set. So a set
   whose member is such a complex is **dissolved into that complex's
   components** (F10).
3. **Catalysts and regulators** (`append_regulators`).
   - A complex becomes its plain stId.
   - A bare set becomes one `LNG_SET_POOL` node per pathway, fed by its
     members (specs/033).
   - The node is found by `stid_to_existing_uuid`, which takes the **first
     registry entry carrying that stId**. That is a global-by-stId join, not
     a positional one. It is how the July catalyst silo was closed.

**Traced in the current build (FGFR2, R-HSA-5654738):**
- "Activated FGFR2 binds SHC1" (R-HSA-5654404) has 68 copies. Each copy's
  required inputs are four separate free species: `HS`, one FGF, one
  p-FGFR2 isoform, and one SHC1 isoform.
- "Activated FGFR2" (R-HSA-5654152) is a DefinedSet of DefinedSets of
  complexes, and each of those complexes contains the set `FGFR2b-binding
  FGFs ×2`. Under F10 that dimer has become separate molecules.
- The same `HS` node (R-ALL-190915) is produced by R-HSA-190408 (FGFR2b),
  R-HSA-190413 (FGFR2c) and R-HSA-5654157, and it feeds every copy. A FGFR2b
  ligand therefore feeds FGFR2c copies through `HS`.
- The SHC1 isoforms (R-HSA-1169451/1169488) are roots, which is correct:
  nothing produces them.
- In the autophosphorylation step R-HSA-190408, the input is dissolved but
  the catalyst, the SAME entity R-HSA-192615, is a set-pool node. The two
  roles of one molecule are two unrelated nodes.

What the design must replace:
- **Key:** one naming function instead of three rules.
- **Pairing:** the HGNC-level Hungarian pairing.
- **Catalyst and regulator emission.**

---

## 1. The variant key

### 1.1 Variant trees

For every entity `X`, `V(X)` is the ordered list of its variants. A variant
is a **choice tree**, built bottom-up from Neo4j structure only:

| `X` is | `V(X)` |
|---|---|
| simple (EWAS, SimpleEntity, Polymer, OtherEntity, Drug, GEE, Cell), a complex containing no set, a modifier-isoform set (`modifier_isoform_set_ids()`), or a set with no members | `[Leaf(X)]` |
| a set (`DefinedSet` or `CandidateSet`; hasMember **and** hasCandidate both count as alternatives, as `get_set_members` does today) | the concatenation of `V(m)` over its members `m`, de-duplicated by key. **The set is transparent:** a variant of a set IS a variant of one of its members. |
| a complex containing a set, with components `(c_i, k_i)` | the product over components of each slot's choice: a component with no set is fixed; a varying component chooses from `V(c_i)` (stoichiometry rule in §1.4). Each variant is `Cx(X, {c_i: [chosen subtrees]})`. |

**Entity cap.** If the raw count `|V(X)|` (computed from the children's
already-capped counts) exceeds `LNG_MAX_VARIANTS`, `X` is **entity-capped**:
`V(X) = [Leaf(X)]`, an opaque generic node keyed by its plain stId.
- The cap is applied at every level, bottom-up, so whether `X` is capped
  depends on `X` alone. It never depends on the reaction or on whether `X` is
  nested.
- This is the property the key needs to be context-free.
- It replaces the reaction-scoped `CAPPED_IDS` for naming. A reaction-level
  cap still exists for copy counts (§2.4), but it never changes a node's key.

### 1.2 The key function

```
vkey(Leaf(X))                = X
vkey(Cx(X, slots))           = X + "::variant::{" +
                               ";".join(f"{c}=" + "+".join(sorted(vkey(t) for t in chosen[c]))
                                        for c in sorted(varying slots)) + "}"
```

**Sets do not appear in a key. A set slot inside a complex does.** The slot
name is the component stId: the set stId, or the stId of the nested complex
that contains sets.

Examples (FGF3 = R-HSA-189886, SHC1-3 = R-HSA-1169451):

```
R-HSA-190411::variant::{R-HSA-189967=R-HSA-189886}
    Activated FGFR2b long homodimer with FGF3 (the slot 189967 x2 holds FGF3 x2 under the homo rule)
R-HSA-5654279::variant::{R-HSA-1169480=R-HSA-1169451;R-HSA-5654152=R-HSA-190411::variant::{R-HSA-189967=R-HSA-189886}}
    Activated FGFR2:SHC1, with the FGFR2 slot holding the variant above
```

Why each rule is there:
- **Transparency** is what makes producer and consumer agree.
  - R-HSA-190408 outputs the set R-HSA-192606.
  - R-HSA-5654404 reads its superset R-HSA-5654152 (a set of sets).
  - Both resolve to the member complex variant
    `R-HSA-190411::variant::{...}`, so the names coincide without any
    matching heuristic.
- **Keeping the slot name inside a complex** keeps the key unambiguous when
  two different sets of one complex share members. R-HSA-5672718 has both
  `activated RAF/KSR1` and `'activator' RAFs`, which share ARAF, BRAF and
  RAF1 (as different modified forms, but the same principle applies).
- **The outer stId is kept**, so the existing consumer
  `split(id, "::variant::")[1]` in `reaction_model.jl` still returns the
  parent stId. The parsers that read the tail as `_`-joined leaves have to
  change; they are listed in §4.6.

**Determinism.**
- The key is a pure function of: Neo4j structure (stIds, the three
  relations, stoichiometry), the release's modifier-set list, and
  `LNG_MAX_VARIANTS`.
- It involves no uuid, no hash seed and no traversal order: every list is
  sorted by key string before use.
- Variant enumeration order is "sorted by key", so the copy order is
  reproducible too.
- **[Optional]** Mint VR uids as `uuid5(ns, reaction stId + sorted binding
  keys)` and node uuids as `uuid5(ns, key + position)`. Rebuilds would then
  be byte-identical, which today's uuid4 minting prevents. This is not
  required by I1-I5.

**Side table.** A `variants.csv` row per key records:
- `key`, `outer_stid`;
- `chosen` (slot → member keys);
- `leaf_multiset` (the chosen leaves with cumulative stoichiometry);
- `entity_capped` (bool).

Every exporter that today string-parses `::variant::` reads this table
instead.

### 1.3 Nested sets, sets of complexes containing sets, CandidateSets, modifier sets

- **Nested sets** (max depth 5 at R97, no cycles, per `set_resolution.py`).
  Flattening through transparency is exact. A leaf reachable by two routes
  (S → S1 → x and S → S2 → x) is one alternative after de-duplication by key.
- **A set whose member is a complex containing sets** (the F10 shape) gets
  that member's own variants as alternatives. It is never dissolved. F10
  disappears by construction, because `_matching_leaves` is no longer used
  for identity.
- **CandidateSet:** members and candidates are both alternatives. This
  matches today. **[DECIDE]** Excluding candidates would shrink RAF scaffolds
  from 8 to 1 and is not what the curators' "any of" means for RAF.
  Recommend keeping both and recording the relation in `variants.csv`.
- **Modifier-isoform sets stay atomic** (Ub, SUMO, NEDD8, ATG8 …). They are
  leaves on both sides, so they never enter a key choice. This is unchanged
  and safe: the list is derived per release.
- **Drug members** (a set of RAF inhibitors, say) produce variant keys whose
  leaves include a drug. `drugs.csv` must flag them (§4.6).

### 1.4 Stoichiometry > 1 of one set inside one complex

Measured at R97:
- 605 `hasComponent` edges with stoichiometry > 1 point at a set, across
  498 complexes;
- 263 of those complexes are reaction participants, in 416 reactions.

Examples:
- FGF2:FGFR dimers (`FGFR2b-binding FGFs ×2`): biologically homodimeric
  ligand pairs.
- Histone H2B ×2 (14 members) and collagen α-chains ×3 (12 members).
- **R-HSA-5674136 "activated RAF homo/heterodimer"** = `activated RAF
  monomer ×2`. Here the curator explicitly means heterodimers too.

Two rules, measured on the catalog's reactions with my `stoich.py`
(modifier sets atomic, cap 512):

| rule | copies, inputs only | copies, full expansion | reactions over cap (full) | entities over cap | distinct uncapped variants |
|---|---|---|---|---|---|
| **homo**: one choice per slot, stoichiometry carried on the edge | 84,763 | 178,334 | 222 | 94 | 44,780 |
| **multiset**: k independent choices, unordered | 124,103 | 235,053 | 321 | 142 | 52,702 |

Further sizes:
- `Activated FGFR2` (R-HSA-5654152): 34 variants under homo, **186** under
  multiset. Multiset invents FGF heterodimer pairs nobody curated.
- RAF dimer R-HSA-5674136: 3 under homo, 6 under multiset (homo loses
  BRAF:RAF1, the dominant RAF dimer).

**Recommendation: `homo` as the default.**
- It is I1 read literally ("exactly one member of every set").
- It is Adam's "n × m" counting.
- It is today's matching-layer behaviour, so it changes nothing silently.
- It is smaller.

`multiset` would be an arm (`LNG_SET_STOICH=homo|multiset`). Under
multiset, I1 is read per stoichiometric unit.

**[DECIDE]** Neither rule is right everywhere. FGF wants homo; RAF dimer and
GABA-A receptors want multiset. Reactome has no field that distinguishes
them. The only alternative to a global rule is a per-set curated list, which
I do not recommend.

The same rule applies to a bare set that is a reaction participant with
stoichiometry > 1, because the reaction behaves like a complex of its inputs.

---

## 2. Copies (virtual reactions)

### 2.1 Participants and coupling

For a reaction `R`:
- inputs `I`, outputs `O`, catalysts `K`, positive regulators `P`, negative
  regulators `N`;
- all from `get_reaction_input_output_ids` and the three Cypher maps used
  today.

The participants form **choice groups**:
- **One stId in several roles is one group.** Input == catalyst occurs in
  309 catalog reactions where the entity contains a set. Examples:
  autophosphorylation R-HSA-190408 (FGFR2b), R-HSA-5674373 (MAP2Ks
  phosphorylate MAPK, input == catalyst R-HSA-5674360), R-HSA-1963581
  (trans-autophosphorylation of ERBB2 heterodimers).
  - It is one molecule, so a copy uses one variant in both roles.
  - Independent expansion would create "FGFR2b:FGF3 phosphorylated, catalysed
    by FGFR2c:FGF10" copies.
  - With divide-form inhibitors it would also square the cross-terms.
- **Different participants that share a nested set** are independent
  groups. This happens 371 times among catalyst/regulator participants; for
  example, RAF binding's input and its negative regulator R-HSA-5675413 both
  contain `mature p21 RAS`. They are different molecules.
  **[DECIDE]** Coupling them would shrink copies but asserts the two RAS
  molecules are the same isoform. I recommend against it.

Copies **below the cap** are the product over choice groups of `V(group)`.
This is Adam's "one copy per combination over all participants", catalysts
and regulators included. It is enumerated in sorted-key order.

### 2.2 How a copy's choices reach its outputs (replacing F9)

Each output `O` has its own **varying occurrences**: every set or
set-containing complex in `O`'s structure, and `O` itself if it varies. Each
occurrence must get a choice. Within one copy, the choices already made are
the chosen trees of the input group variants. They are visible occurrence by
occurrence, each with its **ancestor chain** (the stIds of the enclosing
entities up to the participant).

Bind each output occurrence `o` (stId `s`, chain `A_o`) by the first rule
that applies:

1. **Identity.** Some input occurrence in this copy has the same stId `s`.
   The output takes that occurrence's chosen subtree, the whole subtree when
   `s` is a complex.
   - If several input occurrences carry `s`, keep those whose ancestor chain
     shares the longest common prefix (read from the occurrence upward) with
     `A_o`.
   - If they still disagree on the chosen member, **fan out**: one output
     variant per distinct choice. Never choose by position.
   - Under transparency, an output set whose stId equals an input set's
     binds even when they sit at different depths.
2. **Correspondence by reference entity.** Otherwise, for every input
   occurrence `t` in this copy:
   - compute the overlap between the reference entities of the input's
     CHOSEN member and those of each output member;
   - use the **isoform-specific** reference (`variantIdentifier` when
     present, else `identifier`, else ChEBI or the stId);
   - exclude cofactors and modifier sets;
   - keep the output members with maximal overlap over all `t`. One means
     bound; ties mean fan out.

   The isoform level is not optional. In FGFR2b autophosphorylation the
   input set R-HSA-192615 is {long 190227, short 192591}. The output set
   R-HSA-192606 is {long 190411, short 192597}. At the HGNC level (today's
   matcher) every pairing ties: FGFR2, the FGFs and HS are shared by all
   four. At the isoform level, long P21802-3 matches long, and short
   P21802-18 matches short, uniquely.

   Inside the chosen output member, the nested slot `FGFR2b-binding FGFs`
   then binds by rule 1.
3. **Unbound.** Nothing in the copy corresponds, so fan out over every
   alternative of `o`. The curator says the output is "any of these" and
   nothing says which. This matches what today's surplus pairing achieves.

**Where each rule applies.** Every output set occurrence in the catalog's
reactions, measured (`bind.py`):

| | occurrences |
|---|---|
| all | 4,584 |
| identity (same set stId in an input) | 3,915 (85.4%) |
| reference entity, bijective | 383 (8.4%) |
| reference entity, covering (one input member → several outputs) | 194 (4.2%) |
| unbound | 92 (2.0%) |

- The 736 "ambiguous" identity hits were counted over the generic
  structure. Most are the same set under ALTERNATIVE members of an outer set,
  which never co-occur in one copy. Binding per copy, against chosen trees,
  removes those.
- **[UNSURE]** How many remain genuinely ambiguous after per-copy binding is
  not measured. Example to trace: R-HSA-1963581.

A copy with `f` fan-out choices becomes `f` copies with identical inputs and
different outputs. "Output variants of one copy" is never a set of
simultaneous outputs. Fanning out keeps the meaning "this reaction yields one
of these".

**Worked RAF example (R-HSA-5672972, "MAP2Ks and MAPKs bind to the activated
RAF complex").**
- Inputs: R-HSA-5672718 (80 variants: 4 × 4 RAS × 5), MAPKs (2), MAP2K
  homo/heterodimers (3), RAF/MAPK scaffolds (25: the Focal Adhesion
  candidate itself contains sets). 80 × 25 × 2 × 3 = 12,000, as in specs/045.
- Output: R-HSA-5672720, whose components are **those same four entities**.
  Every output slot binds by identity, so with 12,000 copies there would be
  exactly one output variant per copy.
- The negative regulator R-HSA-5675413 (RAF1:PEBP1, with `'activator' RAFs`
  and `mature p21 RAS` inside, 20 variants) is an independent group.
- Full expansion is therefore 12,000 × 20. R-HSA-5672720 itself has 12,000
  variants, more than 512, so it is **entity-capped**. Its key is the plain
  stId, and §2.4 applies.

### 2.3 Catalyst and regulator edges per copy

Each copy carries one edge per catalyst or regulator group:
- from the chosen variant node;
- `catalyst` / `regulator` edge types as today;
- pos/neg as today;
- `and` for positive and `or` for negative, as today.

The copies are OR'd where their outputs are consumed. Output edges get
`and_or = "or"` when more than one copy produces the key, as now.

- `LNG_SET_POOL` (specs/033) is **replaced below the cap**, because a set
  catalyst is now expanded.
- It survives only as a fallback past the cap (§2.4).
- The spec's pre-registered semantic change applies: a member knockout stops
  only that member's copies.

### 2.4 The cap, and how it composes

There are two caps:
- **Entity cap** (§1.1): changes the key, so it must be identical on every
  side.
- **Reaction cap:** never changes a key. It only changes how a reaction
  reads its participants.

**One observation decides the fallback.** Under DeltaSignal's defaults
(AND = product with `hill_sat`, OR across copies = `DS_OR_MODE=mean`), take a
participant group whose choice does not reach any output, such as a catalyst
or a positive regulator. Expanding it into copies is
`mean_c(A · x_c) = A · mean_c(x_c)`. That is **exactly** one per-copy OR pool
over its variant nodes, as long as the pool is combined with the SAME
operator as the copies. `max` distributes the same way.

It is not exact in two cases:
- at `hill_sat`'s 100-fold saturation, where the product is capped;
- for **negative regulators** under `divide`, where
  `mean_c(A·h(r_c)) ≠ A·h(mean_c r_c)`.

This matches Adam's specs/045 framing ("the same problem solved once per
member and recombined").

Fallback order when a reaction's copy count exceeds the cap:

1. **Pool the positive groups whose choice binds no output** (catalysts,
   positive regulators, non-binding inputs). Each becomes one pool node per
   group: one per catalyst/regulator stId per pathway, as specs/033 has now,
   with its members now being VARIANT NODES. The other groups are still
   expanded. This is near-exact.
2. **Factorize by output.** Index copies by the output-variant tuples. Each
   copy reads, for every input group, a pool over the group variants
   compatible with that output tuple, and a direct edge when only one is
   compatible.
   - The output count is at most the cap for each output that is not
     entity-capped.
   - Negative regulators stay expanded if that still fits. Otherwise they
     are pooled, with the inexactness above logged.
3. **Single copy.** If the output tuples alone exceed the cap, or every
   output is entity-capped, emit one copy that reads one OR pool per varying
   participant and produces every output variant. RAF R-HSA-5672972 lands
   here: its output is entity-capped. This is specs/045's shape, with pools
   over variant nodes instead of leaf nodes.

**[DECIDE] I4 wording.**
- "Past the cap, one OR pool per set" holds literally for a bare-set
  participant, where the participant pool IS the set pool.
- For a complex participant containing sets, the faithful pool is over its
  VARIANT nodes, one pool per participant, not one per inner set.
  - Per-inner-set pools would need a generic complex assembled from slot
    pools (the F7 boundary shape).
  - The producers make variant nodes, not that generic node, so I2 would
    break.
- I recommend rewording I4 to "one OR pool per varying participant, the
  pool over its variant nodes". Today a complex participant past the cap is a
  single plain-stId node, so this is no worse.

**[DECIDE] Pool operator.**
- The copy-equivalence needs the pool combined like the OR across copies
  (`DS_OR_MODE`, default `mean`).
- specs/033 pools are combined by `DS_SET_POOL_MODE=product`.
- Two options:
  - emit fallback pools with a distinct edge type (say `variant_member`), so
    the solver can give them the copy-OR operator;
  - accept that the fallback reproduces 033's product semantics past the
    cap.

  The first makes "past the cap" an approximation of the expansion rather
  than a different model. That is the property specs/045 was after.

**Option F, an equivalent emission below the cap.**
- Step 1 alone, applied below the cap too, would cut full expansion
  (178k copies under homo, `stoich.py`) back toward the inputs-only count
  (~85k).
- It would leave predictions unchanged under product/mean, except at the
  saturation cap.
- The spec asks for full expansion, so I list this as a size lever for the
  implementation arm, not as the design.

---

## 3. Emission wiring

### 3.1 Node identity

The registry is unchanged in shape: `(eid, vr, role)`. Now `eid` is the
`vkey`.
- **Phase 1** (`_register_phase1`): variant sharing (specs/020) keys on
  `(vkey, reaction stId, role)`. It works unchanged and becomes more
  meaningful: the copies of a catalyst-expanded reaction share their input
  variant nodes.
- **Phase 2**: for each connected pair (precedingEvent, plus diagram-drawn
  pairs merged by `augment_reaction_connections`), a producer copy's output
  and a consumer copy's input with the same `vkey` are unioned.
  - That is I2 for connected pairs, automatically, because both sides call
    the same `vkey`.
  - Join by a key index per connection, not the current all-pairs
    `p_vr × f_vr` loop: 512 × 512 per connection is too slow.
- **Catalysts and regulators:** add a third role, `modifier`, to Phase 2.
  Union a copy's catalyst/regulator variant node with connected producers'
  output nodes of the same key.
  - Where no connected producer exists, keep today's global fallback (the
    July silo fix), now keyed by `vkey`.
  - Make the fallback deterministic: prefer a produced copy, then a root,
    then a fresh node, ordered by registry insertion.
  - **[UNSURE]** Insertion order follows `pd.unique` over the
    `reaction_connections` query, and Cypher gives no order without
    `ORDER BY`. Add `ORDER BY` to `get_reaction_connections` or sort the ids.

**[DECIDE] I2's scope.**
- Read literally ("same key → same node"), I2 makes node identity global per
  pathway. That abandons positional identity, which specs/020 kept on
  purpose: bridging copies across unconnected reactions was measured harmful
  four times.
- Those failures were LEAF-to-functional-node bridges with fan-out of about
  112. A same-key whole-variant join is a different object, but it has not
  been measured.
- **Recommendation:**
  - make I2 hold for connected pairs and for catalysts/regulators (as the
    July fix does);
  - measure "global by key" as a separate arm;
  - report the count the arm would change: consumer input variant nodes
    whose key some unconnected reaction in the pathway produces.

### 3.2 A consumed variant that nothing produces (a root)

- It is registered once per key through the root-input cache, as today.
- The boundary layer decomposes it. The decomposition **must follow the
  variant's choices**:
  - a set slot contributes its CHOSEN member, not an F7 set node;
  - a nested varying complex contributes its chosen sub-variant, built or
    joined by key.
- Otherwise the FGF3 variant of a root FGF dimer is fed by a set node of all
  five FGFs, and an FGF1 knockout moves the FGF3 variant. That reintroduces
  the excess I1 forbids, one layer down.
- F7 set nodes remain where a set is genuinely unchosen: entity-capped roots,
  and the plain stId of a capped complex.
- `_is_complex` returns False for `::variant::` ids today, so variant roots
  would not be decomposed at all. Under this design it must consult
  `variants.csv`.
- `_existing_upstream` joins on `vkey`, and the specs/018 downstream rule is
  unchanged.

### 3.3 A produced variant that nothing consumes (a terminal)

- It gets dissociation sinks per variant, built from the variant's chosen
  `leaf_multiset`. They are fresh, as today.
- It is not pruned. It stays addressable as a readout.
- A copy whose outputs are all terminal keeps its edges.

### 3.4 What I3 does and does not guarantee

- `excess_classify3.py` counts only required `input` edges.
- Without dissolution, no copy has a subunit as a required input, so I3
  passes by construction.
- What I3 does not catch is a consumer reading the whole variant node as a
  ROOT because no **connected** reaction produces it. The boundary layer
  then assembles it from leaves, which is the spec's sentence ("never read a
  complex as separate root copies of its subunits") failing in
  `boundary_edges.csv` instead.
- Proposed **I3b**: no consumer input variant node is a root when a reaction
  connected to its consumer produces that key. This is the test of I2.
- Report the unconnected-producer count as a number, not a gate (§3.1).

---

## 4. Interaction with existing layers

1. **`LNG_COMPLEX_AS_NODE` (July, LNG #47).**
   - Superseded under the new flag (`LNG_VARIANT_NODES=1`). It existed only
     because catalysts were plain stIds and products were variants.
   - With one key function, the silo it fixed cannot arise for the same
     variant.
   - Make the two flags mutually exclusive at startup in `env_flags`.
     `_map_annotated_entity_to_nodes`, `_matching_leaves` (for identity),
     `_complex_variant_leafsets` and `_expand_complex_variants` all retire
     under the flag.
2. **Matching layer** (`decompose_by_reactions`, Hungarian).
   - No longer decides identity or pairing. §2.2 replaces it.
   - Keep writing `decomposed_uid_mapping.csv` for provenance only, or
     retire it. **[UNSURE]** whether anything downstream reads it.
3. **Set pools (specs/033).** Replaced below the cap, kept as the
   step-1 fallback past it (§2.4).
4. **Variant sharing (specs/020).** Unchanged mechanism, keyed by `vkey`.
5. **Boundary layer (specs/044, F7).**
   - Kinds are unchanged, which satisfies I5.
   - Content changes for variant roots: chosen members instead of set nodes
     (§3.2).
   - I5 should be restated as "kinds unchanged, and a variant root's
     boundary edges follow its choices".
6. **Exporters that string-parse `::variant::` tails as `_`-joined leaves**
   all move to `variants.csv`:
   - `_parse_variant_members`, `_node_leaves`;
   - `export_nodes` (`member_leaves` must be the CHOSEN leaves, because
     DeltaSignal pins by `member_leaves` under `root_cycle`);
   - `_derive_sets_and_chosen` (gets exact choices from the table);
   - `find_pools.parts` (specs/039 pools);
   - `export_containment` and `export_containment_structure` (keyed by
     stId; a variant key needs its own rows with its chosen leaves, or
     DeltaSignal's self-inhibitor rule sees every alternative as contained);
   - `export_cofactors`, `export_drugs`, `export_node_resolution`.

   **[UNSURE] DeltaSignal side.** A curator readout on a complex containing
   sets now resolves to many variant nodes. The benchmark must aggregate
   across them; `node_resolution.csv` needs a `variant` relation. It does
   not today, because each complex is one node.

---

## 5. Failure modes and tests

### 5.1 Failure modes

| # | failure | guard |
|---|---|---|
| F-a | Two code paths compute keys differently: the silo reborn. | One `vkey`. A test asserts that every node id in `stid_to_uuid_mapping` is either an annotated stId or a key present in `variants.csv`. |
| F-b | Nondeterministic choice order. | Sorted-key enumeration. Rebuild twice with different `PYTHONHASHSEED` values and require identical node-key multisets and edge-key multisets (uuid-free). |
| F-c | Positional tie-break survives somewhere. | No `argmax`/`linear_sum_assignment` in the flag path. Fan-out on ties is the rule. |
| F-d | Stoichiometry rule mismatched between the key and edge stoichiometry. | The `leaf_multiset` in `variants.csv` equals the cumulative stoichiometry `_resolve_to_terminal_reactome_ids` would give. |
| F-e | Entity cap computed per reaction. | The key of a capped entity is its plain stId in every pathway. |
| F-f | Size blow-up. | The cap. `stoich.py`: 178k copies under full expansion against 85k inputs-only (homo); `size_census.py`, which counts slightly differently: 221k against 113k. The mitigation is Option F (§2.4). |
| F-g | Coupled stId gets independent choices. | Test for input == catalyst. |
| F-h | Boundary decomposes a variant root through set nodes. | §3.2 test. |
| F-i | Readouts and pins aggregate wrongly over many variant nodes (DeltaSignal side). | `member_leaves` is the chosen leaves, and the readout aggregation is specified before the arm. |

### 5.2 Test fixtures

Synthetic fixtures need no database; they mock `get_labels`,
`get_complex_components` and `get_set_members`, as `set_resolution` does.
The real ones are marked `database`.

- **I1.**
  - Synthetic: complex `C = {A, S1={a1,a2}, D}` with `D = {S2={b1,b2,b3}}`
    nested. Expect 6 keys. Each copy's required inputs and its catalyst
    contain exactly one of `a*` and one of `b*`.
  - Same with `S1×2` under homo: still 2 variants, stoichiometry 2 on the
    leaf. Under multiset: 3.
  - Real: R-HSA-5654404 has 68 copies. Each reads ONE input node
    `5654152`-variant and ONE SHC1 node. No HS, FGF or FGFR2 leaf is a
    required input.
- **I2.**
  - Synthetic: R1 outputs set `S = {X, Y}` with `X = {s=(p,q)}`. R2 reads
    superset `T = {S, Z}`, precedingEvent R1 → R2. Every `X` key R2 reads
    is the uuid R1 produces.
  - Real (FGFR2): R-HSA-190408 → R-HSA-5654404. The node
    `R-HSA-190411::variant::{R-HSA-189967=R-HSA-189886}` has an `output`
    edge from a 190408 copy and an `input` edge into a 5654404 copy.
  - Catalyst variant: R-HSA-5654407's catalyst (R-HSA-5654279) is the node
    R-HSA-5654404 produces, per key.
- **I3.** `excess_classify3.py` on the arm build reports
  `PRODUCED, leaf node is a ROOT (hand-off cut) = 0`, `same-set = 0` and
  `dissolved* = 0`. Also report the I3b count.
- **I4.**
  - Synthetic: a reaction with an entity-capped output and three set
    inputs. Expect exactly one copy and one pool per varying participant.
    Every pool member is a variant node that also exists as a producer's
    output.
  - Real: RAF R-HSA-5672972 gives one copy, 4 pools (80/25/2/3 variant
    members), and the regulator expanded or pooled per §2.4.
- **I5.** `boundary_edges.csv` edge types ⊆ {assembly, dissociation,
  set_member}. A synthetic variant root `C{S=a1}` gets `a1 → C` (assembly),
  and no set node for `S`.
- **Binding (F9 replacement).**
  - R-HSA-190408 (FGFR2b autophosphorylation) gives 10 copies.
    - The long-isoform input `190227{FGF3}` produces only
      `190411{189967=FGF3}`, never the short form.
    - The catalyst node equals the input node (coupling).
  - Synthetic tie: two output members equally similar. Expect two copies,
    with no dependence on member order (shuffle the mock's return order and
    assert identical output).
- **Determinism.** Two builds of FGFR2 and RAF under different hash seeds.
  The sets of `(source key, target key, edge_type)` triples are equal, where
  a VR's key is its reaction stId plus its sorted binding.

---

## 6. Measurement plan (the arm, pre-registered)

1. **Before code: census per §1.4.**
   - homo vs multiset;
   - full vs Option F;
   - number of reactions in each §2.4 fallback step;
   - **the I3b / unconnected-producer count**, which decides how much of the
     win depends on the [DECIDE] in §3.1.
2. **Build.** Generator flag `LNG_VARIANT_NODES=1` (plus
   `LNG_SET_STOICH=homo`) against the canonical control at the same
   generator commit, through `scripts/catalog.sh`. Fresh cache dirs, or the
   cache-adoption trap applies.
3. **Structural gates before scoring.**
   - I1-I5 tests green;
   - `excess_classify3` gives 0 cut hand-offs;
   - the number of pathways in which a catalyst or regulator variant node is
     unconnected;
   - node and edge growth per pathway;
   - number of SCCs and the largest SCC per pathway (the expansion must not
     re-weld specs/018 loops through variant catalysts).
4. **Score both axes, held-out and tuning,** with McNemar, pathways and
   distinct readouts moved, and convergence. Pre-registered directions:
   - **RAF:** experimental up. R-HSA-5672972's neighbours lose dissolved
     scaffolds, and catalyst/input coupling stops the squaring of specs/036.
   - **FGFR1-4:** held-out up. HS and FGF cross-feeding between the b and c
     isoforms disappears, and hand-offs are restored.
   - **PIP3:** within noise.
   - **Pathways with set catalysts:** a member knockout weakens
     (pool-product → copy-mean). That is the spec's stated direction, so
     list them per pathway first.
5. **Trace RAF, PIP3 and FGFR2** against MP-BioPath's networks to the first
   divergence, as `spec.md` order of work step 4 already says.
