# specs/046: variant nodes end to end (complexes, catalysts and regulators expanded by their sets)

**Requirement (Adam, 2026-10-04).** A complex that contains sets breaks into
its variants: sets of sizes n and m give n × m complexes. Each variant gets
its own virtual reaction (complex_n → reaction_n), multiplied by the
alternatives of the reaction's other set inputs. Catalysts and regulators are
broken apart by their sets the same way. A downstream reaction must never read
a complex as separate root copies of its subunits when, in Neo4j, that complex
is the output of a reaction.

## What the generator does today (code-read and measured, 2026-10-04)

1. **Matching layer: as required.** `break_apart_entity` decomposes a complex
   containing sets into its combinations, and virtual reactions are made per
   combination, capped at 512 (`LNG_MAX_VARIANTS`). Past the cap there is one
   OR pool per set (`LNG_CAP_POOLS`, default since specs/045).
2. **Emission: collapsed.** `_map_annotated_entity_to_nodes` maps every
   variant of a complex onto ONE node, the plain stable id, under
   `LNG_COMPLEX_AS_NODE=1` (LNG #47, July 2026). Every copy reads the same
   complex node, so complex_n → reaction_n is lost. It was done because
   catalysts and regulators were not expanded into variants: a complex
   produced as variant n was a different node from the same complex used as
   a catalyst (the "catalyst silo"). Collapsing joined them.
3. **Set members that are complexes containing sets: dissolved** (code review
   F10). The set branch maps through `_matching_leaves`, which splits such a
   member into its components.
   - Measured on build `20261004-1425_1491276_cappools`
     (`~/deltasignal-catalogs/analysis/046/excess_classify3.py`): 13,048
     required inputs beyond the curated ones, in 3,933 copies.
   - Of these, 6,657 are inside a complex the pathway produces. In 2,504 of
     them the required subunit node is a ROOT: **the hand-off is cut.** The
     rest are fed by a reaction (3,774), fed only by boundary edges (379), or
     are dissolved roots (6,391).
   - The cuts are worst in FGFR2 (544), R-HSA-194315 (528), FGFR3/4/1,
     Mitotic G1 (162), DDX58 (142) and EPH-ephrin (120).
4. **Catalyst and regulator sets: one pool node** (`LNG_SET_POOL`,
   specs/033), combined by `product` in the solver.

## Design

Every participant that contains a set (input, output, catalyst, positive or
negative regulator) is expanded into its variants. A reaction has one copy
per combination over all its participants, and the copies recombine as OR
wherever their outputs are consumed. Each variant is one node.

**The variant key is the crux.** A variant is named by the outer stable id plus
the member chosen in each set inside it, nested sets included. It is computed
by one function, used on the producing side and the consuming side alike, so
the same variant gets the same node wherever it appears. Every break today
(F10, copy splits, the catalyst silo) is a place where the two sides name a
variant differently.

**Past the cap:** OR pools (specs/045), unchanged. **The boundary layer**
(specs/044): unchanged.

## Invariants (the tests)

- **I1.** Every copy uses exactly one member of every set it contains,
  including nested sets, catalysts and regulators. (Fails today: 156
  same-set excess inputs, F1.)
- **I2.** A producer's output variant node and a consumer's input variant node
  with the same key are the same node.
- **I3.** Zero cut hand-offs: no copy requires a root subunit node of a
  complex the pathway produces (2,504 today). `excess_classify3.py` becomes
  the regression check.
- **I4.** Past the cap, one OR pool per set.
- **I5.** `boundary_edges.csv` unchanged in kind.

## Stated in advance: predictions will move

**Catalyst and regulator set semantics change.**
- Today a set-valued catalyst is one pool combined by `product`, so a
  knockout of one member stops every copy.
- Under expansion it stops only that member's copies, and the OR over copies
  keeps the rest. That is close to `mean`, which specs/033 and 038 measured
  as neutral (`max` was refuted).
- Pre-register the direction per pathway before the arm.

## Order of work

1. **Size census.** Count the copies and nodes catalog-wide under full
   expansion, and how often the cap is hit, before writing code.
2. **Blind derivations** (Opus and Fable) of the variant-key function and of
   the catalyst/regulator expansion, against the invariants.
3. **Implementation in the generator,** behind a flag, with I1–I5 as tests.
4. **A catalog arm against the canonical build**, pre-registered, both axes.
   Trace RAF, PIP3 and FGFR2 against MP-BioPath's networks to the first
   divergence.

## Merged design (Adam, 2026-10-09: agreed with all six recommendations)

From `derivation-opus.md` (A) and `derivation-fable.md` (B), written blind.

**Common to both, adopted as derived:**
- Variant key `stId::variant::{slot=member;…}`: sorted, a pure function of
  Neo4j structure and the copy's choices. A bare set takes its chosen
  member's key; modifier-isoform sets and capped entities are their plain
  stable id.
- One shared `parse_variant_key` / `variant_leaves` replaces every tail parser.
- An entity used in several roles (for example input and catalyst) is one
  choice per copy.
- Output binding, in order: identity, then shared reference entities at the
  isoform level, then one copy per alternative. Never by position (replaces F9).
- Root variants are decomposed through the copy's chosen members (no F7 set
  node for a resolved slot). Terminal variants get sinks over their chosen
  leaves.
- Set pools (specs/033) are replaced below the cap and kept past it.
  `LNG_COMPLEX_AS_NODE` is superseded and mutually exclusive with the new
  flag `LNG_VARIANT_NODES`.
- `reaction_connections` gets a deterministic order.
- The DeltaSignal side gets per-variant rows in `containment.csv`, `drugs.csv`,
  `cofactors.csv`, `pools.csv` and `node_resolution.csv`, and pins and
  readouts aggregate over a gene's variant nodes.

**Decisions:**

| | Decision |
|---|---|
| D1 | I2 is scoped to **connected pairs** (precedingEvent or diagram), plus catalysts and regulators. Global identity by key is a separate arm. |
| D2 | **One choice per set per reaction**, across inputs, outputs, catalysts and positive regulators (B). |
| D3 | A set appearing more than once in one complex uses the **same member each time by default** (homo). Heteromers (for example BRAF:RAF1) are an arm, `LNG_SET_STOICH=multiset`. |
| D4 | **Negative regulators are expanded like every other participant** (A). |
| D5 | **A pool past the cap combines by mean**: OR, exactly equal to expanding (A). It needs its own edge type and an arm, because today's cap pools use `product` and the RAF +14 was measured under product. |
| D6 | **A capped reaction falls back in this order** (A): pool the participants whose choice reaches no output; then one copy per output variant; then a single copy reading one pool per varying participant. |

**Tests:** I1 to I5, plus:
- I3b (A): no input variant is a root when a connected reaction produces its key;
- B's split of a "root subunit" into a total cut against "another leaf conducts";
- B's T1–T13, among them determinism under two `PYTHONHASHSEED` values,
  catalyst equal to input, ties, stoichiometry 2, and flag-off byte
  equivalence;
- A's fixtures: FGFR2 R-HSA-190408 → R-HSA-5654404 → R-HSA-5654407, and RAF
  R-HSA-5672972.
