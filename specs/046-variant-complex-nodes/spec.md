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
