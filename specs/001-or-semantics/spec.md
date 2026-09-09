# Feature Specification: Faithful OR semantics for set-valued catalysts

**Status:** investigated, prototyped, NOT enabled by default
**Repos:** logic-network-generator (representation) + deltasignal (processing)
**Feature branches:** `fix/adversarial-review` in both (LNG PR #56, DS PR #12)

## Problem

The generator derives an edge's logic operator from its **sign**, not from curated
semantics — `and_or = "and" if pos_neg == "pos" else "or"`. The shipped 92-pathway
catalog therefore contains **zero** `pos/or` regulatory edges:

| edge | count |
|---|---|
| `catalyst / pos / and` | 33,531 |
| `regulator / pos / and` | 3,849 |
| `regulator / neg / or` | 2,409 |

That inverts curator intent for **set-valued** regulators. A Reactome `DefinedSet`
catalyst means "any one of these plays this role", but flattening it to N members
all marked `and` asserts every isoform is simultaneously required. Example:
catalyst `R-HSA-5622009` "ITPR tetramers" (a DefinedSet of the ITPR1/2/3 tetramers)
on reaction `R-HSA-169680` emits three `pos/and` edges.

Compounding it on the solver side, the per-edge inhibitor AND/OR flag was parsed
into `Reaction` and then **dropped** (`IndexedReaction` had no field for it), so
the only correctly-OR edges in the catalog were aggregated conjunctively anyway.
Net effect before this work: **no regulatory OR semantics survived the pipeline
end to end.**

## Scope

EntitySets are normally **split** — each alternative becomes its own virtual
reaction — so within a VR the member genuinely *is* required and `and` is correct.
The exception is the sets we do **not** split on: `append_regulators` bypasses the
VR machinery and flattens a set-valued regulator's members onto the reaction's
single existing VR.

Verified on a 3-pathway Release97 regeneration: all **117** OR-marked catalyst
groups have **≥2 members co-located on the same target VR** (2–26 members each),
**zero** single-member groups.

- **In scope:** catalyst and positive-regulator edges derived from an EntitySet.
- **Out of scope:** input/output edges (splitting already encodes their OR-ness —
  marking them `or` would be wrong), and nested complex-containing-a-set boundary
  assembly edges (needs per-member provenance from `_decompose_regulator_entity`).

## Architectural principle (project rule)

LNG represents pathways **as curators intended to design them**; DeltaSignal
decides how best to process them. Corollary: when a faithful LNG representation
regresses the benchmark, fix DeltaSignal's *processing* rather than encoding an
interpretation into LNG. Any change to regulator/catalyst connections must be
validated through DeltaSignal.

## User Scenarios & Testing

### User Story 1 — Faithful representation (Priority: P1)
A curator marks a reaction's catalyst as a DefinedSet of three isoforms. The
generated logic network should record those members as **alternatives**, not as
three co-required inputs, so a downstream consumer can choose how to treat
redundancy.

### User Story 2 — Perturbation of one redundant paralog (Priority: P1)
A single-gene knockout of one isoform of a set-valued catalyst should reduce the
reaction's throughput proportionally to the lost capacity — not be ignored
(`max` semantics) and not abolish the reaction outright (`and` semantics).

### User Story 3 — No accuracy regression (Priority: P1)
Adopting the faithful representation must not degrade agreement with the
MP-BioPath experimental ground truth.

### Edge Cases
- A single-member decomposition: operator is a no-op; must not be mis-marked.
- Modifier-isoform sets (ubiquitin/SUMO/NEDD8) are treated atomically and
  decompose to one member — unaffected.
- Sets bundled by the `LNG_MAX_VARIANTS` cap collapse to one opaque node.
- A reaction with an OR cluster **and** AND-clustered inputs — the combination
  rule decides whether the OR loss is visible at all (see plan.md).

## Requirements

### Functional Requirements
- **FR1** LNG must be able to emit set-derived positive regulator members with
  `and_or = "or"` (`LNG_SET_MEMBERS_OR`).
- **FR2** The flag must participate in the per-pathway cache fingerprint so
  flipping it invalidates stale caches.
- **FR3** DeltaSignal must forward the per-edge inhibitor AND/OR flag into
  `IndexedReaction` and honor it (`DS_INHIBITOR_OR`).
- **FR4** DeltaSignal must offer an aggregator for redundant alternatives that
  is neither masking (`max`) nor abolishing (`and`) — `DS_OR_MODE=capacity`
  with a tunable redundancy weight `DS_OR_REDUNDANCY`.
- **FR5** DeltaSignal must offer a combination rule under which an OR cluster can
  **gate** a reaction (`DS_OR_COMBINE=gate`), because the default
  `max(and, or)` discards OR-cluster loss entirely.
- **FR6** All of the above default to current behaviour; none may change shipped
  results without an explicit decision.

### Key Entities
- **OR cluster** — the set of activator slots on one reaction whose `is_and` is
  false; sourced from one set-valued regulator's members.
- **Redundancy weight `w`** — `fold = w·mean(fold_i) + (1−w)·min(fold_i)`.
  `w=1` credits full redundancy (≡ existing `mean`); `w=0` makes any lost
  alternative decisive.

## Success Criteria

### Measurable Outcomes
- **SC1** Faithful representation costs **no** significant accuracy loss against
  the MP-BioPath experimental set (223 scored cases). — **MET** (gate, `w=0.5`:
  169/223 vs control 167/223, macro-F1 0.7346 vs 0.7276, McNemar p=0.688).
- **SC2** A KO of one of N redundant alternatives produces a detectable decrease
  (fold < 0.85 DOWN cutoff) for realistic N up to 26. — **MET** under gating;
  **NOT met** under `max` (fold stays 1.0 at every `w`).
- **SC3** Input/output edge semantics unchanged. — **MET** (change confined to
  `append_regulators`).
- **SC4** No default behaviour change. — **MET** (all four knobs default off).

## Assumptions
- Reactome `DefinedSet`/`CandidateSet` membership expresses genuine functional
  redundancy. In a given cell line only a subset may be expressed, which is why
  full redundancy credit (`w=1`) under-predicts real knockout effects.
- The 0.85/1.15 DOWN/UP cutoffs are the detection thresholds that matter.
- Only 3 pathways (R-HSA-1257604, R-HSA-453279, R-HSA-69620) are eligible in the
  MP-BioPath harness, so all measurements here are on a narrow slice.
