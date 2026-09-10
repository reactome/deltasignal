# Research: Handle loops properly, and compare honestly

Phase 0. All three unknowns are resolved by measurement, not assumption.
Everything below is on the shared catalog build (`$S/cat`, ten pathways, LNG
main 4ff0408, Reactome Release97) and MP-BioPath's published networks in
`~/gitroot/mp-biopath-pathways/pathways/`.

## R1 — The acyclicity asymmetry is real and larger than expected

**Decision**: treat the MP-BioPath networks as acyclic for all purposes.

Directed-graph SCC analysis of both network sets on the same ten pathways:

| pathway | MPB largest SCC | LNG largest SCC | LNG % nodes cycle-resident |
|---|---|---|---|
| Cell_Cycle_Checkpoints | **0** | 643 | 52.1% |
| HDR_through_Homologous_Recombination | **0** | 834 | 40.0% |
| Signaling_by_ERBB2 | **0** | 465 | 37.0% |
| Signaling_by_WNT | **0** | 211 | 9.0% |
| Transcriptional_Regulation_by_TP53 | 4 | 836 | 36.2% |
| S_Phase | 5 | 52 | 14.0% |
| PIP3_activates_AKT_signaling | 7 | 97 | 12.4% |
| Mitotic_Prophase | 8 | 28 | 16.9% |
| RAF_MAP_kinase_cascade | 11 | 443 | 44.3% |
| Mitotic_G1_phase_and_G1_S_transition | (not published) | 86 | 10.5% |

Zero self-loops in any MP-BioPath network. Four of the nine published are
fully acyclic; the largest component anywhere is eleven nodes. Catalog-wide
the LNG networks hold 4,561 cycle-resident nodes.

Adam's account is confirmed by the artifacts. **Consequence: the published
365-versus-407 comparison is not like-for-like and never has been.**

## R2 — NEGATIVE, and a bug I nearly reported that is not one

**Decision**: no identifier-mapping change is needed. The control can be run
directly, because MP-BioPath's network node identifiers are Reactome database
identifiers and the benchmark cases already key on those.

`load_dbid_to_uuids` in `bench/benchmark_mpbiopath_cases.py` derives what it
calls a `dbid` by taking the numeric suffix of a **stable id**
(`R-HSA-975977` → `975977`). Those are different identifier spaces, and
Reactome-wide only **73,334 of 410,271 PhysicalEntities (17.9%)** have
`dbId == stId` suffix. That looked like a live defect affecting every number
we have published.

It is not. Restricted to the entities that actually appear in the ten-pathway
catalog, the identity holds for **4,022 of 4,022 (100%)**, and of the 151
distinct `key_output_dbid` values in the scored cases, zero resolve to a
different entity under the two interpretations and zero are missed. The
Reactome-wide divergence is in entities that never enter these networks.

So the naming is misleading and worth fixing for the reader's sake, but the
mapping is correct and no result depends on changing it. Recording this
because the global statistic is alarming and the next person will find it too.

For the control specifically no mapping is needed at all: 211 of 212 ERBB2
node ids in the MP-BioPath file appear in `db_id_to_name_mapping.txt` as
`Database_Identifier`, which is the same key the cases already use.

## R3 — The giant components are recycling artifacts, and the discriminator is countable

**Decision**: classify a component as a recycling artifact by its
**nodes-per-distinct-reaction-stable-id ratio**, not by polarity alone.

A genuine feedback loop is a chain of distinct curated reactions, so its node
count and its reaction count are of the same order. A recycling artifact is
one or two curated reactions instantiated as hundreds of virtual reactions,
so the ratio is large.

| pathway | SCC nodes | intra edges | negative | distinct reaction stIds | nodes per reaction |
|---|---|---|---|---|---|
| Cell_Cycle_Checkpoints | 643 | 1120 | **0** | **4** | **160.8** |
| RAF_MAP_kinase_cascade | 350 | 1099 | 5 | 8 | 43.8 |
| Signaling_by_WNT | 211 | 498 | **0** | 6 | 35.2 |
| Transcriptional_Regulation_by_TP53 | 836 | 2390 | **6** | 34 | 24.6 |
| Signaling_by_ERBB2 | 465 | 1686 | 121 | 19 | 24.5 |
| HDR | 834 | 4689 | **0** | 40 | 20.9 |
| — smaller components — | | | | | |
| Cell_Cycle_Checkpoints | 18 | 28 | 5 | 6 | 3.0 |
| PIP3_activates_AKT_signaling | 50 | 83 | 7 | 12 | 4.2 |
| Signaling_by_WNT | 23 | 40 | 8 | 4 | 5.8 |

The separation is clean: large components have ratios of 20–161, small ones
3–9. Cell Cycle Checkpoints' 643-node component is **three curated reactions
contributing exactly 320 intra-component edges each**, plus 160 diagram
bridges — one reaction replicated 320 times.

Polarity alone is a weaker discriminator and would have misled us: 21 of the
34 components contain no negative edge, but ERBB2's 465-node artifact
contains 121 negative edges. Use the ratio; report polarity alongside.

**Alternatives considered**: (a) polarity-only — rejected, see ERBB2;
(b) component size alone — rejected, it is a symptom not a cause and would
delete genuine large feedback if any existed; (c) trusting `precedingEvent`
annotation to orient the recycle edge (the 2026-06-12 loop-taxonomy proposal)
— still viable and preferable where the annotation exists, but it requires a
Neo4j round trip per edge and the ratio test is a cheaper first filter.

## R4 — MP-BioPath's four-column format

**Decision**: column 3 is polarity (`1` = positive, `-1` = negative), column 4
is the conjunction flag with **`0` = AND and `1` = OR**.

Confirmed against MP-BioPath's own reader, which documents "if there is only
one parent it will be an AND relation": of children with exactly one parent,
5,148 carry flag `0` and only 47 carry flag `1` (99.1%). Distribution across
all ten pathways is 17,513 `(pos, AND)`, 7,234 `(pos, OR)`, 440 `(neg, AND)`,
9 `(neg, OR)` — so **only 1.8% of MP-BioPath's edges are negative at all.**

That last figure is worth carrying into the comparison: our networks encode
far more inhibition than theirs, which is a second structural asymmetry
alongside acyclicity.

## R5 — The deficit is not non-convergence (carried from feature 002)

Of the 81 cases DeltaSignal gets wrong where MP-BioPath is right, **63 are in
solves that converged**. On the 179 non-converged cases the two score exactly
equally (143 each). Whatever loops are doing, they are doing it to the
converged answer. This is why feature 004 is separate from 003 and why fixing
the solver objective is not expected to close this gap on its own.

## Open questions carried into implementation

1. Does the control (DeltaSignal on MP-BioPath's networks) land near 407 or
   near 365? Unrun; it is the whole point of User Story 1.
2. Is `diagram_bridge` removal accuracy-neutral now, as it was when added?
   If it is neutral *and* removes 2,559 cycle-resident nodes, it is a free
   simplification even if it does not move the score.
3. Does breaking recycling artifacts help, hurt, or do nothing once the
   readout is still reachable?
