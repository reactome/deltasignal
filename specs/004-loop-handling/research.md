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
3–9. At the default threshold of 15 this classifies **12 of the 34
components as recycling artifacts, holding 4,205 of the 4,561 cycle-resident
nodes (92%)**; the other 22 components hold 356. The partition is identical
at threshold 10 and moves two components at threshold 20. Cell Cycle Checkpoints' 643-node component is **three curated reactions
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

---

# Phase 3 result — the control (T011–T013)

## R6 — The propagator is not the problem. The networks are.

DeltaSignal's propagator run on MP-BioPath's own published networks, same
cases, default solver config. Common case set = scored in all three arms
(562 of 847).

| arm | correct | accuracy | macro-F1 | F1 DOWN | F1 NO_CHANGE | F1 UP |
|---|---|---|---|---|---|---|
| **DeltaSignal on MP-BioPath's networks** | **403/562** | **0.7171** | 0.6443 | 0.812 | 0.320 | **0.801** |
| MP-BioPath's published predictions | 405/562 | 0.7206 | 0.6701 | 0.820 | 0.402 | 0.788 |
| DeltaSignal on our networks | 364/562 | 0.6477 | 0.5787 | 0.753 | 0.277 | 0.706 |
| shortest signed path on our networks | 391/562 | 0.6957 | 0.6171 | 0.812 | 0.280 | 0.759 |

**Given the same networks, DeltaSignal and MP-BioPath are two cases apart.**
The 41-case deficit is not the propagator; **it is entirely network
structure**. That answers SC-001 and it redirects the feature: the question
is no longer "is our maths worse" but "what about our representation costs
39 cases".

DeltaSignal's macro-F1 is still 0.026 behind on the same networks, and all
of it is `NO_CHANGE` (F1 0.320 against 0.402) — we over-call change. Our
`UP` F1 is *better* (0.801 against 0.788). That residual is a real but
second-order finding.

### Adversarial checks — the result survives all six

This looked too clean, so it was attacked before being believed.

1. **Not reading back a pinned input.** 12 of 845 control cases had a fully
   pinned readout and all 12 were excluded — a *lower* rate than our own
   arm's 49 of 564.
2. **Not just reproducing MP-BioPath.** The two agree on only **477 of 562
   (84.9%)**. DeltaSignal is right where MP-BioPath is wrong on 36 cases and
   wrong where it is right on 38. Same score by a genuinely different route,
   not a reimplementation.
3. **Convergence: 0 of 845 non-converged**, against 179 of 564 on our
   networks. Direct confirmation that our non-convergence is caused by the
   cycles and not by the solver.
4. **No prediction bias.** DOWN/NO_CHANGE/UP is 200/135/227 on their
   networks and 203/131/228 on ours — near-identical, so the gain is not the
   distribution shift that caught us in feature 002. Both over-predict
   NO_CHANGE against an actual 71.
5. **Per pathway, DeltaSignal-on-their-networks beats DeltaSignal-on-ours in
   9 of 10 pathways**, so it is not one pathway carrying the result.
6. **It also beats MP-BioPath's own published predictions on four
   pathways** — S_Phase 14/15 vs 12/15, ERBB2 18/23 vs 16/23, WNT 25/35 vs
   21/35, TP53 159/232 vs 154/232.

| pathway | DS on their nets | DS on our nets | MP-BioPath published |
|---|---|---|---|
| Cell_Cycle_Checkpoints | 42/44 | 35/44 | 42/44 |
| HDR | 25/32 | 24/32 | 28/32 |
| Mitotic_G1-G1_S_phases | 31/73 | 21/73 | 31/73 |
| Mitotic_Prophase | 18/20 | 17/20 | 18/20 |
| **PIP3_activates_AKT_signaling** | 67/84 | **75/84** | 79/84 |
| RAF_MAP_kinase_cascade | 4/4 | 2/4 | 4/4 |
| S_Phase | 14/15 | 10/15 | **12/15** |
| Signaling_by_ERBB2 | 18/23 | 13/23 | **16/23** |
| Signaling_by_WNT | 25/35 | 23/35 | **21/35** |
| Transcriptional_Regulation_by_TP53 | 159/232 | 144/232 | **154/232** |

**PIP3 is the counter-example and should not be buried**: our network is
better there, 75 against 67. Whatever our representation does right, it does
it in PIP3, and a blanket "their networks are better" claim is wrong.

### The second half of the deficit is coverage, not accuracy

| | cases scored |
|---|---|
| on MP-BioPath's networks | **845 of 847** |
| on our networks | **564 of 847** |

283 cases are scoreable on their networks and not on ours. The breakdown of
why, on our side, is not what "our networks are worse" would predict:

- **204 are `proxy_available_not_enabled`** — the readout *is* reachable
  through a producing-reaction proxy and the benchmark is configured not to
  use it. That is a harness flag, not a network defect, and it is the single
  largest bucket.
- **76 are `exact`** — the readout mapped fine, so the loss is upstream: no
  perturbable root input, or no path.
- **3 are `absent_from_network`.**

So a third of the case set is being discarded before the propagator is
consulted, and most of that is a switch. This needs measuring before any
more effort goes into loop interventions.

## What this does to the plan

- **US1 is answered and the answer is unambiguous.** Do not spend further
  effort attributing; spend it on the network.
- **Loops are implicated but are not the whole story.** They own the
  convergence difference outright (179 → 0). Whether they own accuracy is
  still US3's question, and the control does not settle it, because their
  networks differ from ours in more than acyclicity: they are also a third
  the size, carry only 1.8% negative edges, and resolve readouts our harness
  declines to resolve.
- **A new candidate has appeared that is cheaper than any loop
  intervention**: the 204-case proxy bucket. It belongs in this feature's
  reporting because it changes the denominator of every comparison, but the
  decision to enable proxies is a benchmark-methodology change and is Adam's
  call, not one to make silently mid-feature.
