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

---

# Why the accuracy is low — the decomposition (post-005)

Measured on the corrected 742-case set, same propagator on both network sets.

## R10 — RETRACTION first: connectivity is 18.7% vs 15.9%, not 41% vs 16%

An earlier pass in this analysis reported that our networks fail to connect
41% of perturbation–readout pairs against MP-BioPath's 16%, and that 215
cases had no path in ours but did in theirs. **That was an artifact of the
analysis script, not of the networks**: it resolved a readout only through
`nodes.csv`, so every set-valued readout counted as "readout absent" and
inflated our no-path rate by the whole 178-case set population. Corrected,
with set members resolved:

| | has a directed path | no path |
|---|---|---|
| our networks | 603 / 742 (**81.3%**) | 139 (18.7%) |
| MP-BioPath's | 711 / 845 (**84.1%**) | 134 (15.9%) |

Connectivity is **comparable**. The tell I should have caught immediately:
178, the inflation, is exactly the number of set-resolved cases recovered in
feature 005.

## R11 — Loops are NOT the cause of the accuracy gap

Feature 004's hypothesis is answered, and the answer is no.

| | cases we get wrong that their networks get right | cases both get right |
|---|---|---|
| readout is cycle-resident | **7.5%** | 9.7% |
| solve did not converge | **20.6%** | 33.4% |

The cases we lose are **less** cyclic and **less** often non-converged than
the cases we win. Loops demonstrably cause the non-convergence — 179 of 564
on our networks against 0 of 845 on their acyclic ones — but non-convergence
is not what costs accuracy. Accuracy by convergence makes the same point
backwards: non-converged cases score 0.796 and converged ones 0.614.

So the loop interventions in US3 should not be expected to move accuracy.
They remain worth doing for correctness and to make the solver honest; they
are not the accuracy lever.

## R12 — The two real causes

**Cause 1, the no-path wall — 139 cases (18.7%), and we get 30 of them.**

| | DOWN | NO_CHANGE | UP |
|---|---|---|---|
| ground truth says | 50 | 30 | 59 |
| we predict | 0 | **139** | 0 |

With no directed path the model correctly answers "no change" on every one,
and the ground truth says the readout moved in **109 of them**. This is a
hard ceiling for any directed causal model, not a modelling error:
MP-BioPath scores **0.281** on these same 139 cases and also answers
NO_CHANGE 105 times. It is the co-regulation wall — perturbation and readout
are co-descendants of a shared hub, or the truth reflects feedback and
indirect regulation a directed model does not produce.

**Cause 2, accuracy where a path does exist — 0.765 against their 0.816.**

Here we over-call change badly: on 603 path-having cases we predict
NO_CHANGE **7 times** against 50 actual. Nearly every path is treated as a
conduit that must transmit.

**Decomposition of the 0.6617 → 0.7358 gap:**

| counterfactual | accuracy | gain |
|---|---|---|
| ours as measured | 0.6617 | — |
| give us their **path-accuracy**, keep our connectivity | 0.7038 | **+4.2 pts** |
| give us their **connectivity**, keep our accuracies | 0.6777 | +1.6 pts |
| their measured | 0.7358 | — |

**Path-accuracy is ~2.6× the lever connectivity is**, and connectivity is
where nearly all the effort has historically gone.

## What this says to do next

1. **Stop treating no-path cases as losses to fix by adding edges.** 18.7% of
   the case set is unreachable and MP-BioPath is barely better there. Past
   attempts to close it by adding connectivity regressed the benchmark
   (cross-pathway stitching, diagram descend). The honest move is to report
   this subset separately as the directed-causal ceiling.
2. **The lever is discrimination on path-having cases**, specifically the
   refusal to ever say NO_CHANGE when a path exists (7 of 603). A path is
   not an obligation to transmit. This is where the remaining ~31 cases are.
3. **004's loop interventions are correctness work, not accuracy work**, and
   should be framed and measured as such.

## R13 — Overfitting audit (Adam: "make sure we are not overfitting")

Three checks passed, one correction to R12, and three risks that stand.

### PASSED: feature 002's defaults generalise

The AND/clamp default change was chosen by A/B on this case set, so it is
the most exposed decision. Split by the pre-declared pathway split:

| axis | development | held out |
|---|---|---|
| experimental | +8 cases, macro-F1 **+0.0045** | **+23 cases, macro-F1 +0.0375** |
| curator | +125 cases, +0.1275 | +24 cases, +0.0297 |

The experimental gain is **five times larger on held-out than on
development**, and the curator axis — a different ground truth — moves the
same way. That is the opposite of the overfitting signature.

### PASSED: the classification cutoffs are not tuned to this data

| cutoffs | correct | accuracy | macro-F1 |
|---|---|---|---|
| 0.60 / 1.40 | 481 | 0.6482 | 0.5720 |
| 0.80 / 1.20 | 487 | 0.6563 | 0.5747 |
| **0.85 / 1.15 (default)** | **491** | 0.6617 | 0.5771 |
| 0.90 / 1.10 | 492 | 0.6631 | 0.5782 |
| 0.95 / 1.05 | 494 | 0.6658 | **0.5803** |

The curve is smooth and monotonic toward tighter cutoffs and **our default is
not the peak**. A tuned threshold would sit on the maximum; ours sits below
it, because it is inherited from MP-BioPath's published convention rather
than fitted. The available gain is +3 cases, inside noise — do not chase it.

### PASSED: the decomposition survives leave-one-pathway-out

| dropped | n | accuracy | % with path | acc on path | acc no-path |
|---|---|---|---|---|---|
| (full set) | 742 | 0.662 | 0.813 | 0.765 | 0.216 |
| TP53 | 493 | 0.688 | 0.854 | 0.765 | 0.236 |
| PIP3 | 542 | **0.605** | 0.744 | 0.739 | 0.216 |
| Mitotic_G1 | 667 | **0.702** | 0.850 | 0.787 | 0.220 |
| …the rest | | 0.653–0.667 | | 0.756–0.774 | 0.202–0.224 |

**Accuracy on no-path cases is 0.202–0.236 whichever pathway you remove, and
on path-having cases 0.739–0.787.** The R12 mechanism is not an artifact of
one pathway. The *headline rate*, by contrast, swings 0.605–0.702 purely on
composition — so the decomposition is trustworthy and the single number is
not.

### CORRECTION to R12

R12 said we "over-call change" and predict NO_CHANGE only 7 times in 603.
That figure is right for path-having cases, but the framing was too broad.
Across all 742 cases we predict NO_CHANGE **146 times against 80 actual** —
we **over**-predict it overall, because every one of the 139 no-path cases
gets NO_CHANGE by construction.

And it is not a uniform property. Predicted vs true NO_CHANGE per pathway:
PIP3 **0 vs 7**, HDR **0 vs 3**, Mitotic_Prophase 1 vs 1, RAF 1 vs 1, ERBB2
4 vs 4, WNT 6 vs 2, CCC 11 vs 4, S_Phase 14 vs 2, Mitotic_G1 40 vs 28, TP53
69 vs 28. So "a path is treated as an obligation to transmit" is a **PIP3 and
HDR** property, not a global one, and the proposed lever is narrower than
R12 claimed.

### RISKS THAT STAND

1. **Two pathways are 60.5% of the case set** — TP53 249 cases (33.6%) and
   PIP3 200 (27.0%). Every headline number is essentially a weighted average
   of two pathways, which is why leave-one-out moves it by ±5 points.
2. **The held-out split is not clean**, and the manifest now says so:
   `held_out_was_used_for_configuration: True`. Held-out results are
   therefore better than development results at establishing generalisation,
   but they are not a virgin test set.
3. **Accumulated researcher degrees of freedom.** Dozens of configurations
   have been A/B'd against this same 10-pathway set across the project. No
   individual A/B accounts for that multiplicity, and a +3-case result on
   this set means very little on its own.

**Mitigation to take before acting on R12:** the 92-pathway catalog already
exists and was used for the 2026-07 holdout (13,429/16,696 on 71 pathways).
Any change motivated by the path-accuracy finding should be validated there,
not on these ten — particularly since the finding is now known to be
concentrated in two of them.

## R14 — The "no-path wall" is an EXPERIMENTAL-ground-truth phenomenon, not a model failure

R12 called the no-path subset "a hard ceiling for any directed causal model".
That is right about the experimental axis and **wrong as a general claim**,
and the correction matters because it changes what is worth building.

Checked at 8.5× scale on the 92-pathway catalog, against curator ground
truth (72 pathways with both a network and a truth table, **19,984 scoreable
cases**). With no path the model answers NO_CHANGE deterministically, so
accuracy there is exactly the share of no-path cases whose truth is also "no
change" — computable without the solver, which is what made this cheap.

| ground truth | pathways | cases | no-path rate | truth is NO_CHANGE on those |
|---|---|---|---|---|
| curator, 92-catalog | 72 | 19,984 | **61.6%** | **11,050/12,304 = 0.898** |
| curator, 10-catalog | 8 | 3,418 | 45.8% | 1,377/1,564 = 0.880 |
| **experimental**, 10-catalog | 10 | 742 | 18.7% | **30/139 = 0.216** |

**On the curator axis the directed model and the curators agree ~90% of the
time when no path exists.** There is no wall. The 0.216 is specific to
experimental data: in a real cell the readout moves in ~78% of cases where
Reactome has no directed route from the perturbation.

So the no-path subset does not measure a deficiency of our propagator, or of
MP-BioPath's. It measures the distance between Reactome's curated directed
causality and what a cell actually does — compensation, indirect regulation,
off-pathway effects, and whatever experimental confounding is in the 2019
measurements. MP-BioPath scores 0.281 on the same 139 cases, which is the
same wall from the other side.

### The apples-to-apples comparison, on our own path-having subset

| model | path-having (603) | no-path (139) |
|---|---|---|
| DeltaSignal | 461 = **0.7645** | 30 = 0.2158 |
| shortest signed path | 488 = 0.8093 | — |
| MP-BioPath | 507 = **0.8408** | 39 = 0.2806 |

Stated on identical cases, the gap where a path exists is **46 cases, 7.6
points** — larger than the 5.1 points R12 inferred from unequal subsets, and
it is the whole of the actionable deficit. We also lose to the model-free
traversal here (488 vs 461), on the subset where our model is supposed to
have an advantage.

### What to do with this

1. **Report the two subsets separately and say why.** A single accuracy
   figure over the experimental benchmark averages a tractable modelling
   problem with an intractable one, and 18.7% of it is not about the model at
   all. This is a manuscript framing point, not a metric trick — the curator
   numbers are the evidence that the split is principled.
2. **The target is 0.7645 → 0.8408 on path-having cases.** That is where
   every remaining winnable case lives.
3. **Do not chase the no-path subset with more edges.** Every past attempt
   regressed the benchmark, and the curator result now explains why: the
   edges are not missing, the effect is not directed-causal.

## R15 — Found it: de-repression is the single biggest systematic error

Decomposing the 142 errors on path-having cases:

| error | n |
|---|---|
| outright sign flip (UP↔DOWN) | **85 (60%)** |
| called a change on a NO_CHANGE readout | 50 |
| called NO_CHANGE on a real change | 7 |

**78% of path-having cases have both a positive and a negative route** to the
readout, so the network is usually ambiguous about direction and the model
has to weigh routes rather than follow one. 79 of the 85 sign flips are in
that ambiguous group.

The outlier that gave it away — accuracy split by the signs available:

| routes available | n | accuracy |
|---|---|---|
| positive only | 106 | 0.896 |
| both (ambiguous) | 472 | 0.758 |
| **negative only** | **25** | **0.320** |

Sixteen of those 25 are knockouts in `Mitotic_G1` whose truth is NO_CHANGE
and where we answer UP with values of 1.43, 9.997, 25.66 and 100.0 —
confidently wrong. MP-BioPath scores 21 of 25 on the same cases.

### It is systematic, not a 25-case curiosity

Across **all** knockout cases:

| our call | n | correct | |
|---|---|---|---|
| knockout → DOWN | 241 | 213 | **0.884** |
| **knockout → UP (de-repression)** | **78** | **28** | **0.359** |

MP-BioPath on those same 78: **60 correct (0.769)**, more than twice our
rate, and it calls them UP only 37 times against our 78. Truth on the 78 is
UP 28 / NO_CHANGE 22 / DOWN 28 — so the answer is not "never de-repress"
either; it is that the magnitude is not being graded.

**This is ~32 of the 46-case actionable gap in one error mode**, and it
spreads across six pathways (TP53 24, Mitotic_G1 20, PIP3 19, WNT 6, ERBB2 3,
RAF 3), so it survives the concentration risk that R13 flagged.

### Mechanism, confirmed in the code

`reaction_model.jl:1247` — `h_k = (bl + eps) / (x_inh + eps)` with baseline
`bl = 0.01` and `DS_INHIBITOR_EPS` defaulting to **0.001**. A knockout sets
`x_inh = 0`, giving `0.011 / 0.001 = 11×` de-repression from removing one
inhibitor, clamped to 10× per reaction but **compounding across a path** —
which is where the 25.66 and 100.0 values come from.

That 10× is an assumption that every inhibitor was suppressing its target
tenfold at baseline. Nothing measures it, and `DS_INHIBITOR_FLOOR` caps how
far an inhibitor can *suppress* while nothing caps how far removing one can
*raise*.

The principled statement: de-repression should be bounded by what the
activators can supply. With activators at baseline, removing an inhibitor
should return the target to its uninhibited baseline — not to ten times it.

### Status

Probing `DS_INHIBITOR_EPS` (0.01 gives ~2× maximum de-repression, 0.003 gives
~4.3×) as a single-parameter test of the mechanism. **This is a probe of the
diagnosis, not a tuning exercise**: per R13, anything adopted must be
validated on the curator axis and the 92-pathway catalog rather than chosen
on these 742 cases.

## R16 — The de-repression probe is confirmed on BOTH ground truths

| axis | baseline (ε=0.001, 11×) | ε=0.01 (2×) | delta |
|---|---|---|---|
| experimental (742) | 491, macro-F1 0.5771 | **506, 0.5986** | **+15 cases, +0.0215** |
| **curator (4,346)** | 3,016, macro-F1 0.6711 | **3,075, 0.6867** | **+59 cases, +0.0156** |

Two independent ground truths, both positive. On the experimental axis the
held-out gain (+9 cases, +0.0259) **exceeds** the development gain (+6,
+0.0150), which is the generalisation direction, not the overfitting one.

Per-pathway on curator: TP53 +30, PIP3 +14, WNT +11, Cell Cycle Checkpoints
+6, Mitotic_G1 −2. Four pathways up, one marginally down.

And the targeted error mode behaves exactly as the diagnosis predicts:
spurious knockout→UP calls fall **78 → 68** while the correct ones hold at
**28**, so accuracy on that subgroup goes 0.359 → 0.412. It is removing wrong
calls, not trading them.

**Honest caveat on magnitude**: only 9 of 29 changed experimental predictions
and 34 of 239 curator ones converged in both arms, so some of the delta is
the known TP53 sweep-order noise. The *direction* is supported by two ground
truths, a monotonic ε sweep (11× → 4.3× → 2× giving 491 → 495 → 506) and a
mechanism confirmed in the code. The *exact magnitude* is not yet established.

### Recommendation — and what NOT to do

Do **not** simply move `DS_INHIBITOR_EPS` to 0.01 and call it done. ε is a
tuning constant, this sweep was run on the evaluation set, and R13 says
precisely that is the thing not to trust. Two things follow:

1. **Spec the principled version.** The defect is that de-repression is
   unbounded above while suppression is bounded below. The fix is to bound
   de-repression by what the activators can actually supply — removing an
   inhibitor should return a target to its uninhibited level, not to ten
   times baseline. ε=0.01 approximates that by accident, at 2×.
2. **Validate on the 92-pathway catalog**, which now carries the resolution
   exports, before any default changes. The finding spans six pathways, which
   is encouraging, but 60% of this case set is still two pathways.

Also unanswered deliberately: whether the optimum is beyond ε=0.01. Sweeping
further on these 742 cases is exactly the move R13 warns against; answer it
on the larger catalog or not at all.
