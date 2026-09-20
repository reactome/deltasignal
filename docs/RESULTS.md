# DeltaSignal: current results

**Last measured**: 2026-09-16, Reactome Release97, catalog build `cat_fresh`,
solver defaults as merged in `specs/009-solver-defaults` (`hill_log`,
assembly-limiting on, `DS_INHIBITOR_EPS=1e-12`).

Keep this file current. Every number here is reproducible from the commands at
the bottom; if you change a default, re-run and update the tables in the same
commit.

---

## 1. The headline: same networks, same ground truth

The fair comparison against MP-BioPath is on **MP-BioPath's own published
networks**, because our networks and theirs are built by different pipelines
from different Reactome versions. `bench/analysis/mpbiopath_network_adapter.py`
converts theirs into a DeltaSignal catalog; both tools are then scored against
the same Reactome curator ground truth on the same 72 pathways.

| | cases | accuracy | macro-F1 |
|---|---|---|---|
| MP-BioPath (published) | 23,009 | 84.58% | 0.8063 |
| **DeltaSignal, same networks** | 22,106 | **93.38%** | **0.9209** |

DeltaSignal scores 903 fewer cases (coverage). **Counting every one of those as
wrong**, it still leads: 20,642/23,009 = **89.71%**, a margin of **+5.14pp**.

### It is not carried by a few pathways

| | pathways | DeltaSignal | MP-BioPath | margin |
|---|---|---|---|---|
| transcriptional | 6 | 92.32% | 72.88% | **+19.4pp** |
| everything else | 66 | 93.79% | 89.09% | **+4.7pp** |

DeltaSignal wins in 37 pathways, loses in 21, ties in 14. **MP-BioPath degrades
16 points on transcriptional regulation** (median 64.3% against 93.9%
elsewhere); DeltaSignal holds 92–94% across both groups.

Worth stating plainly because a reviewer will find it: the *aggregate* margin is
concentrated. Excluding the 8 pathways most favourable to DeltaSignal it falls
to +2.95pp, and the median per-pathway difference is +0.008.

**That concentration has a known cause.** MP-BioPath solves a non-convex
optimisation — products for AND, reciprocals for inhibition — with a local
solver (Ipopt) from a fixed start, so it is not guaranteed to reach the
biologically correct root. Its own stored per-node output shows it landing on
**direction-inverted** solutions: in the pluripotency pathway a `FOXP1`
knockout drives targets to 23.03 and `FOXP1` over-expression drives them to
0.37, the reverse of the curator calls, on a network that is 98.6% activating
edges.

That is why that pathway scores **19.1%**. Its curator ground truth is balanced
— 164 down, 164 up, 132 unchanged — so answering "no change" to everything would
score **28.7%**. MP-BioPath scores below that floor, which is not what
under-calling looks like; it is predicting backwards. On the identical network
DeltaSignal called 156/164 knockdowns and 156/164 over-expressions correctly.

Measuring that inversion rate across the 66 pathways where both its raw output
and its published accuracy exist:

| direction-inversion rate | pathways | mean MP-BioPath accuracy |
|---|---|---|
| ≥ 25% | 10 | **64.2%** |
| < 10% | 38 | **93.0%** |

Pearson r = **−0.579**. The pathways where MP-BioPath scores badly are the
pathways where its optimiser inverts, and that is where our margin comes from.
DeltaSignal propagates causally in topological order (SCC condensation), so a
knockout cannot raise a downstream node across an activating edge; on the same
pluripotency network it called 156/164 knockdowns and 156/164 over-expressions
correctly.

So the defensible claim is **"comparable where MP-BioPath's solver converges,
decisively better where it does not"** — a mechanism, not a uniform +8.8pp.

### Failure profile on the same networks

| category | our networks | MP-BioPath's networks |
|---|---|---|
| `no_path` | 1,484 | 328 |
| `false_positive_change` | 1,441 | 597 |
| `propagator_missed` | 808 | 417 |

---

## 1b. The experimental axis, where the advantage does NOT hold

Everything above is **curator** ground truth — what a curator expects the
pathway to do. The other axis is **experimental evidence** — what actually
happened in a cell — and it tells a different story.

Same networks (MP-BioPath's), same 712 cases, the 8 pathways present in both
network sets:

| | cases | accuracy |
|---|---|---|
| MP-BioPath | 564/712 | **79.21%** |
| DeltaSignal | 555/712 | **77.95%** |

A 9-case difference: **level, not better.** DeltaSignal wins 5 of the 8 pathways
(RAF 100% vs 91.8%, Cell Cycle, Prophase, S Phase, WNT) and loses on PIP3, which
is 200 of the 712 cases and where it is 7.5 points down.

| axis | same networks | verdict |
|---|---|---|
| curator reasoning | **+9.03pp** | decisively better |
| experimental evidence | **−1.26pp** | level |

**This is the most important caveat in this document.** DeltaSignal reproduces
curator *reasoning* far better than its predecessor and is no better at
predicting experimental *outcomes*. A model that tracks the representation
better than the biology is exactly what a reviewer will suspect, and the honest
answer is that on this evidence we cannot rule it out.

Two structural limits on ever settling it:

- **Experimental truth exists for ten pathways only — the tuning ten.** There is
  no held-out empirical test available in this dataset, for either tool. That is
  a data-collection problem, not an analysis one.
- 712 cases across 8 pathways is small, and one pathway (PIP3) carries the
  difference.

On our own networks against experimental evidence: 627/849 = 73.85%, macro-F1
0.6599, against MP-BioPath's published 643/849 = 75.74%.

## 2. On our own networks

### 2026-09-20: the welded-loop fix (specs/018) — new canonical numbers

The generator's boundary expansion reused, as a root complex's subunit leaf, a
node that the root complex itself reaches downstream, welding cycles Reactome
does not have (1,994 of 2,077 cycle-carrying assembly edges; TP53's component
836 nodes, DSB's 1,127). reactome/logic-network-generator#91 (`79feca7`) refuses
exactly those reuses: stable-id edge multisets are identical in all 92
pathways, +336 fresh leaf nodes; TP53 → 56, DSB → 290. Regenerated catalog,
production solver, no solver change:

| split | pathways | cases | accuracy | macro-F1 |
|---|---|---|---|---|
| tuning (the paper's ten) | 11 | 5,100 | 0.770 | — |
| **held-out — quote this** | 71 | 19,000 | **0.870** | — |
| all pathways | 82 | 24,100 | **0.8471** | **0.8137** |

Against the previous production catalog on identical keys: held-out **+173**
(199 fixed / 26 broke, p < 1e-4), tuning +109, false change 1,437 → 1,150,
experimental +3 (conditioned), the 102 AKT-KO → TP53 cases unchanged at
100. **Stated plainly: 94% of the held-out gain is DSB Repair** (its welded
component was collapsing to zero under PARP1/2, FEN1, POLQ, RTEL1, XRCC5/6,
MUS81 knockouts) and 187 of the 199 fixes are DOWN → NORMAL where the truth
is NORMAL. It is adopted on correctness grounds; the accuracy gain is one
pathway's. The tables below this section predate the fix and describe the
previous catalog. The production catalog must be regenerated with the merged
generator for these numbers to reach the API.

Everything else measured in September 2026 — loop elasticity, the objective
minimiser, the loop pool, composition edges, catalyst/own-product rules — is
recorded in `specs/013`–`018` as measured and not adopted; the re-evaluation
of each on the fixed catalog is in `specs/018-derived-edge-loops/research.md`.


**Report the held-out split.** The MP-BioPath paper tuned on ten pathways and
reported on the rest. This project drifted off that protocol — the ten-pathway
set kept reversing decisions that held at scale, so decisions moved to the wide
curator set, which is the set we then report on. `bench/analysis/holdout_report.py`
restores the split at reporting time, over an existing dump.

| split | pathways | cases | accuracy | macro-F1 |
|---|---|---|---|---|
| tuning (the paper's ten) | 11 | 5,100 | 0.7322 | 0.7185 |
| **held-out — quote this** | 70 | 18,808 | **0.8552** | **0.8155** |
| all pathways | 81 | 23,908 | 0.8289 | 0.7919 |

The all-pathways figure is dragged down by the tuning ten (TP53, WNT, PIP3, cell
cycle), which are hard for both tools. **That is why the naive headline looks
like a loss**: all-pathways DeltaSignal 82.89% against MP-BioPath's published
83.61%. On held-out pathways the same comparison is:

| held-out | DeltaSignal | MP-BioPath | margin |
|---|---|---|---|
| each on its **own** networks (70 pathways) | **85.52%** | 85.28% | **+0.24pp** |
| **same** networks, MP-BioPath's (63 pathways) | **93.92%** | 84.90% | **+9.03pp** |

Other ground truths, all pathways:

| ground truth | result | macro-F1 |
|---|---|---|
| Reactome curator, valid-only | 19,289/23,022 = 83.79% | — |
| experimental evidence | 616/846 = 72.81% | 0.6521 |

### The tuning caveat, stated

Today's solver defaults (`specs/009`) were chosen by measuring on the wide
curator set, so the 70 "held-out" pathways were not held out from *that*
decision, even though they were from the paper's. Tested directly: the chosen
config wins **+31 on the tuning ten and +228 on the held-out rest**, so tuning
on the ten alone would have produced the identical choice. The bias is real,
small, and points the right way — but future decisions should be made on the
tuning set and reported on the rest.

For reference: MP-BioPath vs curator 83.61%, vs experimental 75.74%; curator vs
experimental (the human ceiling) 81.98%.

Per class against curators: NORMAL F1 0.877, UP 0.761, DOWN 0.737. **The entire
gap to MP-BioPath's headline is DOWN** (their F1 0.789); we are ahead on UP.
Closing DOWN alone would put macro-F1 at 0.8093 against their 0.8063.

### Why the number is lower on our networks, and why that is the point

Our networks are a harder substrate, deliberately:

- **Loops are retained.** MP-BioPath removed them by hand, curator-reviewed, a
  process that cannot be reproduced or maintained. We keep them and the solver
  handles them (SCC condensation).
- **EntitySets are decomposed** to member species, so a perturbation must reach
  the right variant rather than a generic placeholder.
- **Regenerated per Reactome release** from the graph database, not a one-off
  2019 MySQL extraction.

The same propagator scores 93.38% on their networks and 82.89% on ours. **The
deficit is network construction, not the model.**

---

## 3. The known, quantified gap

Of the 1,484 `no_path` failures (`bench/analysis/nopath_anatomy.py`):

| | cases | |
|---|---|---|
| route exists once positional variants of one stable id are collapsed | **1,094** | 73.7% |
| gene or readout absent from the network | 177 | 11.9% |
| MP-BioPath's network connects it, ours does not | 87 | 5.9% |
| genuinely unconnectable in either | 80 | 5.4% |

So `no_path` is mostly **not** missing biology — it is our node identity
splitting one curated entity across positional variants.

### The 87 "MP-BioPath connects it" cases are not lost curator edges

Traced one end to end (AKT1 → NR5A2 in *Regulation of beta-cell development*).
MP-BioPath routes AKT directly into the PDX1-synthesis reaction R-HSA-211272.
**That edge does not exist in Release97**: the reaction's participants are input
PDX1 gene, output PDX1, and regulators FOXA2, FOXO1, PAX6, MAFA. Their 2019
network carries a link the current database does not support, so this bucket is
evidence our networks are *more* faithful, not less.

The trace did expose a real gap of ours, though. The genuine route is
AKT ⊣ FOXO1: R-HSA-211164 "AKT phosphorylates FOXO1A" consumes FOXO1 and
produces p-FOXO1, and it is unphosphorylated FOXO1 that regulates PDX1
synthesis. We model AKT → reaction → p-FOXO1 correctly but nothing represents
FOXO1 being **depleted**, so no signal reaches the readout.

The generator's rule is explicit — `PI_STID = "R-ALL-29372"`, *"emit
catalyst→input depletion edges for PHOSPHATASE reactions only"*. A phosphatase
outputs Pi and gets depletion edges; a kinase outputs ADP and gets none, though
it consumes its substrate just as completely. That asymmetry has no biological
justification.

**Before anyone acts on it**: extending depletion to kinases is 2,267 catalysed
reactions in nucleoplasm alone (9,048 more in cytosol) against the current
rule's 3,446 — a 3-4x increase in the depletion footprint, not a small targeted
fix. Blanket substrate consumption has been tested twice and regressed hard
(-14pp; macro-F1 0.663 -> 0.529). Worth one careful A/B, not a confident
prediction.

**But collapsing them is measured net-negative, and this is why four attempts
have failed:**

    GAIN    1,094 real changes currently missed
    AT RISK 4,316 NORMAL cases that are correct ONLY because unreachable

A ratio of **1 : 3.9 against**. Any future attempt needs a rule that admits the
1,094 while rejecting most of the 4,316; a blanket bridge cannot win, and the
recorded failures (-73, -204, -4pp) are consistent with this ratio.

**Path length does not provide that rule** (`--separate`):

| | n | median hops | quartiles |
|---|---|---|---|
| recoverable | 1,094 | 12.0 | 9 / 18 |
| at risk | 4,316 | 11.0 | 8 / 14 |

Controlled within pathway across the 16 pathways with at least 5 of each, the
median difference is **−0.5 hops**, closer in 9 and further in 6. And the
intuitive rule — reconnect only short paths — makes the ratio *worse*: 1:5.0 at
4 hops, 1:7.9 at 6, 1:5.6 at 8, against 1:3.9 unrestricted.

So the silo is closed unless someone finds a discriminator that is not distance.
The next candidate worth an A/B is the kinase depletion asymmetry above, not
this.

---

## Generator changes: measured, and re-read on the held-out split

Two real generator defects, found while tracing the `no_path` bucket. Measured
first across all pathways, where all three arms looked negative — then re-read
on the held-out split, which reverses the reading.

| arm | all pathways | tuning ten | **held-out (70)** |
|---|---|---|---|
| phosphatase detection in all compartments | −9 | −7 | **−2** |
| depletion edges excluded from root detection | −54 | **−55** | **+1** |
| both together | −26 | −19 | −6 |

The root-detection row is **conditioned on the perturbation set**: a case
counts only where both arms resolved the same `n_gene_uuids`/`n_ko_uuids`, so
the two arms answered the same question. 112 of 23,908 cases (0.5%) fail that
and are excluded; unconditioned the held-out figure reads +2 rather than +1,
and the tuning figure is unchanged at −55. The conclusion does not move. The
other two rows predate the conditioning and their dumps were not retained, so
they are unconditioned — treat them as approximate, and note that a
cross-network comparison, where resolution differs almost everywhere, moved
from −1459 to −151 under the same correction.

**The root fix's entire −54 is the tuning ten, and almost all of it is TP53.**
On the 70 pathways outside the tuning set it is **+2 — neutral**. Rejecting a
correctness fix on that evidence would have been precisely the overfitting the
paper's protocol exists to prevent.

### The defects

*Phosphatase detection* keyed on `R-ALL-29372`, Pi **[cytosol]** alone, while
the rule's stated criterion is "the outputs include Pi". Nucleoplasmic (475
reactions), mitochondrial (343) and extracellular (186) phosphatases were
invisible. Derived by ChEBI now, 12 species at R97.

*Root detection* computed `sources - targets` **after** depletion edges were
appended, and counted them. A root is "produced by no reaction in this pathway";
a depletion edge is not production, it is our own inference. So a depletion edge
landing on a boundary complex silently deleted its curator-derived subunit
decomposition. The fix restores 46 assembly edges, among them
`p-MAPK1 → MAPK1 dimer` and the MAPK3/MAPK7 equivalents.

### What TP53 was doing

The root fix changes 312 predictions there, 239 with a decidable truth, and they
are direction flips between UP and DOWN landing on the truth 41.1% of the time —
worse than a coin toss. TP53 is dense enough that direction is unstable to a
small structural change, and it is a tuning pathway, so that instability should
not veto a fix that is neutral everywhere else.

### Two method errors, recorded

- "+28 for the compartment fix" was inferred by subtracting arms. Wrong —
  effects are not additive and the scored denominators differ (23,022 vs
  22,910). **Measure each arm.**
- The three arms were first reported as clear negatives from the all-pathways
  figure alone, before the held-out split was applied. **Split first.**

## AND mode: the specification vs the empirical winner

The AND default changed from `hill_log` to `hill_sat` on 2026-09-17.
`specs/010-and-multiplication-fidelity/research.md` has the full record; the
numbers that matter for a write-up:

**`hill_log` deviated from the specification below baseline, without bound.**
AND is specified as multiplication of fold-changes capped at 100. Measured
against the product: 0.5x2 = 0.0%, 2x2 = −0.9%, 0.5x0.5 = +0.9%,
0.25x0.25 = **+7.2%**, 0.1x0.1 = **+35.2%**, 0.001x1 = **+167.6%**. The error
is systematically UPWARD, so every down-regulated value was lifted toward
baseline. And `0 x anything` never reached 0 — it read 0.0007 to 0.0135 and
**rose with the co-input**, so an abundant partner partially rescued a
knockout.

**The cost of correctness is 16 held-out cases.** Release97, one freshly
regenerated catalog, conditioned pairing, 0 of 23,908 cases dropped. 81
pathways scored — 11 tuning (5,100 cases), 70 held out (18,808 cases) — from a
93-entry list; the rest lack a curator file or, in one case, do not exist in
Release97:

| arm | held-out net | held-out macro-F1 | all-pathway macro-F1 |
|---|---|---|---|
| `hill_log` (was default) | — | 0.8152 | 0.7846 |
| `hill_sat` (now default) | **−16** | 0.8138 (−0.0014) | 0.7835 |
| `hill_log_asym` (hybrid arm) | −12 | 0.8143 (−0.0010) | 0.7840 |

**The gain is depth-invariant magnitudes, and this is the publishable part.** A
chain of SINGLE-input reactions has nothing to combine, so it must be the
identity. It was not:

| source fold | `hill_log` d1 | d5 | d10 | `hill_sat`, any depth |
|---|---|---|---|---|
| 100x | 74.1 | 33.5 | **19.1** | **100** |
| 10x | 9.6 | 8.3 | 7.2 | **10** |
| 0.1x | 0.104 | 0.120 | 0.139 | **0.1** |

`hill_log` lost 81% of a 100x signal over ten hops, and the loss depended on
path length — median path length here is ~10 hops. Two readouts with identical
biology and different path lengths therefore received different predicted
folds. That is what made those continuous outputs unusable as a scale, which
was the stated reason `specs/002` adopted `hill_log` in the first place.
`hill_sat` is exact at every depth.

**The previous justification for `hill_log` was void.** `specs/002` and
CLAUDE.md recorded `hill_sat` as 90 cases worse. That was measured against a
`hill_sat` carrying two then-unknown defects: it could not represent a knockout
(a hardcoded `eps=1e-6` floor against baseline 0.01), and a catastrophic
cancellation in its 100x cap inverted the saturation on wide reactions, so a
maximally elevated wide reaction read as ~0 and
`Class_I_MHC_mediated_antigen_processing_presentation` did not solve at all.
Roughly 74 of the 90 were bugs. Corrected, the gap is 16.

**Method note.** The cancellation was found because the broken arm scored
**23,788 cases against the others' 23,908** — a 0.003 macro-F1 gap on a smaller
set is indistinguishable from a real result. Compare denominators before
metrics.

**Measured on the experimental axis too**, after an adversarial review found
the decision had been justified on curator ground truth alone. Same catalog,
conditioned, 849 of 849 comparable: `hill_log` macro-F1 0.6506 / 617 correct,
`hill_sat` 0.6487 / 615 — **net −2, both in TP53, which is 0 fixed and 2
broke and McNemar p = 0.50: no detectable difference.** The curator cost by
contrast is 9 fixed / 25 broke on the held-out split, p = 0.0090 — small but
real; on the tuning ten it is 2 / 5, p = 0.4531, not detectable. As percentages
(−0.085pp and −0.24pp) the two axes look similar; on the paired test one is a
signal and the other is a coin flip. So: correct AND arithmetic costs a
measurable amount of agreement with curator *reasoning* and has no measurable
effect on predicting experimental *outcomes*. Note that
every experimentally-evidenced pathway is inside the tuning ten, so there is no
held-out empirical split for either tool.

**Still unexplained.** Why downward compression helps classification at all.
It is a 16-case effect, not 90, but it is a consistent paired signal across two
independent formulations. The likeliest explanation — that lifting suppressed
values toward baseline masks false DOWN calls in over-coupled regions, making
it a symptom of over-coupling rather than a property of AND — is untested.

## Reproducing

```bash
# Both arms need the server's catalog mount to match DS_CATALOG_ROOT, or every
# observation silently matches nothing (guarded since PR #27/#29).
export PATHWAY_CATALOG=<catalog>
docker compose -f docker-compose.dev.yml up -d --force-recreate julia-api

DS_CATALOG_ROOT=<catalog> python bench/benchmark_vs_mpbiopath.py \
  --ground-truth curator --dump-cases cases.tsv --report report.tsv

# MP-BioPath's own networks as a DeltaSignal catalog
python bench/analysis/mpbiopath_network_adapter.py \
  --pathways ~/gitroot/mp-biopath-pathways/pathways \
  --catalog <catalog> --out mpb_catalog --cases cases.tsv

# Anatomy of the no_path bucket, and its gain/risk ratio
python bench/analysis/nopath_anatomy.py --cases cases.tsv \
  --ours <catalog> --mpb mpb_catalog [--exposure]
```
