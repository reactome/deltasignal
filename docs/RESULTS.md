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

## 2. On our own networks

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
