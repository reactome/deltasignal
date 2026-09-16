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
to +2.95pp, and the median per-pathway difference is +0.008. The defensible
claim is **"comparable on signalling, decisively better on transcriptional
regulation"** — not a uniform +8.8pp.

### Failure profile on the same networks

| category | our networks | MP-BioPath's networks |
|---|---|---|
| `no_path` | 1,484 | 328 |
| `false_positive_change` | 1,441 | 597 |
| `propagator_missed` | 808 | 417 |

---

## 2. On our own networks

| ground truth | result | macro-F1 |
|---|---|---|
| Reactome curator, end-to-end | 19,818/23,908 = **82.89%** | 0.7919 |
| Reactome curator, valid-only | 19,289/23,022 = **83.79%** | — |
| experimental evidence | 616/846 = **72.81%** | 0.6521 |

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

**But collapsing them is measured net-negative, and this is why four attempts
have failed:**

    GAIN    1,094 real changes currently missed
    AT RISK 4,316 NORMAL cases that are correct ONLY because unreachable

A ratio of **1 : 3.9 against**. Any future attempt needs a rule that admits the
1,094 while rejecting most of the 4,316; a blanket bridge cannot win, and the
recorded failures (-73, -204, -4pp) are consistent with this ratio.

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
