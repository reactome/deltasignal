# A held-out empirical axis: feasibility

**Status**: FEASIBILITY SIZED. **The transcript framing below is REJECTED** —
see "Protein, not transcript" for the design that replaced it. No claim made.

## Why this and not more accuracy work

`docs/RESULTS.md` states the problem against us. On MP-BioPath's own networks
we beat their published curator result decisively (93.29% against 84.58%,
re-measured on the current solver 2026-09-21) and we are **level on
experimental evidence** (77.67% against 79.21%). On our own networks we are
*behind* their published experimental figure (70.20% against 75.74%).

A model that reproduces curator *reasoning* far better than its predecessor
while being no better at predicting what actually happened in cells is exactly
what a reviewer will suspect: something that tracks the representation rather
than the biology. On the current evidence we cannot refute that, and **every
accuracy point won this month was on the curator axis**, which strengthens the
critique rather than answering it.

The reason the axis cannot arbitrate is structural, not analytical:

- Experimental truth exists for **ten pathways, and all ten are the tuning
  ten**. There is no held-out empirical test for either tool.
- 849 cases, 156 readouts, one pathway (PIP3) carrying the differences.
- A regeneration with **no change at all** moves this axis by 13 cases
  (measured 2026-09-21), so it cannot resolve anything smaller than that.

MP-BioPath's experimental truth is a **hand-curated 0/1/2/-999 matrix per
pathway** (`experimental_results/*.tsv`, rows = key outputs, columns =
`GENE_0` knockdown / `GENE_2` over-expression). That is why it stopped at ten
pathways: extending it their way means manual curation. Their `expression/`
(E-MTAB-2836) and `cellline_metadata/` (CCLE drug response) directories serve
a different purpose and are not the source of that matrix.

## The alternative: derive the axis from public perturbation data

Sized by `bench/analysis/empirical_axis_feasibility.py` (tested in
`test_empirical_axis_feasibility.py`), on catalog `cat_prod`, solver
`5032771`:

| | cases | pathways | triples |
|---|---|---|---|
| curator axis | 24,100 | 82 | 12,056 |
| experimental axis today | 849 | 10 | 556 |
| curator truth, **no** experimental truth | 19,132 | **72** | 9,572 |

The naive read of that last row — "9,572 triples available" — is wrong, because
most readouts cannot be measured by a transcriptional assay at all. The funnel,
restricted to the 72 pathways with no experimental truth:

| readout kind | cases | share |
|---|---|---|
| `simple_complex` | 8,658 | 45.3% |
| `simple_entity` (incl. dissociation sinks) | 7,088 | 37.0% |
| unmapped | 2,670 | 14.0% |
| `reaction` | 716 | 3.7% |

Of the single-entity readouts, only those whose name is a **bare gene symbol**
survive. Rejected deliberately: phospho-forms (`p-S317_S345-CHEK1`) and
ubiquitin-forms, which are state changes and not abundance changes; lipid
modifications (`4xPalmC-CD36`); proteolytic fragments (`ACIN1(1-1093)`), which
are not distinct transcripts; and unnamed numeric ids.

**Addressable: 4,048 cases / 2,024 triples / 136 readout genes / 162 perturbed
genes / 19 pathways — every pathway outside the tuning ten.**

That is roughly 4.8x today's case count and, more importantly, it would be the
**first held-out empirical test available to either tool**.

## The gating question, not yet answered

Does public perturbation data cover it?

- knockdown and over-expression signatures for **162** perturbed genes;
- measured transcript levels for **136** readout genes.

Candidate sources are LINCS L1000 consensus signatures (~12,300 measured
genes, several thousand perturbagens) and DepMap CRISPR gene effect. Checking
this requires fetching third-party data, which has not been done.

**The caveat no dataset removes**: mRNA abundance is a proxy for protein
activity, and our predictions are about activity. A derived axis is therefore
weaker evidence than the hand-curated matrix per case — its value is breadth
and held-out status, not per-case fidelity. Any result from it must be reported
as "transcript-proxy agreement in held-out pathways", never as "experimental
accuracy" alongside the 849.

## Pre-registration (to be committed BEFORE any arm runs)

Not yet written — there is nothing to pre-register until external coverage is
known. If coverage clears, say, 60% of the 2,024 triples, the arm is: score
the existing solver against the derived axis on the 19 held-out pathways, and
report agreement, the confusion matrix, and the per-pathway spread, with **no
tuning against it whatsoever**. Its first job is to be a held-out witness, not
a target.

## Defect found and fixed while sizing this

The first version of the sizing script resolved a readout's `node_kind`
catalog-wide, first match wins. **2,785 of 15,616 stable ids appear as
`simple_entity` in one pathway and `dissociation_sink` in another**, so the
answer depended on directory iteration order: two runs over identical data
reported 1,691 and 1,812 triples. Now resolved per pathway, and pinned by two
tests that fail if the catalog-wide lookup is restored.


---

## Protein, not transcript

**Direction taken 2026-09-21: only protein-level data is acceptable here, and
the test has to be both statistically significant and defensible as an
argument.** That voids the transcript plan above. Reactome models **proteins,
complexes and their modification states**. An mRNA readout inserts a proxy step
between what we predict and what is measured, and a reviewer would say exactly
that — it would be the *weakest* possible evidence against the "tracks the
representation, not the biology" critique, not an answer to it.

**The inversion this creates.** The readouts the transcript framing discarded
are the ones where our claim is most specific and most directly checkable.
Re-classified by what a protein assay could measure, across all 82 pathways:

| readout class | cases | triples | readouts | perturbed genes | pathways |
|---|---|---|---|---|---|
| complex / other modified form | 14,724 | 7,368 | 558 | 490 | 79 |
| unmodified protein | 6,192 | 3,096 | 197 | 194 | 25 |
| proteolytic fragment | 1,642 | 821 | 80 | 125 | 11 |
| **phospho, named site** | **836** | **418** | **36** | **134** | **16** |
| phospho, site not parsed | 706 | 353 | 19 | 178 | 20 |

### The argument to make, and why it is not dismissible

Take the named-site phospho readouts: `p-S317_S345-CHEK1`, `p-S133-CREB1`,
`p-S166_S188-MDM2`, `p-S216-CDC25C`, `p-T308-AKT1`-style entries — a named
residue on a named protein.

- Reactome asserts that a specific kinase phosphorylates that specific residue.
- DeltaSignal predicts a direction for that site under perturbation of the
  responsible kinase.
- Site-level phosphoproteomics, or a phospho-specific antibody, measures
  **that exact site**.

The modelled entity and the measured entity are *the same object*. There is no
proxy step to attack. That is the difference between "our numbers correlate
with something" and "we predicted this phosphosite would fall and it fell."

Only 418 triples — but 418 three-class calls has ample power against any
sensible null (marginal-distribution or all-NORMAL baseline), so **significant
and arguable are not in tension here**; the small n is what makes each case
strong rather than weak.

### The design that answers BOTH open questions at once

The other standing question is whether the *magnitude* of our predicted
changes means anything. Dose-resolved phosphoproteomics under kinase inhibition
(the decryptM-style design: many drugs, several doses, site-level readout)
supplies both in one dataset:

1. **direction agreement** at named sites → the validity test above;
2. **dose-response** → the magnitude test, because a graded perturbation should
   produce a graded, monotone predicted response, and depth along a cascade
   should attenuate it in a specific way.

That is the clever version: one experimental modality, protein-level and
site-specific, that tests the two things we cannot currently defend.

### Candidate sources, coverage NOT yet checked

| source | perturbation | readout | fit |
|---|---|---|---|
| dose-resolved phospho under kinase inhibition | drug, several doses | phosphosite | best: direction **and** magnitude |
| LINCS P100 | drug / genetic | ~90 phosphopeptides | clean perturbation, narrow readout panel |
| RPPA panels (cell lines / tumours) | genomic alteration | ~200-450 antibodies, phospho-rich | large n, observational |
| deep tumour proteogenomics | driver event | proteome + phosphoproteome | deep, observational, confounded |

The gating question is now: **how many of the 36 named sites (and the 197
unmodified-protein readouts) appear in a dataset where the responsible
upstream gene is also perturbed?** That needs third-party data.

### Two things to decide before any arm runs

1. **Complexes are 61% of the curator axis and have no direct protein readout.**
   Arguing from subunit abundance imports a modelling assumption ("a complex
   tracks its limiting subunit") which is itself one of the things under test.
   Either exclude complexes and say so, or test that assumption separately —
   do not quietly fold it in.
2. **Observational designs are confounded.** A tumour with a loss-of-function
   event carries co-occurring alterations and lineage effects. A perturbation
   design (drug or knockout, measured before and after) does not. Prefer the
   perturbation design even at much smaller n, because the argument is the
   point.

---

## Correctness and utility are separate claims

**Decided 2026-09-21.** Tumour-level work (predicting which pathways are
affected, clustered against clinical outcomes such as survival) is a genuine
result about **utility**. It cannot establish **correctness**, because cancer
data is messy enough that a correct prediction and a confounded one look
identical. The two claims must not argue for each other:

| claim | evidence | why it fits |
|---|---|---|
| the tool is **correct** | named phospho-sites under perturbation of the responsible kinase | mechanistically specific, protein-level, no proxy step |
| the tool is **useful** | TCGA pathway-level calls clustered against survival | clinically meaningful, and messiness is tolerable because the claim is weaker |

The failure mode is using the second to argue the first. Tumour data cannot
establish correctness — co-occurring alterations, lineage effects and selection
mean a correct prediction and a confounded one look the same. Survival
association is a genuine result about *utility* and should be presented as
exactly that, downstream of correctness rather than as a substitute for it.

## The target list (committed: `phosphosite-targets.tsv`)

836 rows — every case whose readout is a named residue on a named protein,
with the perturbed gene, our prediction, the curator call and whether they
agree. **34 sites, 418 triples, 134 perturbed genes, 16 pathways.**

The list is favourable, because the canonical kinase-substrate pairs are in it
and each has a well-characterised inhibitor that phosphoproteomics studies
routinely profile:

| site | perturbed upstream | cases | curator agreement |
|---|---|---|---|
| `p-S317_S345-CHEK1` | ATM, ATR, CHEK1, BRCA1/2 … (48 genes) | 138 | 119/138 |
| `p-S133-CREB1` | AKT1/2, MAP2K1, MAPK1, EGFR … (33) | 68 | 52/68 |
| `p-S216-CDC25C` | ATM, ATR, CHEK2, PKN1-3 … (18) | 36 | 32/36 |
| `p-T369_S640_S964_S975-RBL1` | CDK4, CCND1, CCNE1, MYC … (17) | 34 | **16/34** |
| `p-S183_T246-AKT1S1` (PRAS40) | AKT1/2, PIK3CA/B, PDPK1 … (16) | 32 | 30/32 |
| `p-S939_T1462-TSC2` | AKT1/2, PIK3CA/B, PDPK1 … (16) | 32 | 30/32 |
| `p-S166_S188-MDM2` | AKT1/2, PIK3CA/B, PDPK1 … (16) | 32 | 30/32 |
| `p-S99-BAD` | AKT1/2, PIK3CA/B, PDPK1 … (16) | 32 | 30/32 |
| `p-T210-PLK1` | AURKA, PLK1, BORA, TPX2 (6) | 24 | **24/24** |
| `p-T308_S473-AKT1` | PIK3CA, AKT1, KDR, SRC (7) | 14 | 10/14 |
| `p-S123-CDC25A` | ATM, ATR, CHEK2, WEE1, TP53 (7) | 14 | **14/14** |
| `p-Y705-STAT3` | EGFR, ERBB2/3/4, PTK6, DOK1 (7) | 14 | 12/14 |

`p-S345-CHEK1` is the standard readout in every ATR-inhibitor paper.
`p-T246-PRAS40`, `p-T308/S473-AKT`, `p-S99-BAD` and `p-TSC2` are standard
PI3K/AKT phospho-panel members. `p-T210-PLK1` is the canonical Aurora-A→PLK1
readout. These are not obscure assertions we would be asking anyone to take on
trust.

The weak sites are informative too, not embarrassing: `p-T369…-RBL1` at 16/34
and `p-S95-PHLDA1` at 6/12 are where we already disagree with the curator, so
a phospho measurement there discriminates between "we are wrong" and "the
curator expectation is wrong", which no curator-only comparison can do.

### The one assumption to state out loud

Our benchmark perturbation is gene **knockdown or over-expression**;
phosphoproteomics perturbation is usually a **small-molecule inhibitor**.
Treating kinase inhibition as loss-of-function is defensible *for downstream
phosphorylation specifically* — an inhibited kinase does not phosphorylate its
substrate — but it is an assumption, not an identity: inhibitors have
off-targets, and inhibition does not remove scaffolding or non-catalytic
functions. State it, and prefer genetic perturbation where the data exists.

### Next step

Check, for these 34 sites and 134 perturbed genes, which appear in a dataset
where the upstream gene is perturbed and the site is measured. That requires
third-party data and has not been done.

---

## A matching public dataset exists: decryptM 2.0

Searched 2026-09-21. The design proposed above — dose-resolved
phosphoproteomics under kinase inhibition — has a public dataset built for
almost exactly it.

**decryptM 2.0**: **17 million peptidoform dose-response curves for 133
clinical kinase inhibitors** across 5 cell lines (A204, A431, MESSA, SKES1,
SKLMS1), phosphoproteome and full proteome, with fitted dose-response curves
per peptidoform. Zenodo `10.5281/zenodo.17533475`; the earlier decryptM
(Science 2023, `10.1126/science.ade3925`) is the 31-drug predecessor.
Per-site browsing without a bulk download is possible through ProteomicsDB's
decryptM explorer, which makes a cheap pilot feasible before committing to the
full pipeline.

Why it fits the argument we need:

- **site-level** phosphorylation, which is the entity Reactome models;
- **graded doses**, so the same dataset tests direction *and* magnitude;
- **clinical kinase inhibitors**, so the perturbed node is a named kinase that
  our networks contain;
- dose-response curves are already fitted, so the comparison is against a
  curve parameter rather than a single noisy measurement.

### What must be checked before this is a plan

1. **The drug-to-target table.** 133 inhibitors, but the list is in the paper's
   Table S1 and was not enumerated on the record page. The gating question is
   how many of our **134 perturbed genes** are the primary target of one of
   them. The classes we need are well represented among clinical kinase
   inhibitors — EGFR, PI3K, AKT, MEK, ATR, CHK1, WEE1, PLK1, Aurora — but
   "well represented" is a prior, not a measurement.
2. **Cell-line coverage is the real risk.** Five lines, and they are sarcoma
   and epidermoid (rhabdoid, epidermoid carcinoma, uterine sarcoma, Ewing,
   leiomyosarcoma). A site is only measurable where the protein is expressed
   and phosphorylated in those lines. DNA-damage and PI3K/AKT readouts are
   plausible; several of our 16 pathways almost certainly are not represented.
   **This will shrink 418 triples substantially, and the shrinkage must be
   reported as coverage, not hidden by quoting only the surviving cases.**
3. **Inhibition is not knockdown.** Already stated above; with dose-resolved
   data it is partially testable, since a dose series should look like graded
   loss of function rather than a switch.
4. **Direction mapping.** Our benchmark direction codes are knockdown (0) and
   over-expression (2). An inhibitor series only supplies the knockdown arm, so
   **half of each triple pair is unavailable** and the over-expression
   predictions stay untested by this route.

### Pre-registration sketch (to be written properly before any arm)

Score only sites present in both, report coverage first and accuracy second,
use the knockdown arm only, and make **no tuning decision** against this axis —
its first outing must be as a held-out witness. The magnitude test is separate
and stricter: monotonicity of our predicted fold against dose, and attenuation
with path depth, both pre-specified.

Sources: [decryptM 2.0 dataset](https://zenodo.org/records/17533475) ·
[Science 2023](https://www.science.org/doi/10.1126/science.ade3925) ·
[ProteomicsDB decryptM explorer](https://www.proteomicsdb.org/decryptm)

---

## Can we say the magnitudes mean something?

An open goal for the project: state something defensible about the *size* of a
predicted change, not only its direction.

### First, the blocker: there is barely any magnitude to validate

Measured on `cat_prod`, all 24,100 curator cases, baseline = 1.0 on the UI
scale:

| predicted value | cases | share |
|---|---|---|
| exactly baseline (no change) | 15,527 | **64.4%** |
| at the floor (0) | 2,660 | 11.0% |
| at the ceiling (100x) | 1,528 | 6.3% |
| interior, below baseline | 1,602 | 6.6% |
| interior, above baseline | 2,783 | 11.5% |

**82% of predictions are baseline, zero, or full saturation**, and even the
18.2% interior is pushed to the top (median 10x, 2,326 of 4,385 above 10x).
Cases carrying a genuinely graded prediction between 0.5x and 10x are
**1,036 of 24,100 — 4.3%.**

So the model is close to three-valued in practice. That is a consequence of the
AND semantics we deliberately chose: `hill_sat` multiplies fold-changes capped
at 100, which saturates quickly along any cascade, and a knockout propagates to
exact zero. specs/010 adopted that for depth-invariance, and the cost is
dynamic range.

**Therefore a quantitative magnitude claim ("a predicted 5x means 5x") is not
available and should not be attempted.** Any correlation computed over these
values would be dominated by the rails.

### What IS supportable today, from data already in hand

Magnitude is a **calibrated confidence signal**. Direction accuracy against the
curator, by predicted fold:

| predicted fold | cases | direction correct |
|---|---|---|
| knockout to 0 | 2,660 | **88.6%** |
| 0–0.5x | 1,023 | 75.1% |
| 0.5–0.9x | 428 | 57.7% |
| **1.1–2x** | **204** | **37.3%** |
| 2–10x | 226 | 63.3% |
| 10–100x | 2,326 | **88.3%** |
| railed at 100x | 1,528 | 78.5% |
| ~baseline (0.9–1.1x) | 15,705 | 85.8% |

The shape is a **U**: confident at both extremes, worst in the middle. A
predicted 1.1–2x change is right 37.3% of the time on a three-class problem —
barely above chance, and *worse than saying nothing*. And full saturation
(78.5%) is worse than 10–100x (88.3%), which matches the earlier finding that
railing destroys information.

That is a real claim about magnitude meaning something, it is defensible
without any external data, and it is **useful for shipping**: a small predicted
change should be suppressed or greyed out in the Reactome view rather than
displayed as a call. The publishable form is a reliability curve plus a
precision/coverage trade-off — "discarding the smallest predicted changes
raises precision from X to Y at Z% coverage".

### The test that makes the claim falsifiable

Replace every predicted magnitude with the mean magnitude of its own predicted
class and re-run the analysis. If the calibration survives that, the magnitude
carries **no** information beyond the direction call and the claim is empty.
This must be run and reported; without it "magnitude means something" is
unfalsifiable.

### The stronger claims, and what each needs

| claim | test | needs | verdict |
|---|---|---|---|
| ordinal | Spearman of predicted magnitude rank against measured effect-size rank, change cases only | decryptM | achievable; needs no unit calibration |
| calibrated | binned predicted fold against mean measured log2 fold, monotone | decryptM | interpretable, stronger than a correlation |
| **dynamics** | **EC50 ordering along a cascade: sites nearer the inhibited kinase should saturate at lower dose** | decryptM dose curves | **the best one — see below** |
| quantitative | predicted log-fold ≈ measured log-fold, slope ≈ 1 | decryptM | do not attempt; the rails forbid it |

**Why the EC50 ordering test is the one to aim for.** decryptM fits a
dose-response curve per phosphosite, so it yields a *potency* per site, not
just an effect size. Reactome's topology says which sites are one step from the
inhibited kinase and which are three. A graded loss of kinase activity should
reach proximal sites at lower dose than distal ones, so our networks predict an
**ordering of EC50s** along each cascade. Testing rank agreement of that
ordering:

- needs no calibration of our units, only their order;
- is a statement about *dynamics*, which is what magnitude is supposed to
  encode;
- has no reason to come out right in a model that merely tracks the curator's
  representation, which makes it the sharpest available discriminator;
- and it is testable on the knockdown arm alone, which is all an inhibitor
  series provides.

### Necessary internal checks, to run before any external work

Cheap, local, and if either fails there is no point going further:

1. **Dose monotonicity.** Perturb one input at 2x, 5x, 20x, 80x and confirm
   every readout's predicted fold is monotone in the input. Non-monotonicity
   would mean the magnitude is an artefact of the solver, not a response.
2. **Depth behaviour.** specs/010 chose `hill_sat` *for* depth-invariant
   magnitudes, so predicted fold should **not** decay with path length. That is
   a design claim currently untested at catalog scale, and it is in tension
   with the EC50-ordering test above, which expects proximal and distal sites
   to differ. Resolving that tension is a prerequisite, not a detail: if our
   magnitudes are genuinely depth-invariant, we predict *no* EC50 ordering and
   the test becomes a falsification of the design rather than a confirmation
   of it.


---

## Method recovered from a deleted script

A feasibility-scan script under `bench/analysis/` was removed because its
docstring named a collaborator, their lab and their unpublished data. It is not
named here either, since the filename carried the disease area. The *method* in it is
dataset-independent, was not recorded anywhere else, and is the same shape as
the sizing above, so it is written down here rather than lost.

**The method.** Given a gene set of interest and a pathway catalog, classify
every member that appears in a pathway as:

- **input-only** — feeds reactions but is not produced in-pathway (a root);
- **output-only** — produced by a reaction but feeds nothing (a terminal);
- **both** — an intermediate.

Then find pathways that contain input-only members *and* output-only members
with a **directed path** between them. Those are the testable cases: pin the
inputs to their measured values, predict the outputs, and compare the predicted
direction against the measured one. It is a way of turning any differential
gene or protein list into a set of "pin these, predict those" cases without
needing a designed perturbation.

**No code was actually lost, contrary to a first draft of this note.** The
script carried a forward reachability BFS, and `bench/analysis/_common.py` does
only hold backward helpers (`bfs_upstream`, `shortest_path_back`) — but the
repo has upwards of a dozen forward implementations elsewhere, including
`check_silo_bug._bfs_reach`, `nopath_anatomy.reaches`,
`version_skew.reachable` and `benchmark_vs_mpbiopath.reachable_from`. The
Cypher gene-to-stable-id query was byte-identical to `build_gene_cache.py` and
the loader duplicated `_common.load_one`. Only the *method* above was unique,
which is why it is recorded here.

Nothing else was lost: the Cypher gene-to-stable-id query was byte-identical to
`build_gene_cache.py`, and the loader duplicated `_common.load_one`.

---

## Pre-registration: does magnitude carry case-level signal? (committed before running)

**The claim under test.** Direction accuracy is U-shaped in predicted fold:
88.6% for a predicted knockout to zero, 88.3% at 10–100x, 37.3% at 1.1–2x. The
proposed reading is that magnitude is a *calibrated confidence signal* — a
small predicted change is less trustworthy than a large one.

**Why the obvious falsification test is useless.** "Replace each magnitude with
its class mean and show calibration collapses" cannot fail: it removes all
within-class variation by construction, so the calibration always collapses.
That test would have been reported as a success whatever the data said.

**The real threat is composition.** The U-shape could arise entirely from which
cases end up where:
- *predicted class* — extreme folds may be mostly one direction, and that
  direction may simply be easier;
- *pathway* — some pathways saturate more AND score higher. This is the
  between-pathway confound that killed nine candidate levers in this project,
  and saturation, path length and readout in-degree have all previously FAILED
  to discriminate errors within a pathway. So the prior is against this claim.

**Test.** Strength = |log10(predicted fold)|, fold on the UI scale where 1 is
baseline, with a predicted zero floored at 1e-6. Cases predicted NORMAL are
excluded — the claim is about the size of a predicted *change*. Within each
(pathway, predicted class) stratum containing at least one correct and one
incorrect case, compute the AUC of strength for predicting correctness, and
combine strata weighted by their number of correct-incorrect pairs (a
stratified Mann-Whitney). Null: permute strength within strata, 2,000 times.

**Primary split: held-out.** Tuning reported alongside, not used to decide.

**Pre-registered predictions.**
- **P1** (reproduces the signal): unstratified AUC > 0.5.
- **P2** (the test that matters): stratified AUC > 0.5 with permutation
  p < 0.01 on held-out.
- **Decision.** If P2 holds, the magnitude claim survives composition and may be
  stated as case-level calibration. If P2 fails, the U-shape is composition and
  **no magnitude claim is made**; this is recorded as a negative result, and the
  "suppress small predicted changes" recommendation is withdrawn as unsupported
  at the case level.

**Data.** Catalog build `20260925-1039_d4f4f64`, solver `ae84de9`, file
`results/ae84de9/curator_cases.tsv` (24,100 cases).

### Result: P2 fails as pre-registered. Case-level calibration is NOT ESTABLISHED.

Run 2026-09-25 against build `20260925-1039_d4f4f64`, solver `ae84de9`, by
`bench/analysis/magnitude_calibration.py` (7 tests, including one that pins the
statistic removing pure composition). 2,000 within-stratum permutations.

**Held-out (decides)** — 5,891 changed predictions:

| stratification | AUC | strata | pairs | permutation p |
|---|---|---|---|---|
| none | 0.5601 | 1 | 4,508,248 | — |
| within predicted class | 0.5475 | 2 | 2,258,774 | — |
| **within (pathway, class)** | **0.5134** | 85 | 72,564 | **0.1154** |

**Tuning (reported only)** — 2,509 changed predictions:

| stratification | AUC | strata | pairs | permutation p |
|---|---|---|---|---|
| none | 0.5741 | 1 | 1,078,858 | — |
| within predicted class | 0.5717 | 2 | 541,808 | — |
| within (pathway, class) | 0.6791 | 21 | 98,770 | 0.0005 |

**P1 holds; P2 fails as pre-registered**, and withdrawing the "suppress small
predicted changes" recommendation under the pre-registered rule is correct.
**What P2's failure does not show is that magnitude carries no case-level
signal.** An earlier revision of this section said it did. Adversarial review
showed that overstated the test, in three places.

**1. The pre-registered weighting is the one standard choice that hides the
signal.** Pair weighting (correct x incorrect cases per stratum) gave 30% of all
weight to the two DSB Repair strata, whose within-stratum AUCs are 0.461 and
0.530 — flat. The same 85 held-out strata under other standard combinations:

| weighting (exploratory, post hoc) | AUC | p |
|---|---|---|
| pairs — pre-registered | 0.5134 | 0.11 |
| van Elteren, the textbook stratified Wilcoxon | 0.5475 | 1.7e-7 |
| equal weight per stratum | 0.5806 | 6e-8 |
| pathway sign test, weighting-free | 29 of 42 pathways positive | 0.0098 |

A logistic regression with (pathway x class) fixed effects gives a positive
strength slope (likelihood-ratio statistic 152). So "almost all of it is
composition" was also wrong: under van Elteren, adding pathway to the
stratification leaves the AUC at 0.5475, unchanged. The drop to 0.513 came from
the change of weighting, not from removing composition.

**2. "Power is not the explanation" was false.** Simulating on the real stratum
structure, 80% power at p < 0.01 needs a pair-weighted AUC of about 0.536. The
test could not have resolved an effect of the size the other analyses find.

**3. The honest uncertainty is wide.** Strata are heavily overdispersed (sum of
z-squared 530 on 79 df), so every within-stratum permutation p here, the
pre-registered one included, is anti-conservative. Resampling whole pathways
gives van Elteren AUC **[0.494, 0.618]** with P(AUC <= 0.5) = 0.06.

**What can be said.** There is moderate, not decisive, evidence of a small
within-pathway effect: cases predicted to change only slightly are less often
right than strong ones in the same pathway and direction (exploratory: held-out
cases below 5x are 60.2% correct against 88.6% for strong cases in the same
strata). But the practical gain is small. Suppressing predictions below 5x
raises held-out precision from 84.65% to 85.80% (+1.15pp), and below 2x by
+0.27pp. The recommendation is not reinstated: it was not supported by its own
test, and it buys about a point.

**Why a regeneration is not a replication.** Rerunning on a freshly built
catalog would not settle this: held-out predictions differ by only ~15 cases
across relabellings of identical content, so it would re-find the same result
from the same cases. A real replication needs new ground truth, which is what
the phosphosite axis is for.

**The tuning result is one pathway, not an overfitting signature.** An earlier
revision read the tuning AUC (0.679, p = 0.0005) as "the shape overfitting
leaves". It is not. The Transcriptional Regulation by TP53 strata carry more
than 100% of the excess over 0.5 — the TP53 UP stratum alone has AUC 0.875 and
45% of the pairs — and dropping TP53 takes tuning to 0.459. Within a single
perturbation, tuning is 0.508 (p = 0.39). So it is between-perturbation
composition inside TP53, the pathway already known to be dominated by the MDM2
loop basin. The p = 0.0005 is also just the permutation floor (1/2001).

**Corrections to the motivating numbers.** The U-shape figures that motivated
this test (88.6% / 88.3% / 37.3%) came from the lost `cat_prod`, not the tested
build. On the tested build the held-out 1.1–2x bin holds only 24 cases in 7
pathways, against 175 on tuning, so the 37.3% was largely a tuning phenomenon.
The tuning set is also 11 pathways in `TUNING_PATHWAYS`, not ten.

**Design limitation, recorded.** A monotone score (|log10 fold|) applied to a
U-shape that turns down at the very top — for up-predictions of 10x or more,
10–100x beats railed 100x within stratum (Mantel-Haenszel OR 2.17) — dilutes
the signal. Combined with pair weighting, the pre-registered test was
structurally tilted toward a null. It stands as the record; the lesson is to
pre-register van Elteren or a pathway sign test next time.

**Status.** Not refuted, not established. Unlike the saturation, path-length and
in-degree levers, this is an unresolved result rather than a dead one, and it
should not be counted among them.

## Pre-registration: dose-response (committed before running)

**The question**, as asked: does moving an input further from baseline move the
outputs further from baseline? The calibration test above could not answer it,
because every benchmark perturbation is the same size — knockdown to UI 0 or
overexpression to UI 80 — so input magnitude never varies there.

**Design.** Rerun the identical benchmark pipeline at graded strengths, changing
only `DS_PERTURB_UI_DOWN` / `DS_PERTURB_UI_UP`, so gene resolution, readout
aggregation and every other step are unchanged:

| run | knockdown to | overexpression to |
|---|---|---|
| 1 | 0.5 | 2 |
| 2 | 0.2 | 5 |
| 3 | 0.05 | 20 |
| 4 (existing) | 0 | 80 |

Unit: each (pathway, gene, direction, readout) where both the gene and the
readout resolve in the network. Build `20260925-1039_d4f4f64`, solver = the
commit these runs are made at.

**M1 — monotonicity (the pass/fail one).** Across the four strengths in one
direction, a readout's output should move consistently: non-decreasing or
non-increasing (an inhibited readout legitimately falls as its input rises).
Among readouts that move at all, **P1: at least 95% are monotone.** A readout
whose output *reverses* as its input strengthens would be a solver defect, not a
modelling choice, and anything below 95% means the magnitudes cannot be read as
a response at all.

**M2 — graded or switched (descriptive, with a stated reading).** Among moving
readouts, the fraction already at a rail (output 0 or at/near the 100 cap) at
the mildest input (2x up, 0.5x down), and the fraction whose output still
changes between the mildest and strongest input. **Reading rule, fixed now:** if
more than half of moving readouts are already railed at a 2x input, magnitude is
effectively binary for them and cannot carry dose information.

**M3 — transfer (descriptive).** For moving, unrailed readouts, the slope of
|log output fold| against |log input fold|. The `hill_sat` design claims
depth-invariance, so a single-input chain should pass a fold through unchanged
(slope 1); multi-input reactions are expected to damp it (slope below 1).

**No accuracy claim is made from these runs.** The ground truth is defined at
full strength, so accuracy at a partial strength is not comparable to anything.

### Result: P1 passes as registered; on non-trivial readouts it straddles the gate — NOT ESTABLISHED

Build `20260925-1039_d4f4f64`, solver `3328b5f` (runs 1–3); run 4 is the
`ae84de9` production scoring of the same build, **copied, not rerun**. The
solver is identical between the two commits (only the health endpoint changed),
so the knob being inert at its default is established by reading the code (its
one use is the pin value), not by an experiment.
`bench/analysis/dose_response.py`.

Two rounds of adversarial review (PR #69) replaced two earlier versions of this
section. The first said "P1 passes, 97.0%, M3 median exactly 1.00,
depth-invariance holding on real networks". The second said "P1 FAILS, 94.2%,
depth-invariance withdrawn". Both overstated. What follows is what the data
supports.

**Readouts that copy the input.** 3,423 of the 8,604 moving readouts (2,858
held-out) equal the pinned input at every step, to a relative 1e-4. They are the
perturbed node itself (a key-output uuid set can include the gene's own pinned
node) or an exact pass-through. The case table records uuid counts, not uuids,
and resolving a `key_output` dbId needs Neo4j, so the two are not separated
here. Either way they are monotone by construction and cannot test M1. The
pre-registration did not exclude them. **Excluding them is itself a post-hoc
choice**, made after seeing the data, exactly like the print tolerance below.
So neither column is "the" result.

| held-out | as registered | non-identity | non-identity, print tolerance |
|---|---|---|---|
| moving readouts | 6,063 | 3,205 | 3,205 |
| M1 monotone | **97.0%** | **94.4%** | **95.4%** |
| further from baseline, same side | — | 91.4% | 94.0% |
| M2 railed at the mildest input | — | 29.0% | |
| output changes mildest → strongest | — | 70.3% | |

Tuning, non-identity: M1 94.2% (95.1% at print tolerance), further from
baseline 77.8% (79.7%), M2 30.1%.

"Print tolerance" means one 1e-6 print unit on a step, and 1e-5 on the |log10
fold| step. `pred_ui` is written to six decimals, so the registered 1e-9 is
below the file's own resolution.

**M1 — not established.**
- As registered, P1 passes (97.0%).
- Removing the trivial readouts gives 94.4% (95.4% at print tolerance), which
  straddles the 95% gate.
- The binomial standard error at n = 3,205 is about 0.4pp, and the cases are
  correlated within a perturbation, so the gate lies inside the noise.
- The data says a stronger input moves a readout monotonically in about 94–95%
  of non-trivial cases. It does not say which side of 95% that falls.

**M1 is also weaker than the question.** `monotone` checks the raw direction,
not distance from baseline. On "further from baseline, same side of 1", the
result is 91.4–94.0% held-out but only 77.8–79.7% on tuning. The tuning
pathways (TP53, PIP3, cell cycle) are the loopy ones.

**M2 — graded, not switched.** 29% of non-identity moving readouts are at a
rail at a 2x / 0.5x input, under the 50% reading rule. 70% still change between
the mildest and the strongest input.

**M3 — descriptive only, and biased low.** The slope is fit on finite-input
steps only. The knockdown to 0 has no finite log, and flooring it at 1e-6
dominated the fit in the first version.
- Held-out median 0.57, IQR 0.07–1.00, over 1,395 never-railed non-identity
  readouts; tuning median 0.12.
- That sample is 44% of the non-identity movers, and it excludes by design the
  strongest transmitters (e.g. an OE readout reaching the 100 cap). It
  describes the unrailed subset, not the typical readout.
- The identity readouts may include genuine multi-hop pass-throughs, which
  would be evidence *for* depth-invariance.
- Depth-invariance on real networks is therefore **not established**; it is
  not withdrawn as false.
- 10.5% of held-out slopes (21.8% tuning) are negative: those readouts move
  *less* as the input strengthens.
- Why transfer is damped where it is damped is not traced.

**Reversals.**
- 244 non-identity readouts reverse by more than one print unit. Only **92
  reverse by ≥ 1e-3 UI**: 52 held-out, in 12 pathways, led by RUNX2 (34), TP53
  (17) and PIP3 (16).
- In the 10 acyclic pathways with non-identity movers there are 0 reversals
  among 205 readouts, and all 205 are further from baseline. All 92 large
  reversals are in cyclic pathways.
- 205 is a small base, and the zero is partly by construction: an acyclic
  component is solved exactly in one pass, so iterative noise cannot occur.
- The largest reversals are real sign changes. WNT5A knockdown in Signaling by
  WNT holds readouts at 100x at 0.5 and 0.2, then collapses them to ~1e-5 at
  0.05 and 0. CHEK2 knockdown in Cell Cycle Checkpoints spikes once instead
  (0.53 → 100 → 6.3 → 0).
- Both are consistent with the specs/014 loop basin flip, but neither is
  traced.

**What this licenses.** Magnitudes are graded rather than binary. A stronger
input usually (about 94–95% of non-trivial readouts) moves a readout
monotonically. A small, identifiable set of cyclic cases reverse sign. Neither
the pre-registered dose-response bar nor depth-invariance is established.
Case-level calibration also remains not established (above).

**Analysis fixes, each pinned by a test:**
- Identity readouts are recognised with a relative tolerance and excluded from
  M1-non-identity, M2 and M3.
- `away` applies its tolerance to both of its checks.
- Slopes are fit on finite steps only.
- `load_ladder` warns on dropped cases; none was dropped here (23,268 in every
  run).
- `benchmark_selective_joint.py` ignores `DS_PERTURB_UI_*` and now says so.

### Traced: the WNT5A reversal starts in one acyclic reaction, not in the loop

Case: Signaling by WNT (`R-HSA-195721`), WNT5A knockdown, readout 3322393.
Readout 100, 100, 1.3e-5, 1.6e-5 at KD 0.5, 0.2, 0.05, 0. It reproduces
exactly on re-solve, and every solve reports converged. Build
`20260925-1039_d4f4f64`, solver `a22d752`.

1. **The readout is not in a loop.** A 168-node strongly connected component
   upstream of it flips entirely between KD 0.2 and 0.05 (all 168 members), and
   559 further nodes outside every loop flip with it (not checked to be all
   downstream). The flip costs iterations: 90 sweeps
   in the high basin, 266 in the low one.
2. **The loop only amplifies.** Its perturbed entries are 48 edges from
   `WNT:FZD:LRP5/6`, which reads **2, 2, 0.5, 0**, a non-monotone input
   computed outside every loop.
3. **The source is one reaction, *WNT binds to FZD and LRP5/6*.** Its activator
   is the WNT ligand (x = 0.5, 0.2, 0.05, 0). Its two negative regulators,
   `WIF1:WNT` and `WNT3A:sFRP`, are sequestration complexes that *contain the
   ligand*, so each also reads x. Under `divide`, each contributes 1/x, and the
   product is clamped at the de-repression ceiling of 10
   (`reaction_model.jl:1859`):

       output = x * min(x^-2, 10)   ->   2, 2, 0.5, 0      (matches exactly)

   A knockdown therefore reads as **up** (1/x) until the ceiling binds at
   x = 1/sqrt(10) ≈ 0.32. Below that it reads 10x, which falls under baseline
   once x < 0.1.
4. **This is the specs/012 self-contained-inhibitor structure.**
   `self_contained_inhibitor_pairs` flags both regulators of this reaction.
   specs/012 recorded that *one* such inhibitor cancels the signal exactly
   (x · 1/x = 1). With *two*, the substrate is double-counted into an
   inversion, 1/x. Across the catalog, 580 reactions have exactly one
   self-contained inhibitor and **94 have two or more**, in 10 pathways: PIP3 61,
   WNT 16, and one more (`R-HSA-177929`) with 7. PIP3 carries 16 of the 92 large
   reversals. That is suggestive, not traced.

**Consequences for the result above.**
- The "large reversals are loop basin flips" reading is only half right. In
  this case the loop is the amplifier and the cause is acyclic operator
  composition.
- `x * min(x^-2, 10)` gives 2, 2, 0.5, 0, which is non-increasing, so `monotone`
  accepts it. In an acyclic pathway this readout would count as monotone under
  M1, even though it crosses baseline. Only `away` rejects it. The "0 reversals
  among 205 acyclic readouts" therefore cannot exclude this mechanism in acyclic
  pathways.

**A second, separate observation (not assessed).** The 37 uuids pinned for
"WNT5A" include generic complexes whose WNT component is a set (`WLS:WNT`,
`WIF1:WNT`). Pinning them knocks down every WNT ligand (WNT1, 3A, 4, 8A, 8B,
9A), not only WNT5A. That is the benchmark's rule of pinning every entity
containing the gene. Whether it over-perturbs set-containing complexes is its
own question.

**Not done: no fix is claimed.** specs/012 measured *deleting* these edges as
−61 held-out. The bounded alternative, consistent with what has paid off on
these networks, would be for the self-contained inhibitors of one reaction to
contribute **once**: combine them by min rather than product, so x · 1/x = 1
(the specs/012 single-inhibitor behaviour) and the response stays monotone. It
needs its own pre-registration and a held-out A/B on both axes.
