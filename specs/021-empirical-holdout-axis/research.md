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
