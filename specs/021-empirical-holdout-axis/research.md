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

## Protein, not transcript (Adam, 2026-09-21)

> "I don't really trust aligning anything but protein data with these tests as
> that is what we are really modelling in reactome mainly. we need to be clever
> to come up with something that not only comes out with significant results
> but also is arguable."

Correct, and it voids the transcript plan above. Reactome models **proteins,
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

Adam's other standing question is whether the *magnitude* of our predicted
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
