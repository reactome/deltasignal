# specs/179: first-hop engagement — does a protein moving its direct substrate tell us when to trust the rest?

**Status:** registered 2026-11-27, before any engagement call was made or counted.
**Requested by:** Adam, 2026-11-27 (option 1 of the follow-ups to specs/176: "I think 1", then "I like it").

## Why

specs/176 found that the largest factor in Lun 2019 is whether the over-expressed protein is active at all in these
cells (17 of 288 half-A proteins change nothing, 53 change ten or more readouts). The model assumes every one is
fully active. This asks whether one measurement from the same well can supply that: **did the protein move a site
Reactome says it phosphorylates directly?** It is the over-expression counterpart of drug engagement on korkut
(specs/138, specs/161).

The protein's own site cannot serve: it rises with the protein's amount (20 of 24 measured self pairs read UP).

## Definitions (fixed here)

- **Direct pair.** `bench/analysis/179/direct_substrates.tsv`, built from the Release 98 top-level build
  (`20261126-2322_d6d026b_top98`) without reading any Lun value: the protein is in a catalyst of a reaction whose output
  carries a phosphate at one of the readout's registered sites that no input or catalyst carries. 135 pairs, 86
  proteins.
- **Engaged.** A protein with at least one measured direct pair is *engaged* if any of them reads UP by the specs/173
  call, *not engaged* if none does. Proteins with no measured direct pair are *unknown* and are not in the test.
- **Test pairs.** The scored pairs of engaged and not-engaged proteins, **leaving out every direct pair** (they define
  engagement) and self pairs. Predictions: the specs/175 file (the Release 97 build, abundance rule), with and without
  the specs/177 tyrosine-kinase rule. The Release 98 predictions are still running; the test is repeated on them when
  they finish.
- Both halves of Lun 2019 are FIT. The test is run on each half separately, as a replication, and pooled.

## Registered tests

**E1. Does the model do better for engaged proteins?** Precision of the model's change calls (right direction), engaged
minus not engaged, with a protein-bootstrap 95% interval (10,000, seed 179). Also reported: route enrichment (change
rate with a route over without) inside each group.

**E2. Is it the first hop, or just general activity?** specs/176 already showed that a protein changing many readouts
predicts more changes. The same comparison as E1 with a placebo engagement: the protein raises at least one readout it
has **no route to**. If the placebo separates precision as well as the direct substrate does, the first hop carries
nothing specific. Reported: both differences, and E1 within proteins whose placebo engagement is the same.

**E3. A gated model (the gate).** The model's change calls are withheld for not-engaged proteins. As a false-alarm
fix (specs/161's rule):
- fixed > broken, sign test p < 0.05;
- recall of measured changes over the test pairs falls by at most 0.01;
- **on each half separately.** Both halves must pass.

Reported: balanced accuracy of the gated model against the model and against reachability gated the same way.

## Reading

- E3 passes on both halves and E2 shows the direct substrate separates better than the placebo: a measured first hop
  tells when the rest of the prediction can be trusted. Proposed to Adam as an engagement report for over-expression
  data and as a use for soft observations (specs/156); nothing is adopted here.
- E3 passes but the placebo does as well: it is general activity (specs/176 T6), not the first hop. Recorded.
- E3 fails: engagement at the first hop does not predict the rest. Recorded.

**Limits, stated in advance.** Engagement is known for at most 86 proteins. It needs a measurement from the experiment,
so it cannot improve a blind prediction. A direct substrate that is a tyrosine site under a tyrosine kinase is also
covered by the specs/177 rule; results are reported with tyrosine-specific kinases left out as well.
