# specs/177: the tyrosine-kinase rule, and the one scoring of Lun 2019 half B

**Status:** registered 2026-11-26, before half B was opened.
**Decision:** Adam, 2026-11-26. Asked whether to register the rule with a clean definition and score the candidate
protocol once on the sealed half B: "yes go ahead with what ever you think is right".

## Why

specs/176 found that the largest pattern in Lun 2019 half A is not a pathway: an over-expressed tyrosine kinase raises
tyrosine phosphosites across the panel, with or without a curated route. It was found on half A with a crude class.
This spec fixes a clean class in advance and tests the rule once on proteins the analysis has never seen.

## The rule

For an over-expressed protein and a phospho readout:

- **the protein is a tyrosine-specific kinase**: Reactome's reference entity for it carries the UniProt keyword
  "Tyrosine-protein kinase" and not "Serine/threonine-protein kinase" (so the dual-specificity MAP2Ks are out).
  Table: `bench/analysis/177/tyrosine_kinases.tsv`, Release 97, 60 of the 636 Lun proteins;
- **the readout has a tyrosine at a registered site** (canonical sequence): p-STAT3, p-STAT5, p-SRC, p-FAK,
  p-BTK/ITK, p-PLCγ2, p-ERK1/2, p-p38, p-JNK. Table: `bench/analysis/177/readout_sites.tsv`;
- then the call is **UP**, wherever the model itself calls no change. The model's own change calls stand.

Both tables hold identifiers only and are committed with this registration. Nothing in the rule is fitted.

**The candidate protocol** is the frozen model's Lun predictions under the abundance rule (specs/175, file sha256
`21e287c3…`, made before this spec) plus this rule. The route-length and PIP3 candidates of specs/174 are not in it
(specs/176 T5: they cost more recall than they gain).

## Half A (development, already seen)

With the clean class: tyrosine kinase × tyrosine site is UP in 72% of 306 pairs (9–14% elsewhere); 27 of 34 kinases
raise at least half of their tyrosine readouts. Candidate balanced accuracy 0.435, model 0.384, reachability 0.380;
candidate − reachability +0.056 [+0.030, +0.083]. Added calls: 120, 89 right. This is where the rule came from, so
it is not evidence.

## Half B: scored once

`bench/analysis/177/score_177.py B`, on `lun2019_b_perbin.csv` (sha256 `6ae9507d…`, checked before reading), with
the specs/173 measured call unchanged (same bins, same control reference, same thresholds). Population: the P ∩ B
pairs. The positive control (a protein's own site reads UP, at least as in specs/127 step 7) must pass or scoring stops.

**Primary gate.** Three-class balanced accuracy, candidate − reachability, protein bootstrap (10,000, seed 177): the
95% interval's lower end above 0.

**Secondary gates.**
- candidate − model: lower end above 0;
- the added calls: right (measured UP) outnumber wrong, sign test p < 0.05.

**Reported, as specs/173 registered for the model on half B:** model − reachability and model − always-NOT; route
enrichment and lift with intervals (enrichment interval above 1.5 = fit for development, including 1 = not).

**Reported, as replication of specs/176 (no gate):** the class table (share UP for tyrosine kinase × tyrosine site
against the other three cells); the AUROC table (protein × readout, route, model, candidate); the share of basal change
calls that repeat in the EGF well; the arms with the measured call at 2 δ and 3 δ.

## Reading

- Primary gate passes: the rule holds on unseen proteins. It is proposed to Adam as the Lun protocol and as a served
  statement for over-expression queries ("a tyrosine kinase over-expressed at this level raises tyrosine sites beyond
  its curated routes"); adoption is his decision.
- Primary gate fails: the rule was a half-A pattern. Recorded as such; nothing is proposed.

Either way half B is not scored again: after this run `lun2019_b` is registered FIT (seen once), as decryptM and
korkut were after theirs.
