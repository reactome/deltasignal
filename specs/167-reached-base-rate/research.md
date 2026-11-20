# specs/167 research: the measurement (2026-11-20)

Reader `bench/analysis/167/base_rate.py` (fixture-tested, committed at `a0fec811` before this read). Output:
`~/deltasignal-catalogs/analysis/167/base_rate.json`. Records: build `20261113-1311_baf608b_r98fix`, solver
`12b5db3` (curator, experimental), `cff961d` (korkut), `552ba89` (decryptM).

| axis | cases | reached | changed, of reached | changed, of unreached | call rate | precision (any change) | precision (called direction) | lift [95%] |
|---|---|---|---|---|---|---|---|---|
| curator held-out | 18,481 | 6,312 | 0.883 | 0.074 | 0.954 | 0.904 | 0.888 | +0.021 [+0.009, +0.035] |
| curator tuning | 5,030 | 2,832 | 0.931 | 0.037 | 0.927 | 0.929 | 0.872 | −0.003 [−0.008, +0.002] |
| experimental | 845 | 707 | 0.938 | 0.703 | 0.972 | 0.936 | 0.876 | −0.002 [−0.004, −0.001] |
| korkut, own markers excluded | 9,476 | 1,307 | 0.142 | 0.049 | 0.839 | 0.165 | 0.109 | +0.023 [+0.014, +0.034] |
| korkut, own markers included | 10,469 | 1,999 | 0.225 | 0.066 | 0.894 | 0.248 | 0.213 | +0.024 [+0.015, +0.035] |
| decryptM | 628 | 323 | 0.325 | 0.246 | 0.972 | 0.328 | 0.325 | +0.003 [−0.007, +0.012] |

Lift = precision (any change) − the share of reached readouts that changed. The "called direction" column is the
specs/165 precision (0.888, 0.876, 0.109).

Reached readouts the model declined to call:

| axis | declined | of those, changed |
|---|---|---|
| curator held-out | 288 | 0.444 |
| curator tuning | 206 | 0.966 |
| experimental | 20 | 1.000 |
| korkut, own markers excluded | 210 | 0.024 |
| decryptM | 9 | 0.222 |

## The reading

1. **On every axis, a change call is about as good as the base rate among the readouts a route reaches.** The
   model calls a change on 84–97% of what it reaches, so its precision is the share of reached readouts that
   changed, plus at most 0.024. The condition fixed in the spec holds for the lift on every axis (largest
   +0.024). The call-rate half of it (≥ 0.85) holds on four of five: korkut with own markers excluded is 0.839.
2. **What the call adds to reachability is the same size on a gene axis and a drug axis**: +0.021 on the curator
   held-out pathways and +0.023 on korkut, both intervals excluding 0. On curator tuning, experimental and
   decryptM it adds nothing.
3. **The precision gap between the gene axes and the drug axes is the base rate.** 88–94% of reached readouts
   are labelled a change on the curator and experimental axes; 14% on korkut and 33% on decryptM. The model does
   not choose better among reached readouts for gene perturbations than for drugs.
4. **The cutoff is informative on korkut and costly on the gene axes.** The 210 reached korkut readouts the
   model declines changed 2.4% of the time against 14.2% for all reached ones. On curator tuning 96.6% of
   declined readouts are labelled a change, and on experimental all 20. That is why every contraction of
   propagation gained on the drug axes and lost on the gene axes (specs/128, specs/164): where nine reached
   readouts in ten are labelled a change, declining any of them loses.
5. **The gene axes are selected readouts.** Readouts the network cannot reach changed 70.3% of the time on the
   experimental axis, which is 90% changes overall: it is a set of reported effects. On the curator axes the
   label is close to a statement about reach (0.88–0.93 where a route exists, 0.04–0.07 where none does):
   curators derived the expectation by reading the same pathway. korkut measures every readout on its array
   whatever was expected (unreached 0.049). decryptM sits between (0.246; its cases start from regulated sites).

## Correction to specs/165, reading 2

specs/165's gates and figures stand. Its reading that "the low precision is specific to drug readouts", and the
conclusion drawn from it and from specs/164 that drug responses propagate more weakly than gene perturbations, do
not follow from this record. The recorded axes confound the perturbation type with how the readouts were chosen:
every gene-perturbation axis has readouts chosen because a change was expected or reported, and every axis that
measures an unselected panel is a drug axis. Precision by axis is recorded correctly; "by perturbation type" is
not established. `docs/API.md` is reworded accordingly.

What the record supports: the model carries reach and direction. A reached readout is about 3 times as likely
to have changed as an unreached one on korkut (0.142 against 0.049), and a call on a real change has the right
direction 66% of the time on korkut (86% with own markers included), 94% on experimental and 98% on curator
held-out. It does not carry which reached readouts change by more than the cutoff.

## What would separate perturbation type from readout selection (FR-004)

An unselected readout panel under genetic perturbation. The two readings predict different base rates among
reached readouts there: near 0.9 if genetic perturbations do carry through the curated routes, near korkut's
0.14–0.23 if the gene axes' figure comes from how their readouts were chosen.

| dataset | perturbation | readouts | status |
|---|---|---|---|
| Lun 2019 | 649 kinases and phosphatases over-expressed | 36 phospho markers, CyTOF | TEST, sealed; opening it is Adam's decision alone. Best kept as the final test of whichever reading holds |
| genome-wide Perturb-seq, K562 (Replogle 2022; specs/096 shortlist #4) | CRISPRi knockdown of nearly every expressed gene; each perturbation's own knockdown is measured | transcriptome | not downloaded |
| LINCS L1000 CRISPR (specs/096 shortlist #3, its recommended genetic check) | ~7,500 knockouts | 978 measured transcripts | not downloaded |

Perturb-seq fits the open question best. Engagement is measured for every perturbation, which the drug axes
lacked (specs/138, specs/161). Cells are sampled days after the knockdown, close to the steady state the model
solves. It is large enough to split by perturbation into FIT, SELECT and TEST before any readout is read.

**How much of it the served networks could score** (`bench/analysis/167/transcript_feasibility.py`, a structural
count on the served build with Release 98 gene-expression regulations from Neo4j; output
`~/deltasignal-catalogs/analysis/167/transcript_feasibility.json`; no external data):

- Reactome 98 curates 1,446 gene-expression regulations on 1,066 events; 846 genes are the product of one.
- 371 of those target genes are nodes of the 92 served networks, and 3,698 genes have a root copy to perturb.
- That gives 93,853 (perturbed gene, target gene) pairs in the same pathway, 18,479 of them reached, in 79
  pathways. korkut has 9,476 scored cases and 1,307 reached.
- Largest: Transcriptional Regulation by TP53 (3,480 reached pairs), Class I MHC antigen presentation (2,770),
  Mitotic G1 (1,952).

Caveats: the count is before any filter on K562 expression or knockdown strength, and uses a plain root rule
(no in-edge, depletion ignored) rather than the bench's pin protocol. A transcript reads the gene's expression
event, which the networks carry as the translated protein copy, so the readout rule needs registering (the
specs/134 abundance rule is the starting point).

Nothing was downloaded. Under specs/096 a dataset is approved by Adam before any download.
