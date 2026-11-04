# specs/134 research

## R1. korkut arms (2026-11-04)

Build `20261026-1349_8b6527e` (frozen), code `8b8d5fa`, `run_arm.sh --axis korkut`, arms `k134_forms` (control),
`k134_protein`, `k134_abund`. All 1,493 solves converged in each. The control reproduces specs/133's `c0` korkut
predictions case for case (same 11,626 cases, 0 class and 0 value differences). Scores: `score_134.py`
(`~/deltasignal-catalogs/analysis/134/score_*.json`), on cases scored in both arms.

| | cases | balanced acc. | accuracy | macro-F1 | change calls (wrong) | fixed / broken |
|---|---|---|---|---|---|---|
| forms (control) | 11,062 | 0.4133 | 0.6434 | 0.3368 | 3,922 (3,472) | |
| protein − forms | 11,062 | −0.0005 [−0.0011, 0] | −0.0016 | −0.0008 | 3,939 (3,489) | 22 / 39 |
| abundance, control on the same cases | 10,865 | 0.4161 | 0.6512 | 0.3416 | 3,766 (3,317) | |
| abundance | 10,865 | 0.4185 | 0.7003 | 0.3585 | 3,142 (2,735) | 637 / 103 |
| abundance − forms | | **+0.0024 [−0.0067, +0.0104]** | +0.0491 | +0.0169 | | |

On totals alone, abundance minus forms: balanced −0.032 [−0.054, −0.011], accuracy +0.073; on phospho every prediction
identical. Under abundance, 197 control-scored cases leave scoring: 132 become trivial (every abundance copy is a
pinned drug target, e.g. total AKT under an AKT inhibitor) and 65 have no abundance copy (`no_copy`).

**The gate fails**: the interval includes 0.

## R2. Where the remaining total change calls come from

Each abundance readout split by its copies (starting copy = no in-edge except depletion):

| copies | cases | change calls | right |
|---|---|---|---|
| all starting copies | 1,710 | 0 | 0 |
| starting and made copies | 4,848 | 1,472 | 30 |
| made copies only | 742 | 409 | 10 |

Every change call left comes from a copy the pathway makes (translation of an mRNA, or a protein produced by a
reaction with no input of that protein), and 40 of 1,881 are right. Traced: in *Interleukin-4 and Interleukin-13
signaling*, PI3K p85 exists only as the translation product of a STAT3-driven mRNA; under a STAT3 or JAK inhibitor it
reads 0 (DOWN); measured, no change. p21 in *TP53 Regulates Transcription of Cell Cycle Genes* reads 11.7× under
nutlin; measured, no change. The model holds no stock of a made protein apart from what it is making now, so stopping
the synthesis empties it; the cells keep the protein they had over the assay's time.

The correct DOWN calls the control made on totals (76) were mostly bound or modified forms dropping, right for a reason
the antibody does not measure; most are lost.

## R3. Predictions

- Q1 holds: phospho identical in all three arms.
- Q2 refuted: total change calls 2,505 → 1,881 (not under 500); accuracy on totals 0.637 → 0.710 (not above 0.90).
- Q3 holds: `protein` moves 61 of 7,497 total cases (0.8%).
- Q4 refuted: the gate fails; balanced accuracy +0.0024 with an interval across 0, while accuracy gains +0.049.
- Q5 holds: every remaining total change call comes from synthesis copies, and 40 of 1,881 are right.
- Q6: CPPA training, see R4.
