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

**Corrected in R5 (amendment 1): the split below is wrong.** It classed copies by producing reactions only, as the
code did, and missed assembly edges; re-counted, most remaining change calls come from starting complexes.

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
- Q5 withdrawn (R5): most remaining total change calls come from assembled starting complexes, not synthesis; 40 of
  1,881 are right.
- Q6 holds on both parts, see R4.

## R4. CPPA training portion (2026-11-04)

Build `20261027-1346_8b6527e_top114`, code `289cae5`, `DS_TOTAL_READOUT=abundance`, 3 servers; control the frozen CPPA
run (`8276c28/cppa_frozen`, the default readout, whose code this spec does not change). Re-scored here, the control
reproduces the stored frozen training score exactly. Convergence is identical in both runs. Scored on the training
portion only, with `score_121.py` (primary) and a paired secondary (balanced accuracy at |delta-p| ≥ log2(1.5), on
(sample, antibody) pairs both runs predict, interval from resampling drugs). The test portions stay sealed. Numbers are
kept with the private analysis (`~/deltasignal-catalogs/analysis/134/`) under the CPPA terms.

- Primary (per-drug correlation): moves by much less than 0.01 and much less than its interval's half-width. The
  no-harm condition holds.
- Secondary (balanced accuracy): rises, and its interval lies above 0. Phospho predictions are identical, as on
  korkut.

So Q6 holds. The CPPA training portion was not looked at when the rule was chosen, so this is the independent check the
spec names, and it agrees with the rule's direction. It does not change the registered gate, which is on korkut and
fails. **Not adopted** under the registration; adoption is Adam's decision. A real test is one scoring on the CPPA
test portions or Lun 2019.

## R5. R2 re-counted (amendment 1, before its arms ran)

Traced first: total AKT under an AKT inhibitor in *Signaling by VEGF* (`R-HSA-194138`) reads 0 under R1's `abundance`.
Its 12 abundance copies have no in-edge in the reaction graph; 3 are AKT1/2/3 (pinned by the drug), 9 are starting
complexes of AKT with PIP3 and PDPK1, which `boundary_edges.csv` builds from AKT through `assembly` edges (specs/030).
The solver reports no twin, carrier or pool edges in that solve; the complexes fall with AKT through assembly.

Each R1 `abundance` total case (7,300 scored), by its copies (starting = no in-edge in the logic network except
depletion) and by whether the readout's gene is a target of the case's drugs:

| copies | drug target? | cases | change calls | right |
|---|---|---|---|---|
| starting only | no | 5,686 | 906 | 24 |
| starting only | yes | 372 | 372 | 5 |
| starting and made | no | 551 | 181 | 1 |
| starting and made | yes | 36 | 36 | 0 |
| made only | no | 643 | 384 | 10 |
| made only | yes | 12 | 2 | 0 |

1,278 of 1,881 change calls come from readouts whose copies all start the graph: complexes assembled from the protein
and its partners, which move with the partners or, for a drug's own target, with the knockout pin. Synthesis copies
account for about 600. R2's "every change call left comes from a copy the pathway makes" and Q5's "holds" are
withdrawn; the trace of p85 and p21 stands as an example of the synthesis class.
