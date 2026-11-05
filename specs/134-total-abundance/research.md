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

## R6. Amendment 1 arms

korkut: code `f03d0f9`, `run_arm.sh --axis korkut --bench DS_TOTAL_READOUT=abundance`, arm `k134_abund_a1`, frozen
build; all 1,493 solves converged. Against R1's control `k134_forms`, on cases scored in both (`score_134.py`):

| | cases | balanced acc. | accuracy | macro-F1 | change calls (wrong) | fixed / broken |
|---|---|---|---|---|---|---|
| control | 10,493 | 0.4239 | 0.6738 | 0.3545 | 3,394 (2,950) | |
| abundance, amendment 1 | 10,493 | 0.4411 | 0.7961 | 0.4027 | 1,956 (1,576) | 1,503 / 219 |
| difference | | **+0.0172 [+0.0007, +0.0324]** | +0.122 | +0.048 | | |
| RAF set aside: difference | 9,905 | **+0.0212 [+0.0038, +0.0373]** | +0.138 | +0.056 | | 1,491 / 129 |

R1 with RAF set aside, for comparison: +0.0029 [−0.0069, +0.0113].

**The registered gate passes**, with and without RAF. On totals alone, balanced accuracy falls (−0.030
[−0.065, +0.006]; RAF set aside −0.022 [−0.056, +0.015]) while accuracy rises 0.671 → 0.856: the correct DOWN calls
on totals are lost, and the pooled gain is the "not" class recovered on totals while phospho recall is unchanged.
Phospho predictions are identical.

569 control-scored cases leave scoring (504 `trivial`, 65 `no_copy`), 528 of them a drug's own target whose every
abundance copy the drug's knockout pins (total AKT under an AKT inhibitor); the control was wrong on 522 of them. They
are left out of both arms, as registered, so they do not enter the difference; an inhibitor modelled as removal cannot
predict its target's abundance, and here the bench now says so instead of reading the knockout.

Remaining total change calls, classed as in R5 but with boundary edges included (starting = no in-edge at all):

| copies | drug target? | cases | change calls | right |
|---|---|---|---|---|
| starting only | no | 5,572 | 0 | 0 |
| starting and made | no / yes | 578 / 36 | 250 / 36 | 3 / 0 |
| made only | no / yes | 730 / 12 | 407 / 2 | 10 / 0 |

Every remaining total change call (695, 13 right) is on a readout with a copy the pathway makes from something else
(translation; a protein produced by a reaction with no input holding it). R2's account holds once the rule is
implemented as written.

Predictions: Q7 holds. Q8 holds for (a), synthesis copies; (b) is refuted, only 38 of 695 are a drug's own target,
because most such cases leave scoring. Q9 is refuted: the gate passes. Q10: CPPA training, R7.

## R7. Amendment 1 on the CPPA training portion (2026-11-04)

Build `20261027-1346_8b6527e_top114`, code `7cb140c` (the rule as in `f03d0f9`), `DS_TOTAL_READOUT=abundance`, 2
servers; control and scoring as in R4. Convergence is identical to the frozen run. Numbers are private, as in R4.

- Primary (per-drug correlation): rises; the paired difference by drug has an interval across 0. The no-harm
  condition holds. A drug whose only predicted changes were totals of its target's binding partners (read DOWN
  because the complexes they sit in lost the target) now predicts no change, as the rule intends: a partner's own
  stock does not move.
- Secondary (balanced accuracy): rises more than in R4, and its interval lies above 0. Phospho predictions are
  identical.

Q10 holds. Under the decision above, `abundance` is now the korkut and CPPA bench default.

## R8. Amendment 2 on korkut (2026-11-04)

Code `e4809b8`, `run_arm.sh --axis korkut --bench DS_TOTAL_READOUT=abundance`, arm `k134_abund_a2`, frozen build; all
1,493 solves converged. Scored with `score_134.py`.

Case by case against `k134_abund_a1`: phospho identical; nothing moves outside the four readouts named in Q12. In RAF/MAP
kinase cascade, cRAFV and fibronectin (42 each) and SRC and bRAF (30 each) move DOWN → NOT, measured NOT in all 144; the
other 24 (12 SRC, 12 bRAF) leave scoring as `trivial`, a drug's own target whose remaining copies its knockout pins.
That is the 168 change calls of Q12, and none elsewhere.

| | cases | balanced acc. | accuracy | change calls (wrong) |
|---|---|---|---|---|
| amendment 1 (`k134_abund_a1`) | 10,469 | 0.4418 | 0.7980 | 1,932 (1,552) |
| amendment 2 (`k134_abund_a2`) | 10,469 | 0.4469 | 0.8117 | 1,788 (1,408) |
| difference | | +0.0051 [+0.0039, +0.0063] | +0.014 | 144 fixed / 0 broken |

Against R1's control `k134_forms` on the 10,469 cases both score: balanced +0.0224 [+0.0060, +0.0376], accuracy
0.675 → 0.812. RAF set aside (9,905 cases): +0.0212 [+0.0038, +0.0373], identical to amendment 1, as every moved case is
in RAF.

Q11 holds. Q12 holds exactly. Q13 holds (accuracy +0.014; balanced rises with RAF kept, unchanged with RAF set aside).
Q14: CPPA training, R9.

## R9. Amendment 2 on the CPPA training portion (2026-11-05)

Build `20261027-1346_8b6527e_top114`, code `e4809b8`, 2 servers; control the R7 run (`7cb140c/cppa134_abundance`),
scoring as in R4. Convergence identical. Numbers private, as in R4.

- Primary (per-drug correlation): does not fall; the paired difference by drug lies above 0.
- Secondary (balanced accuracy): unchanged overall (interval across 0). Phospho identical.

Q14 holds. Amendment 2 meets its registered criteria on both korkut (R8) and CPPA training; adoption is Adam's.

## R10. Amendment 3 on korkut (2026-11-05)

Code `a490cdc`, `run_arm.sh --axis korkut --bench DS_TOTAL_READOUT=abundance_a3`, arm `k134_abund_a3`, frozen build;
all 1,493 solves converged. Against `k134_abund_a2` (R8), scored with `score_134.py`: **no case moves**, not a class
and not a predicted value. b-Catenin in Signaling by WNT now reads 3 copies (was 2), but in that network all of them
rest at baseline under every korkut drug. Q15 holds.

## R11. Amendment 3 on the CPPA training portion (2026-11-05)

`scripts/run_cppa.sh` at `a490cdc` on the three registered pathways, spliced into the R9 run for every other pathway;
convergence identical to R9 there (the readout rule changes no solve). Numbers and traces private, as in R4.

- Only the registered (pathway, readout) pairs move: Q16's first part holds.
- Primary: no harm. Q16 holds.
- Secondary: the interval spans 0; the private notes record which calls moved and how they scored.

Amendment 3 meets its registered criteria on both axes. Adoption is Adam's.
