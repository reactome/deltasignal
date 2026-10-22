# specs/089 research: arms ulu (gated) and ulu88

Variant build `20261022-1541_d31f95e_ulu` (LNG `d31f95e`, `LNG_UB_LIGASE_UNIT=1`), solver
`569bef6` (main defaults plus specs only), against `results/e4e7192/dflt086b` on canonical
`20261021-2218_f99c4ab`. Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade
--pair-on-valid` (589 curator and 4 experimental cases invalid in both, unpaired). Case lists:
`diff_cases.py`.

| arm | results |
|---|---|
| ulu | `20261022-1541_d31f95e_ulu/results/569bef6/ulu` |
| ulu88 (`DS_CARRIER_LOOP_OPEN=fed_single_ub`) | `20261022-1541_d31f95e_ulu/results/569bef6/ulu88` |

| split | ulu net (fixed / broken) | p | ulu88 net | p |
|---|---|---|---|---|
| curator held-out | **−2** (1 / 3) | 0.63 | −2 (7 / 9) | 0.80 |
| curator tuning | +22 (26 / 4) | 0.0001 | +22 (27 / 5) | 0.0001 |
| experimental (RAF set aside; RAF does not move) | +3 (5 / 2) | 0.45 | +3 (6 / 3) | 0.51 |

**The registered gate fails on both arms**: held-out −2 < 0. No pathway is below −10
(worst −1). Not adopted as registered.

## Predictions

1. `ulu`: at least 2 of the 4 CCNE1 experimental cases fixed: **holds**, all 4 (Mitotic G1
   overexpression at CCNA1 and CCNA2, S phase knockout at R-HSA-187568 and R-HSA-187920).
2. `ulu`: EGFR and MET move by at most 5 curator cases: **holds**, 0 (none changes).
3. `ulu88`: Mitotic G1 does not lose specs/088's 6 / 2, and both MDM2 cases stay fixed:
   **partly**. MDM2 holds (both experimental cases, plus 4 TP53 curator cases in Cell Cycle
   Checkpoints). CDKN2A and PPP2R1A no longer lose, but CCND1 and CDK4 knockout at CCNA1
   (68891) still do (0.85× → 0.93×, curator and experimental). Mitotic G1 nets +2 / 0.
4. No prediction for held-out: −2.

## Where it moves (`ulu`, curator)

- **Gains:** TP53 +10 (AKT1, AKT2, MDM4 knockout, BRCA1), Mitotic G1 +9 (CCNE1, CCND1, CDK4,
  CDKN2A, PPP2R1A, MYC, MAX, RBL1), WNT +3, S phase +1 (CCNE1, CDKN1B knockout).
- **Losses:** Cell Cycle Checkpoints −1, S phase −2 inside its +1, Hedgehog −1, TGF-beta −1,
  RUNX2 0 (1 / 1).

## The losses, case by case

- **Held-out TGF-beta, SMAD7 knockout at p-SMAD2,3:SMAD4** (curators UP): 4.56× → 1.00×. The
  removed edges are exactly the class targeted: TGFBR:SMAD7:SMURF/NEDD4L draining free SMURF1
  and SMURF2. On canonical, SMAD7 knockout reached the readout only through those drains. The
  curated route (SMAD7 holding the receptor) does not carry the knockout. Like specs/078's CBL
  class, the answer was right through the wrong edge. Not traced further.
- **Held-out Hedgehog, SHH knockout at GLI3R** (curators NORMAL): 1.00× → 2.74× UP. The removed
  edge is PTCH1:SMURF draining free SMURF1/2. GLI3R forms when Hedgehog is absent, so UP is
  arguably the biology, and the curated NORMAL is questionable. Not counted either way here.
- **Held-out RUNX2, PPM1D knockout:** a threshold pair, +1 / −1 (0.76× → 0.86×, 1.13× → 1.17×).
- **Tuning Cell Cycle Checkpoints, ATM overexpression at PolyUb-TP53** (curators DOWN;
  curator and experimental): 0.0005× → 50×. **Non-converged in both runs**; the fallback picks
  the value. Not a mechanism of the rule.
- **Tuning S phase, CCNE1 overexpression at the pre-replicative complex and p-FZR1,p-RB1**
  (curators NORMAL): 0.93× (non-converged) → 45× (converged). Cyclin E driving RB1 and FZR1
  phosphorylation UP is the curated direction of those reactions.

## Reading

The rule does what it was built for: every CKS1B drain is gone, the 4 CCNE1 experimental cases
specs/078 fixed are fixed, and CBL's EGFR and MET cases do not move (specs/078 lost 39 there).
It gains on tuning (+22, p 0.0001) and experimental (+3). The held-out shortfall is 2 cases at
p 0.63: one right-for-the-wrong-reason SMAD7 case, one SHH case where the arm reads the biology
and the curators do not, and one threshold pair.

`ulu88` adds specs/088's MDM2 fixes (+6 curator, +2 experimental in Cell Cycle Checkpoints) and
churns EGFR (0, 9 / 9), and loses CCND1 and CDK4 knockout at CCNA1 by a threshold. Its totals
equal `ulu`'s.

Adoption is Adam's decision; the gate as registered fails.

## Decision

Adam, 2026-10-22: "ok Lets go with your recommendation", i.e. adopt `ulu` (the generator rule
alone, not with specs/088) despite the 2-case held-out shortfall. `LNG_UB_LIGASE_UNIT` is on
by default in LNG PR #122. Once it is merged, the canonical catalog is rebuilt from generator
main and checked to reproduce `ulu` on every prediction.
