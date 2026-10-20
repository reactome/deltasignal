# specs/081 research: arm rcp

Solver `0174e2c` (main `fdd759e` plus this branch's multi-base pools, recycle edges and
bench root rule), variant build `20261020-1517_55b0f1b_rcp`, run 2026-10-20 with
`scripts/run_arm.sh --catalog`, against `results/3e96e7c/dflt` on canonical. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`. Results are in
`results/0174e2c/rcp`. Both axes log 68 recycle edge-solves.

## Results

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **0** (0 / 0) | 1 |
| curator tuning | **+70** (70 / 0) | < 1e-4 |
| experimental (RAF set aside; RAF does not move) | **+8** (10 / 2) | 0.039 |

| | before | after |
|---|---|---|
| curator tuning, in release (accuracy, macro-F1) | 79.53%, 0.7822 | **80.92%, 0.7969** |
| curator all, in release | 87.97%, 0.8531 | 88.26%, 0.8565 |
| experimental, in release | 71.04%, 0.6378 | **71.99%, 0.6428** |

Only Mitotic G1 moves, on both axes.

**The registered gate passes:** held-out net 0 ≥ 0, experimental +8 ≥ −5, and no pathway
below −10 (Mitotic G1 is +70 / +8).

## Predictions

1. At least 60 of the kinase-up / phosphatase-down cases fixed: **holds**. CCND1 and CDK4
   overexpression (17 each), PPP2R1A knockout, CDKN2A knockout, and others make 70 in all.
2. Kinase-down / phosphatase-up cases not fixed: **partly wrong.**
   - PPP2R1A overexpression gains 8.
   - CCND1 and CDK4 knockout gain a few cases each at CCNA1/2, from 0.865 to 0.841–0.849.
   - The φ₀ asymmetry still limits them; most stay NORMAL.
3. Mitotic G1 curator net at least +40: **holds**, +70. MP-BioPath's networks score 146 of
   the 187; we now get 70 of them on curator.
4. Held-out 0, nothing outside Mitotic G1 moves: **holds**.

## Experimental cases (traced)

All 12 changed cases are in Mitotic G1:

| | cases | change |
|---|---|---|
| fixed | CCND1 overexpression at CCNE1, DHFR, MYBL2, CDT1, CDK1, TK1 (6) | NORMAL to UP |
| fixed | CDK4 overexpression at CDT1 (1) | NORMAL to UP |
| fixed | CCND1 and CDK4 knockout at CCNA1 (2) | 0.865 to 0.849 / 0.841 (DOWN) |
| fixed | CDKN1B overexpression at CDK1 (1) | to UP; see below |
| broken | CCND1 overexpression at CDC25A and CDC6 (2) | to UP |

In the two broken cases the experiment says NORMAL and the curators say UP; we now give
the curators' answer, as specs/079 found for most experimental errors.

**CDKN1B (p27) reads opposite to the curators**, both before and after. Under `rcp`, p27
overexpression raises the E2F targets and knockout lowers them slightly. The curators
expect the reverse: an inhibitor. The cause is curated. The active kinase complex
(R-HSA-8949896, "p-T172-CDK4,p-T177-CDK6:CCND:CDKN1A,CDKN1B,(CDKN1C)") is a set that
includes forms containing p21/p27, the assembly-factor role of p21/p27 for cyclin D:CDK4/6.
So more p27 is more of that kinase complex.
- p27's curator cases were wrong before (NORMAL) and are wrong now (inverted), so net 0.
- The one experimental "fix" (CDK1) is reached through that curated assembly-factor
  reading.
- Not a defect of this rule. It is how the kinase is curated, and it is a question for
  the curators, not a rule here.

## Reading

Closing the p130/p107 cycle as a conserved pool fixes 70 tuning-pathway curator cases with
none broken, and experimental +8 (p 0.039). The held-out split cannot move: Mitotic G1 is one
of the paper's ten. The gain is concentrated in one pathway and two cycles, the only ones
the rule finds on canonical (`count_recycle.py`). Adoption is Adam's. Adopting means the
generator default `LNG_RECYCLE_POOLS=1` and a rebuilt canonical; the solver and bench
changes are inert without recycle pools.

## Adoption (2026-10-20)

Adam, 2026-10-20: "go for it and continue". The adoption steps:
1. merge DeltaSignal (solver and bench changes, inert without recycle pools);
2. make `LNG_RECYCLE_POOLS=1` the generator default;
3. rebuild canonical and check it reproduces `rcp` prediction for prediction.

Reported with three caveats:
- the curator gain is on a tuning pathway, and the rule was traced from it;
- about half of the 187 cases still fail on the φ₀ asymmetry (a separate question);
- CDKN1B reads inverted through the curated p21/p27 assembly-factor role.
