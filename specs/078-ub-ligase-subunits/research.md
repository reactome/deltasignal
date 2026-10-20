# specs/078 research: arms ubx (gated) and ubc

Solver `144c416` (main `3e96e7c` defaults plus specs only), run 2026-10-20 with
`scripts/run_arm.sh --catalog`, against `results/3e96e7c/dflt` on canonical
`20261017-1246_5e663d9`. Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade
--pair-on-valid`: 589 curator and 4 experimental cases are invalid in both and unpaired.

| arm | build | results |
|---|---|---|
| ubx | `20261020-0756_9770d28_ubx` | `results/144c416/ubx` |
| ubc | `20261020-0717_9770d28_ubc` | `results/144c416/ubc` |

## Results

| split | ubx net (fixed / broken) | p | ubc net | p |
|---|---|---|---|---|
| curator held-out | **−44** (12 / 56) | < 1e-4 | −44 (12 / 56) | < 1e-4 |
| curator tuning | +19 (23 / 4) | 0.0003 | +23 (27 / 4) | < 1e-4 |
| experimental (RAF set aside; RAF does not move) | **+5** (6 / 1) | 0.125 | +6 (7 / 1) | 0.070 |

By pathway (curator, `ubx`):
- **Losses:** EGFR −29 (CBL, EPS15 and others), RUNX2 −11 (ESR1 knockout 10), MET −10 (CBL), Hedgehog −2.
- **Gains:** TP53 +10, Mitotic G1 +7, NOTCH1 +7 (MIB1 knockout), S phase +1, Fanconi +2, WNT +2.

Experimental: Mitotic G1 +4 (CCNE1 overexpression 2, CDK4 knockout 2), S phase +1
(CCNE1 knockout +2, CCNE1 overexpression −1). `ubc` adds a CCND1 knockout.

**The registered gate fails** (amendment 1: gated on `ubx`):
- held-out −44 < 0;
- EGFR −29 and RUNX2 −11 are below −10.

`ubc` fails the same way. **Not adopted.**

## Predictions

1. At least 2 of the 4 traced CCNE1 experimental cases fixed: **holds**, all 4 on both arms.
2. CBL, EPS15 and SH3GL1 lose again, MET/EGFR −20 or worse: **holds**, −39 together.
3. No prediction for held-out: −44.

## The losses

- **EGFR and MET (−39)** are specs/057's mechanism, unchanged by specs/077. Curators expect a CBL-complex
  knockout to raise receptor signalling. On canonical it does, but only through the
  removed drains. Without them CBL depletes only the ubiquitinated receptor, which has
  no free node, so CBL has no effect. Getting these right needs a degradation rule
  for the CBL class, which neither spec registered.
- **RUNX2 ESR1 knockout (−10)** is a threshold shift: the readouts go from 0.81× (DOWN) to 0.856× (NORMAL),
  and 0.905× to 0.928× for the rest. The removed RUNX2 edges are the class this spec targets:
  RUNX2:STUB1, RUNX2:SMURF1 and RUNX2:HIVEP3:WWP1 draining free STUB1, SMURF1, WWP1 and
  HIVEP3, the ligases and adaptor. The drain on RUNX2 itself is kept. Not traced further.

## Reading

Removing ubiquitin drains on a ligase's own subunits fixes what specs/077 exposed
(every CCNE1 case, experimental +5). It also gains in TP53, Mitotic G1 and NOTCH1. But it
still costs the CBL-class negative regulators, which were right through the same wrong
edges, and the held-out net is −44. As in specs/057, the generator fix is more faithful
and loses on the curator set until the CBL class has a rule of its own.
