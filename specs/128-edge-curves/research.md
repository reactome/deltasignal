# specs/128 research: response curves on edges

Frozen canonical build `20261026-1349_8b6527e`, solver `110bae4` (every arm), curator, experimental and korkut axes.
Arms `ec_off`, `ec_id`, `ec_inh`, `ec_sub90`, `ec_all90`, `ec_all75`; tables in `bench/analysis/128/`. Scored with
`bench/analysis/score_128.py` (tested in `bench/test_score_128.py`), which uses specs/124's comparison for curator and
experimental; rail-dependence uses the specs/124 mild-dose run (`ce631f0/a124_mild`). Cases: curator held-out 18,481,
tuning 5,030, experimental 845, korkut 11,062 (always-NOT right on 10,049).

## D0

`ec_off` reproduces the frozen run `05dc867` case for case (24,100 curator, 849 experimental: validity, call, value as
printed, convergence). `ec_id`, every class listed at its identity setting, reproduces `ec_off` exactly on all three
axes, korkut included. **P1 holds.**

## Results against `ec_off`

| arm | curator held-out right (fixed / broke, p) | tuning | experimental | korkut right (fixed / broke, p); change calls |
|---|---|---|---|---|
| `ec_off` | 16,778 / 18,481 | 4,415 / 5,030 | 643 / 845 | 7,117 / 11,062; 3,922 |
| `ec_inh` (inhibition smoothed) | 16,777 (0 / 1, p 1) | 4,413 (10 / 12) | 644 (2 / 1) | 7,124 (22 / 15, p 0.32); 3,905 |
| `ec_sub90` (substrates, slope 0.9) | 16,685 (34 / 127, p 8e-14) | 4,331 (23 / 107) | 641 (6 / 8, p 0.79) | 7,139 (48 / 26, p 0.014); 3,890 |
| `ec_all90` (all classes, 0.9) | 16,653 (46 / 171, p 4e-18) | 4,197 (48 / 266) | 629 (8 / 22, p 0.016) | 7,157 (67 / 27, p 5e-5); 3,870 |
| `ec_all75` (all classes, 0.75) | 16,179 (102 / 701, p 1e-110) | 3,681 (81 / 815) | 529 (26 / 140, p 4e-20) | 7,470 (371 / 18, p 7e-87); 3,527 |

Setting RAF aside changes the experimental losses only (`ec_all75` 26 / 124) and leaves every held-out figure as is.
Pathways with a net change: 3 (`ec_inh`), 25, 32, 61. Under every curve the five held-out solves that did not converge
now do. Readouts at a rail (held-out curator): 3,009 → 2,985 / 1,906 / 1,707 / 1,647; experimental 216 → 203 / 160 /
131 / 120.

| prediction | result | verdict |
|---|---|---|
| **P1.** `ec_id` reproduces `ec_off` exactly | identical on all three axes | **holds** |
| **P2.** `ec_inh` changes < 1% of held-out curator and experimental calls | 0.01% and 0.47% | **holds**: today's `divide` with its 10× cap is, in effect, the k = 9, m = 1 curve |
| **P3.** `ec_all75` loses held-out curator calls at p < 0.05 | −599 (102 / 701), p 1e-110; the milder arms lose too (−93, −125) | **holds** |
| **P4.** rails fall by more than half under `ec_all75` | −45% curator held-out, −44% experimental | **fails** (narrowly) |
| **P5.** `ec_all75` makes fewer korkut change calls and is not less accurate | 3,922 → 3,527 change calls; right 7,117 → 7,470 (+353, p 7e-87) | **holds** |

## What the losses are

Almost every lost case is a correct change call weakened into "no change", not a direction flip. Held-out curator,
`ec_sub90`: 71 UP → NORMAL and 52 DOWN → NORMAL of 127; `ec_all75`: 535 UP → NORMAL, 150 DOWN → NORMAL, 13 UP → DOWN
of 701. Experimental `ec_all75`: 105 UP → NORMAL, 34 DOWN → NORMAL, 1 UP → DOWN of 140.

The arithmetic is the curve compounding along a route: with K = 9 a hop maps 80× to 9×, 9× to 5×, 5× to 3.6×, ..., so a
change decays towards rest with every curved hop. Examples: CCND1 over-expression in Pre-NOTCH reads 13.65× today, 1.10×
under `ec_sub90`; BRCA1 over-expression in TP53 (experimental) 12.8× → 1.12×. I tried to count the curved hops on the
CCND1 route from `logic_network.csv` and found none: the route runs through edges the solver adds (pools, twins), which
that file does not hold. So the hop count is not traced here; the class of the change is what is shown.

Under `ec_all75`, 130 of the 844 changed held-out curator cases were rail-dependent (specs/108) in the control: some
of what contraction removes was decided by a rail, as intended; most was not.

## Reading

**Nothing is adopted**; the default stays `off`. The machinery is exact (P1), so stage 3 can start.

What this tells stage 3:

1. **The two kinds of ground truth pull in opposite directions.** Curator expectations reward changes that pass at
   full strength however far they travel: every contraction loses curator calls, in proportion to its strength. The
   measured drug data (korkut, FIT) reward attenuation: every contraction gains there, and the strongest gains most,
   by making fewer change calls far from the target (specs/115's finding). Experimental (MP-BioPath's measurements) sits
   between: unchanged under substrate curves, losing under the stronger ones.
2. **Even the best arm is far from always-NOT on korkut** (7,470 against 10,049 of 11,062): attenuation by class alone
   does not make the model calibrated; it only makes it call fewer changes.
3. **The gate for learned parameters needs deciding before stage 3.** If fitted curves are judged by curator agreement
   (SELECT), they will be pushed back to the identity, whatever the measurements say. specs/125 already sets the judge as
   cross-validation on FIT data and the locked tests; this result says that choice decides the outcome, so it should be
   confirmed (Adam's decision), and curator agreement reported, not used as a no-harm gate.
4. **Inhibition is already where a curve would put it** (P2): the open parameters are on the activating side.
