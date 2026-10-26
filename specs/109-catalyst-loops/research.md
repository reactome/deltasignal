# specs/109 results

Solver commit 5b63781 (rule at 9991c0f), canonical build `20261024-2211_7f4b28a`. Control: `dflt102`
(defaults, same build). Arm: `cat109` (`DS_SELF_FED_MODE=catalyst_loop`).

## Behaviour checks (G1)

| catalog | arm | converged | rails | monotone | relabel (max fold difference) |
|---|---|---|---|---|---|
| benchmark (16 pathways) | rule off | 476 / 476 | 30,737 / 50,392 (61.0%) | 21,201 / 25,196 | 3.0e-12 |
| benchmark | `catalyst_loop` | 476 / 476 | 30,489 (60.5%) | 21,201 | 3.0e-12 |
| sample (15 pathways) | rule off | 443 / 443 | 6,758 / 16,574 (40.8%) | 8,275 / 8,287 | 7.8e-14 |
| sample | `catalyst_loop` | 443 / 443 | 6,315 (38.1%) | 8,287 | 7.8e-14 |

**G1 passes** on both catalogs: nothing is worse, rails fall, relabelling is exact within 1e-9.

Where it acts, the rails go: PAK1 (R-HSA-373760) 32 → 0 of 32, HNF1B ⇄ ONECUT1 (R-HSA-186712) 40 → 0,
R-HSA-194315 64 → 0, CEBPA (R-HSA-9616222) 208 → 0, TIMM22 (R-HSA-1268020) 15 → 0; partly in caspase
activation (R-HSA-140534 27 → 12, R-HSA-418889 180 → 108), ISG15 (216 → 144) and LHX1/OSR1/PAX2 (128 → 92).

Where it does not: in 10 of the 16 benchmark pathways no selected readout moves, including NHEJ (848 of
848 pairs railed), TGF-beta, DSB repair and the deubiquitination pathways. In those components
the scope test holds (checked by solving one perturbation each: R-HSA-5617472, 5693571, 170834, 177929,
75153 report qualifying components), but the self-fed test finds no node. Signal enters these loops
at several points, so no node is fed only by what it feeds, and nothing is held. The loop still closes
with gain 1. Rule A's dominator test finds the closing edge of a loop with one way in (autocatalysis, a
gene switching itself on) and of no loop with several.

## Benchmark (G2–G5)

| axis | baseline | arm | net | fixed / broke | p |
|---|---|---|---|---|---|
| curator held-out | 0.9079 | 0.9079 | 0 | 0 / 0 | 1 |
| curator tuning | 0.8765 | 0.8765 | 0 | 0 / 0 | 1 |
| curator tuning, RAF set aside | 0.8779 | 0.8779 | 0 | 0 / 0 | 1 |
| experimental | – | – | 0 | – | – |

Curator: 44 edge-solves held across 1,725 solves; 25 valid values move (DSB repair 12, caspase activation
via death receptors 7, via dependence receptors 6), none changes class. Experimental: no solve held an
edge, so the bench's guard stopped the axis ("the arm is the control"); the axis is identical to the
control by construction.

- **G2 holds:** 12 scored cases come off a rail, and all 12 keep their direction.
- **G3, G4 hold:** no held-out loss, so nothing to justify; experimental unchanged.
- **G5:** raw scores above.

**The gate passes.**

## Predictions

- C1 holds: in the selection, the scope test covers 15 of the 17 benchmark pathways and 13 of the 13 sample
  pathways that specs/106 lists. It adds R-HSA-5617472 on the benchmark and two sample pathways, and misses
  R-HSA-74752 and R-HSA-8848021.
- C2 fails: rails fall by 0.8% of their level on the benchmark and 6.6% on the sample, not by a quarter.
- C3 holds: curator held-out 0, experimental 0.
- C4 holds: 12 of 12.

## Reading

The rule is safe and does what it was written for in the loops it reaches. It reaches few: the small
autocatalytic and self-activating loops, not the large catalyst-closed components where most of the rails
are. Making it reach those would need a different choice of which edge closes a loop with several
entries (for example the catalyst edges that point backwards in the solver's flow order, specs/013). That
would be a new rule, pre-registered before it is measured. Adoption of this one is Adam's.

## Amendment 1 results (solver b5d4a44, arm `cat109b`)

**Behaviour checks** on the re-made selections (NHEJ and TGF-beta no longer selected; benchmark 10 pathways,
sample 9):

| catalog | arm | converged | rails | monotone | relabel |
|---|---|---|---|---|---|
| benchmark | rule off | 278 / 278 | 24,888 / 40,104 (62.1%) | 16,233 / 20,052 | 3.0e-12 |
| benchmark | `catalyst_loop` | 278 / 278 | 24,665 (61.5%) | 16,233 | 3.0e-12 |
| sample | rule off | 273 / 273 | 1,161 / 1,576 (73.7%) | 776 / 788 | 2.2e-14 |
| sample | `catalyst_loop` | 273 / 273 | 718 (45.6%) | 788 / 788 | 2.2e-14 |

G1 passes. Five benchmark components that are catalyst-closed with the fixed pair test are still untouched
(Signaling by EGFR, Signaling by WNT, anterior HOX activation, DSB repair, DNA repair): signal enters each at
several points, so no node is self-fed.

**Benchmark:** identical to `cat109`, so the defect never fired on a scored solve. Curator held-out 0, tuning 0,
experimental 0 (no experimental solve held an edge; the bench stopped the axis as the control). 25 values move,
none changes class, and all 12 that come off a rail keep their direction. **The gate passes.**

**Predictions:** A1 holds (NHEJ and TGF-beta drop out of the selection). A2 holds (no class change). A3 holds on
the sample (rails fall by 38%) and fails on the benchmark (by 0.9%).

## What the traces say about the large loops

NHEJ and TGF-beta are not catalyst loops. In each, a complex is the input and the catalyst of its own
reaction: the receptor complex phosphorylating its own SMAD2/3 (R-HSA-170868), MTMR4 dephosphorylating the
p-SMAD2/3 it holds (R-HSA-2187401), and the DNA-PK complexes in NHEJ (R-HSA-5693575, 5687183, 5687360,
5693604). Read with the fixed pair test, they are material cycles: SMAD2/3 is a phosphorylation cycle, and Ku is
released by ligation and rebinds. With specs/106 amendment 2, material-closed positive loops are the largest
unhandled class (benchmark 28 loops in 14 pathways, sample 46 in 20). The large loops that really are closed
only through catalysts are EGFR, WNT, HOX activation and DNA/DSB repair.

**Decision (Adam, 2026-10-26, at the freeze, specs/111): not adopted.** It passes the gate but changes no class on any axis.
