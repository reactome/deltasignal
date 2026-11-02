# specs/124 research: how much rests on the 100× ceiling

Build `20261026-1349_8b6527e` (frozen canonical), both axes, every other setting the default.
- **Ceiling arms:** solver `ce631f0`. `a124_c100`, `a124_c10`, `a124_c30` and `a124_c1000`, plus the mild-dose run
  `a124_mild` (knockdown 0.5×, over-expression 2×, default ceiling) for specs/108 rail-dependence.
- **Dose-only arms** (amendment 1): `f81185a`, which differs from `ce631f0` only in the scorer. `a124_d24` and `a124_d8`.
- **Scoring:** `bench/analysis/score_124.py`, tested in `bench/test_score_124.py`. Cases are those valid in both
  runs: curator held-out 18,481, tuning 5,030, experimental 845.

## D0: the control reproduces the frozen model exactly

`c100` against the frozen run `05dc867`: all 24,100 curator and 849 experimental cases are identical in validity,
call, value as printed, and convergence. The resting-level knob is exactly neutral at its default.

## Results against the registered predictions

| arm (ceiling, over-expression pin) | curator held-out: right (changed, fixed / broke, p) | tuning | experimental |
|---|---|---|---|
| c100 (100×, 80×) | 16,778 / 18,481 (0.9079) | 4,415 / 5,030 (0.8777) | 643 / 845 (0.7609) |
| c1000 (1,000×, 80×) | 16,776 (2, 0 / 2, p 0.5) | 4,406 (11, 1 / 10, p 0.012) | 643 (1, 0 / 0, p 1) |
| c30 (30×, 24×) | 16,762 (21, 1 / 17, p 0.0001) | 4,405 (11, 0 / 10, p 0.002) | 641 (2, 0 / 2, p 0.5) |
| c10 (10×, 8×) | 16,739 (102, 29 / 68, p 0.0001) | 4,325 (103, 5 / 95, p 1e-22) | 638 (5, 0 / 5, p 0.06) |

Setting RAF aside removes one changed case (a c10 tuning case) and changes nothing else. Pathways with a net change:
5 under c1000 (largest Mitotic G1 −8), 8 under c30 (TP53 −10), 14 under c10 (TP53 −84, over-expression). Under c1000,
18 held-out solves that converge at 100× no longer converge (18,476 → 18,458); none of them changes class.

| prediction | result | verdict |
|---|---|---|
| **P1.** c30 and c1000 change < 2% of scored held-out curator and experimental cases | c30: 0.11% and 0.24%; c1000: 0.01% and 0.12% | **holds** |
| **P2.** c10 changes 2–8% on each axis, mostly over-expression | 0.55% held-out curator, 0.59% experimental (tuning 2.05%); over-expression 94 of 102 and 5 of 5 | **fails** (fewer than predicted); the over-expression part holds |
| **P3.** At c30 and c1000, more than half the changed cases are rail-dependent (specs/108) | 0 of 21 + 2 (c30) and 0 of 2 + 1 (c1000) | **fails**, and the measure misses the cases that matter (below) |
| **P4.** The share of readouts at a rail moves < 5 points between c30 and c1000 | held-out 16.1% vs 16.0%; experimental 25.3% vs 24.3%; top rail alone 6.7% vs 6.2% (curator) | **holds** |
| **P5.** Under c1000, > 10% of scored readouts move by > 10% in fold | curator 15.2%, experimental 40.5% (Spearman of log fold 0.999 and 0.991) | **holds** |

The rail-dependence measure itself checks out on this build. 7.6% of correct curator cases and 9.5% of correct
experimental cases are rail-dependent, as specs/108 found (8.0%, 9.5%).

## Amendment 1: dose or ceiling?

c30 and c10 lower the over-expression dose as well as the ceiling. Each changed case is attributed by the dose-only
arm at the same dose and the default ceiling:

| | changes vs c100 | from the dose (dose arm changes the same way) | from the ceiling | ceiling-only contrast (arm vs dose arm): held-out / tuning / experimental |
|---|---|---|---|---|
| c30 (d24) | 32 curator, 2 experimental | 30 curator, 2 experimental | 2 curator | 2 / 1 / 0 |
| c10 (d8) | 205 curator, 5 experimental | 177 curator, 5 experimental | 28 curator | 39 (10 fixed / 26 broken, p 0.011) / 2 / 1 |

The dose alone (at the default ceiling) changes more than the ceiling does. d24: held-out 19, tuning 12, experimental
2. d8: held-out 87, tuning 103 (TP53), experimental 6 (0 fixed / 6 broken). The 80× over-expression pin is as
arbitrary as the ceiling, and it matters more.

## Traced: what the ceiling decides

Every change that comes from the ceiling falls into one of two patterns.

1. **Threshold drift.** A readout near a cut-off moves a little and crosses it. Under c1000: Mitotic G1 ABL1, JAK2,
   PTK6 and SRC knockouts move 0.841 → 0.859 (8 cases); TP53 AKT1/AKT2 knockouts move 1.155 → 1.150; RUNX2 PPM1D
   moves 0.856 → 0.850. These are small shifts from rails elsewhere in the network.
2. **A product or ratio of bounds.** A readout combines a node at the rail with one at another bound, so its value is
   the ratio of two arbitrary constants:
   - **WNT, WNT5A over-expression** (`R-HSA-4411380` p-TCF/LEF:CTNNB1, made by "NLK phosphorylates TCF/LEF").
     Substrate TCF/LEF:CTNNB1 at 0.00258×, catalyst pT298-NLK dimer at the rail. The product reads 0.258 (DOWN) at
     100× and 2.580 (UP, the curators' answer) at 1,000×.
   - **RUNX1, RUNX1 over-expression** (IFNG). Activator RUNX1:CBFB:IFNG gene at 34.8×, inhibitor
     RUNX1:CBFB:FOXP3:IFNG gene at the rail. Divided, it reads 8.49 (UP, right) at 100× and 0.849 (DOWN) at 1,000×.
   - **EGFR, EGFR over-expression at 8×** (`R-HSA-180523`). Substrate at the depletion floor (0.1×), catalyst at the
     rail. It reads 0.1 × 100 = 10 (UP) at 100× and 0.1 × 10 = 1.0 (NORMAL) at 10×. Same for ARHGEF7, CBL, EPS15 and
     CALR under c10.

   In these cases neither answer is credible: the call is the ratio of the depletion floor, or a near-zero substrate,
   to the ceiling. specs/108's test cannot see them, because it asks whether the readout itself rails at a mild
   dose, and here the rail is one step upstream.

Bound arithmetic shows in the control's values too: 1,454 curator readouts read exactly 100.000, 409 read 80.000 (the
pinned dose), 137 read 50.500 ((1 + 100) / 2, a railed copy averaged with one at rest), and 128 read 10.000.

## Reading (fixed in advance, applied)

- **P1 holds: the up/down results do not rest on the 100.** At 1,000× two held-out cases change, and at 30× (once the
  dose is separated) two. The ceiling can enter specs/125 as one more parameter (a resting level per node class)
  without first re-examining every headline.
- **P5 holds: magnitudes do depend on it.** Under 1,000×, 15% of curator and 40% of experimental readouts move by more
  than 10% in fold. Every correlation-type score (CPPA, decryptM magnitudes) is reported with the ceiling named.
- **Not registered, recorded for specs/125:**
  - The over-expression dose (80×) changes more calls than the ceiling, so it should be measured or learned too.
  - A class of predictions is decided by ratios of the model's bounds (ceiling, depletion floor, near-zero
    substrates) one step upstream of the readout. A rail-dependence test that only looks at the readout misses them.
  - These are the cases a learned, contractive model (loop gain below 1; no node pinned at a bound) would change.

Nothing is adopted. `DS_BASELINE` stays 0.01.
