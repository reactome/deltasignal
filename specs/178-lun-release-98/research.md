# specs/178 research: Lun 2019 on Release 98

Run 2026-11-27. Build `20261126-2322_d6d026b_top98` (Release 98, checksum 525378806530, 28 of 28 pathways, generator
`d6d026b`); predictions `results/906698aa/lun178_r98` (967 of 991 solves converged). Release 97 beside it is
specs/175's file (`20261027-1346_8b6527e_top114`, `results/356cb2b2/lun175_abund`). Scoring output is outside the
repo, under `analysis/178/`.

**Half B is FIT** (seen once in specs/177). Nothing below is a held-out test.

## Verdict

Every specs/177 verdict holds on Release 98. The Release 98 figures are the Lun figures quoted from here on.

## How the releases differ for these pairs

| | Release 97 | Release 98 |
|---|---|---|
| prediction rows | 11,028 | 11,051 |
| protein–readout pairs | 6,953 | 6,955 |
| proteins | 428 | 428 |

- The 23 new rows are all in R-HSA-9909396: p-AMPKα (16) and p-CREB/ATF1 (7). They add two scored pairs, both on
  half B (CSNK1D and SNF1LK at p-AMPKα, both called no change).
- The model's call moves on 7 pairs:
  - six go from no change to DOWN at p-MAPKAPK2 (ACVR1B, ACVR1C, ACVR2A, ACVR2B, TGFBR1, TGFBR2, in Signal
    Transduction). Five of the six measure no change and one measures UP, so none is right;
  - SRC at p-MKK3/6 goes from DOWN to no change, and it measures DOWN.
  - Fixed 0, broken 6 (the seventh was wrong before and after).
- Reachability gains 6 pairs, all ZAP70 (p-JNK, p-MAPKAPK2, p-MKK3/6, p-MKK4, p-PLCγ2, p-p38); 3 of the 6 measure UP.
- The tyrosine-kinase class table is identical between releases (60 of 636 proteins; only the header comment
  differs). The readout table does not depend on the release.

## The specs/177 gates, Release 98 beside Release 97

Three-class balanced accuracy; intervals are protein bootstraps.

| | Half A, R97 | Half A, R98 | Half B, R97 | Half B, R98 |
|---|---|---|---|---|
| pairs | 3,127 | 3,127 | 3,633 | 3,635 |
| candidate (model + rule) | 0.435 | 0.433 | 0.420 | 0.420 |
| model alone | 0.384 | 0.383 | 0.358 | 0.357 |
| reachability | 0.380 | 0.381 | 0.354 | 0.354 |
| rule alone | 0.431 | 0.431 | 0.426 | 0.426 |
| **candidate − reachability** | +0.056 [+0.030, +0.083] | +0.052 [+0.027, +0.079] | +0.066 [+0.036, +0.096] | +0.066 [+0.036, +0.095] |
| candidate − model | +0.051 [+0.031, +0.072] | +0.051 [+0.031, +0.072] | +0.063 [+0.039, +0.085] | +0.062 [+0.039, +0.085] |
| model − reachability | +0.005 [−0.011, +0.021] | +0.001 [−0.014, +0.018] | +0.004 [−0.013, +0.020] | +0.004 [−0.013, +0.020] |
| added calls, right of all | 89 of 120 | 89 of 120 | 89 of 103 | 89 of 103 |
| route enrichment | 1.55 [1.22, 1.96] | 1.56 [1.23, 1.97] | 1.50 [1.17, 1.88] | 1.50 [1.16, 1.88] |

- All three gates pass on both halves on both releases.
- The model alone still does not separate from reachability on either half.
- Tyrosine kinase × tyrosine site reads UP on 72% of pairs on both halves; the class table is unchanged, because it
  is built from measurements and the class table only.
- Stricter measured calls (Release 98): candidate − reachability is +0.078 at 2δ and +0.115 at 3δ on half A, +0.127
  and +0.152 on half B, every interval above 0.

## specs/179's engagement test on Release 98

Registered in specs/179 to be repeated here. Model change calls, direct pairs left out.

| | Release 97 | Release 98 |
|---|---|---|
| engaged / not engaged proteins | 36 / 49 | 36 / 49 |
| precision, engaged | 36.2% | 35.9% |
| precision, not engaged | 6.0% | 5.8% |
| difference, pooled | +0.30 [+0.21, +0.40] | +0.30 [+0.20, +0.40] |
| difference, half A / half B | +0.29 / +0.39 | +0.29 / +0.39 |
| withholding gate, fixed / broken | 175 / 12 | 180 / 12 |
| recall, half B, before → after | 0.176 → 0.124 | 0.176 → 0.124 |

The same reading: engagement separates right from wrong change calls on both halves, and withholding calls for
not-engaged proteins still fails the recall guard on half B.

## What this does not change

The three decisions specs/177 and specs/179 left for Adam are unchanged: the Lun protocol of record, serving the
tyrosine-kinase statement on over-expression queries, and reporting Lun split by engagement.
