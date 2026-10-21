# specs/086 research: arms dflt086 and rlp

Solver `e1ba794` (`CYCLE_MAX_STEPS` 10). Run 2026-10-21 with `scripts/run_arm.sh`.

## Build check and baseline

- **Variant build** `20261021-1203_aa3f605_rlp` (LNG `aa3f605`, `--env LNG_RELEASE_POOLS=1`),
  compared with canonical by `compare_build.py`:
  - Signaling by WNT: 1 release pool, 8 state and intermediate rows, 11 transitions, 1
    `released_twin` and 1 `recycle` edge;
  - L1CAM interactions: 1 release pool, 3 rows, 3 transitions, the same two edges;
  - nothing removed, no other pathway differs.

  Only 2 of the 9 counted entities qualify; the other 7 fail the generator's
  modification-signature test. Because the arms ran back to back, this check was run
  before the arm's results were read rather than before the arm.
- **Baseline** `dflt086` on canonical equals `results/1b52a97/cct2` on every prediction
  (24,100 curator, 849 experimental). The path limit change is inert, as registered.

## Arm rlp

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 | 1 |
| curator tuning | **−5** (2 / 7), all WNT (WNT5A, WNT1, WIF1) | 0.18 |
| experimental (RAF set aside) | +2 (2 / 0), WNT (WIF1, WNT1) | 0.5 |

The gate passes on its face (held-out 0, experimental +2, WNT −5).

## Predictions

1. APC, at least 5 of 7 experimental cases fixed: **fails, 0**.
2. WNT curator net at least +15: **fails, −5**.
3. Nothing outside the counted pathways moves: **holds**.

So the rule did not do what it was registered for.

## Why (traced)

The APC values moved towards baseline without crossing it:

| case | before | after |
|---|---|---|
| APC knockout, CTNNB1 [cytosol] | 0.023× | 0.26× |
| APC overexpression, CTNNB1 [cytosol] | 98× | 53× |

The released copy of β-catenin is made by **20 release reaction copies**, one per variant of
the destruction complex (the CK1/GSK3/AXIN member choices). Each consumes its own complex
copy. The generator's route search took the single shortest route, so:
- the pool holds one variant chain;
- the carrier lists 1 of the 20 release copies;
- the released copy still reads the other 19 copies' traffic, inverted as before.

The effect is diluted, not reversed. This is an implementation gap: the spec's M is "the
form the release consumes", and there are 20 such forms.

**Candidate amendment 1 (post hoc, not run):** one release pool per (S, P) holding every variant chain.
- States: each release copy's consumed complex is a base state, as specs/081 does for p130's
  four resting copies.
- Paths: S → … → each M along its own variant route, and each M → S through its release.
- Carrier: P, with all release copies as its release steps.

## Amendment 1: arm rlp2

- **Variant build** `20261021-1334_b751de4_rlp2` (LNG `b751de4`). `compare_build.py`: WNT
  release pool now 141 state and intermediate rows and 220 transitions (all 20 variant
  chains); L1CAM 7 and 9; the same two edges each; nothing removed, no other pathway differs.
- Solver `b30acdd` (`e1ba794` plus spec commits), against `results/e1ba794/dflt086`, results
  `results/b30acdd/rlp2`.

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | **0** | 1 |
| curator tuning | **+73** (84 / 11), all Signaling by WNT | < 1e-4 |
| experimental (RAF set aside) | **+8** (10 / 2), all WNT | 0.039 |

By perturbation:
- **Curator:** APC +28, AMER1 +28, TCF7L2 +24, against losses at WNT1 (−3), WNT5A (−3) and WIF1 (−1).
- **Experimental fixed:** APC 5, AMER1 knockout 2, TCF7L2 overexpression 2, WIF1 overexpression 1.
- **Experimental broken:** APC overexpression at the ubiquitinated phospho-β-catenin:SCF
  complex (UP, right before, now DOWN), and AMER1 overexpression 1.

APC experimental cases:
- **Fixed:** knockout at nuclear β-catenin, cytosolic β-catenin and the canonical targets
  (DOWN → UP); overexpression at nuclear β-catenin and the targets (UP → DOWN).
- **Still wrong:** overexpression at cytosolic β-catenin (now NORMAL, 0.91×; expected DOWN);
  knockout at p-TCF/LEF:CTNNB1 (DOWN → UP, expected NORMAL).
- **Broken:** overexpression at the degradation complex. The pool holds β-catenin bound or
  free; the ubiquitination branch reads the bound phospho form, whose share the pool now sets.

**The gate passes** (held-out 0 ≥ 0, experimental +8, no pathway below −10). Post hoc.

Predictions:
1. APC at least 5 of 7: **holds**, 5.
2. WNT at least +10: **holds**, +73.
3. Nothing outside WNT and L1CAM moves: **holds**. L1CAM moves no scored case.

## Reading

β-catenin treated as a conserved free/bound pool, with the released copy reading the free
state, gets the destruction complex's sign right for APC, AMER1 and TCF7L2: +73 curator
and +8 experimental, nothing broken outside WNT.

Caveats:
- WNT is a tuning pathway, so held-out cannot move.
- The amendment was written after the first arm.
- The rule fires for 2 entities on canonical, and only β-catenin moves a scored case.

Adoption is Adam's: generator default `LNG_RELEASE_POOLS=1`, a canonical rebuild, and
DeltaSignal PR #154 (`CYCLE_MAX_STEPS` 10, inert without release pools).
