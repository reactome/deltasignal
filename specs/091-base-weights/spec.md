# specs/091: a pool's resting copies share one base weight

**Status:** pre-registered 2026-10-23, before the implementation and the arm.

Adam, 2026-10-23: "ok" to candidate I2 of the 2026-10-22 review
(docs/reviews/2026-10-22-review.md, finding 3).

## The defect (review finding 3)

`index_pools` (src/core/cycles.jl) gives every state the baseline weight ρ^d, with ρ = φ/(1−φ)
and d its distance from the nearest base, and every base weight 1 (specs/081: "each resting copy
starts at full weight"). With k bases, a form one step away therefore rests at
ρ/(k + ρ) = 1/(9k + 1), not φ₀ = 0.1, and the largest fold it can reach is 9k + 1. k is the
number of variant copies the generator emits for the resting complex, so a magnitude depends on
isoform multiplicity:

| pool | bases | at-rest share of the next form | largest fold |
|---|---|---|---|
| ordinary pools (122) | 1 | 0.1 | 10× |
| p107 (Mitotic G1) | 2 | 0.053 | 19× |
| L1CAM (R-HSA-373760) | 3 | 0.036 | 28× |
| p130 (Mitotic G1) | 4 | 0.027 | 37× |
| β-catenin (WNT) | 20 | 0.0055 | 181× (above the 100× clamp) |

## The rule (`DS_POOL_BASE_WEIGHT`, default `full`; `split`)

Under `split`, each of a pool's k bases weighs 1/k, so together the resting copies hold what a
single base holds, and the form one step away rests at φ₀ whatever k is. Non-base forms keep
ρ^d. A pool with one base (all 122 ordinary pools) is unchanged bit for bit. The solve reports
`pool_base_weight` and the number of multi-base pools it split.

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `pbw` | canonical `20261022-2011_233665c` | `--server DS_POOL_BASE_WEIGHT=split` | `results/8a94f2e/dflt` |

Read with `holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`, case by case.

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

## Predictions

1. Only Mitotic G1, Signaling by WNT and R-HSA-373760 (L1CAM interactions) move.
2. WNT: APC and AMER1 knockout no longer rail at the 100× clamp at β-catenin targets
   (≤ 10×), and their class does not change (UP either way).
3. Mitotic G1: perturbations of the p107/p130 kinases and phosphatase (CCND1, CDK4, CDKN2A,
   PPP2R1A) move further, because the modified share at rest is larger. No direction is
   predicted for their net.
4. No prediction for any net. This is proposed on faithfulness: a magnitude should not depend
   on how many variant copies the generator emits.
