# specs/090 research: arm psa

Canonical `20261022-2011_233665c`, solver `dccf17d` with `--server DS_POOL_SUPPLY_AND=on`,
against `results/8a94f2e/dflt` (main defaults). Results `results/dccf17d/psa`. Read with
`holdout_report.py --set-aside RAF_MAP_kinase_cascade --pair-on-valid`; cases with
`specs/089-ub-ligase-unit/diff_cases.py`. The solves formed 204 AND groups (state-solves).

| split | net (fixed / broken) | p |
|---|---|---|
| curator held-out | 0 (0 / 0) | 1 |
| curator tuning | **−6 (0 / 6)** | 0.031 |
| experimental (RAF set aside) | 0 (0 / 0) | 1 |

All 6 are in Mitotic G1:
- RBL2 overexpression at CCNA1 (68891), its CDK2 complex (68889) and 389106: DOWN → UP
  (0.36× → 27.6×, 0.005× → 1.25×, 0.52× → 3.34×);
- PPP2R1A overexpression (2) and CDKN2A knockout (1): threshold shifts (0.84× → 0.86×,
  1.15× → 1.13×).

**The gate as registered passes** (held-out 0, experimental 0, worst pathway −6), but the net
is −6 at p 0.03, and the registered prediction about overexpression is refuted.

## Predictions

1. Only Mitotic G1 moves: **holds**.
2. RBL1/RBL2 knockout: the pools read 0; the curator cases stay right; the 20 experimental cases
   stay wrong: **holds** (no knockout case changes class).
3. RBL1/RBL2 overexpression reads lower at the E2F targets: **refuted** for RBL2 (below).
4. Net within ±5 on every axis: **refuted**, tuning −6.

## RBL2 overexpression (traced, `trace_rbl2_oe.jl`)

Pinning the RBL2 root at 80×, the pool that moves is p107's, not p130's:

| | p107 pool forms | CCNA1 (68891 copies) |
|---|---|---|
| off | 0.68× | up to 2.2× |
| on | **0.013×** | up to 100× |

The p107 pool's RBL1 copy is not a starting copy (review finding 18). It is the product of
"RBL1 gene expression is repressed by the DREAM complex" (R-HSA-8964525), and DREAM contains
p130. So p130 overexpression represses RBL1 transcription, which is curated biology. Under `off`
that repression was diluted by the partners in the mean; under `on` it empties the p107 pool,
and the E2F targets that p107 held are de-repressed. The p130 pool's own repression of the same
targets does not offset it here, which is not traced further.

So the rule is more faithful to the recycle pool, and it exposes a curated cross-repression
(p130 → DREAM ⊣ RBL1) that the model then reads more strongly than the curators' answer
(RBL2 overexpression: E2F targets DOWN).

## Reading

Not recommended for adoption: −6 on tuning, nothing gained on any axis, and the
overexpression prediction is refuted. The knockout behaviour it was written for (an RBL1 or RBL2
knockout empties its pool) does not change any class, because those cases were already right
on the curator axis and the experimental axis disagrees with the curators there.

## Decision

Adam, 2026-10-23: "lets go with 1", i.e. merge as a not-adopted record. `DS_POOL_SUPPLY_AND`
stays an opt-in, default off.
