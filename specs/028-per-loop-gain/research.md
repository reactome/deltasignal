# 028 — A minimal pull towards 1 per loop, not per edge

## Why

Adam, 2026-09-25: add *"a little weight to pull values towards one, but it
would need to be extremely minimal for loops"*. A cycle without inhibitors
should still amplify: up above one, and down below one. For loops with
negative interactions, which act as dampeners, the pull matters less.

`DS_LOOP_ELASTICITY` (specs/014, re-tested in specs/024) reads every in-loop
activator edge at fold^ε. Around a cycle of n edges the loop gain is ε^n, so
the same ε is a tiny pull on a 2-node loop and a massive one on TP53's
component (ε = 0.95 cost TP53 −57).

## The rule (`DS_LOOP_GAIN=g`, default 1 = off, byte-identical)

Each in-loop activator edge is read at fold^(g^(1/n)), with n the size of its
strongly connected component. A cycle's total gain is then about g, whatever
its length, and a positive loop driven by U settles at A = U^(1/(1−g)): finite
amplification instead of railing.

This is pinned in `test/test_loop_elasticity.jl` (+23 assertions):
- the exact root;
- equal amplification for 4-node and 8-node loops;
- acyclic networks untouched;
- a guard against combining it with `DS_LOOP_ELASTICITY`.

**Measured limitation.** Each trip round the cycle closes only (1 − g) of the
gap, so a near-1 gain converges slowly. At g = 0.99 the fixture reads 1.76
after 500 sweeps and 2.67 after 5,000. The solver's default budget is 500
(`DS_MAX_ITERS`), so a minimal gain needs a larger budget, and a larger budget
also changes the ungained knife-edge loops.

## Pre-registration (committed before any arm runs)

Canonical build `20260925-2326_d4f4f64`, current defaults (`DS_KO_AGG=mean`):

| arm | settings | paired against |
|---|---|---|
| `gain0.95` | `DS_LOOP_GAIN=0.95` | the new canonical baseline |
| `gain0.9` | `DS_LOOP_GAIN=0.9` | the new canonical baseline |
| `iters5000` | `DS_MAX_ITERS=5000` | the new canonical baseline (what the budget alone does) |
| `gain0.99_iters5000` | both | `iters5000` |

- **Adopt** only if held-out in-release net > +15 with p < 0.05, the gain is not
  one pathway, and experimental is no worse than −15.
- If several arms pass, prefer the gain closest to 1, because Adam asked for
  minimal.
- **Also reported:** the MH within-pathway error excess for positive loops,
  before and after, and how many solves converge.
