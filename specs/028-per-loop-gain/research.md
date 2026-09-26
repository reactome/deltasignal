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
strongly connected component. A cycle through k such edges has gain g^(k/n).
That is g only when the cycle spans its whole component, as in the test
fixtures, where a positive loop driven by U settles at A = U^(1/(1−g)).
Inside a large component with short cycles the pull is much weaker (see the
Result).

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

## Result (arms at `30dd48c`; base `261c94a/baseline`, which has `DS_KO_AGG=mean`)

| arm | paired against | held-out net | pathways / perturbations | tuning | experimental |
|---|---|---|---|---|---|
| `gain0.95` | baseline | +34 (34 / 0), p 1e-10 | 1 / 2 (IFN-γ) | −1 | 0 |
| `gain0.9` | baseline | +40 (40 / 0), p 2e-12 | 2 / 4 (IFN-γ +34, Fanconi +6) | −1 | 0 |
| `iters5000` | baseline | **−89** (3 / 92), p 7e-24 | 3 / 17 (DSB Repair −74, DSB Response −16) | 0 | 0 |
| `gain0.99_iters5000` | `iters5000` | +34 (34 / 0) | 1 / 2 (IFN-γ) | −1 | 0 |

- **Not adopted.** The gain is one pathway (Interferon-γ), which fails the
  pre-registered concentration condition. It is also harmless: 0 cases broken
  outside tuning's −1, and 0 on experimental.
- **The budget result is the more important finding.** Raising
  `DS_MAX_ITERS` from 500 to 5,000 alone costs 89 held-out cases, 74 of them in
  DNA Double-Strand Break Repair. **The default model's answers in DSB depend
  on the iteration budget.** Its positive loops are still on the specs/014
  knife-edge, and more sweeps rail them.
- **Why the per-loop gain does not fix that.** It divides g across the size of
  the strongly connected component. DSB's component is hundreds of nodes, but
  its cycles are short, so each actual cycle keeps a gain of about
  g^(L/n) ≈ 1.
- **The pull has to be per cycle.** For example, scale it by the length of the
  shortest cycle through each edge rather than by the component's size. That
  is the next design; it is not attempted here.


## Corrections after review of PR #74

- **Convergence was not measured.** "How many solves converge" was
  pre-registered but not recorded. The benchmark now logs `Converged: x of y
  solves`. The gain arms ran at the default 500-sweep budget.
- **The iteration-budget result is a caveat on the canonical numbers
  themselves.** DSB's answers are budget-dependent (`iters5000` −89). The
  mechanism, "more sweeps rail the loops", is inferred from specs/014, not
  traced here.
- **MH positive-loop excess** on `261c94a/baseline`: all-positive loops
  +11.3% within pathway (34 pathways); loops with an inhibition +4.2% (16);
  giant loops +13.1% (9); welded loops +16.1% (6).
