# specs/131: learning a strength per edge (specs/125 stage 3b), feasibility first

**Status:** registered 2026-11-03, before any code. Adam, 2026-11-03: "ok" to trying per-edge learning, small and sparse
first, after specs/128–130 found that every way of making the model call fewer or smaller changes trades correct
change calls for correct no-change calls. What the model lacks is which readouts respond; a strength per edge, learned
where data reach, is the part of specs/125 that could supply it.

## Approach: a differentiable surrogate proposes, the real solver judges

The solver has hard switches (rules, caps, iteration fallbacks) and gives no gradient (specs/125 §4). Instead of
rebuilding it first, a **linear-response surrogate** of the propagator is fitted, and fitted weights are only ever
judged in the real solver.

- **Surrogate.** On each network's `logic_network.csv`, log2-fold changes propagate as δ = A δ + s, where s holds the
  pins and A has one entry per edge: +1 for an AND input, catalyst or positive regulator; −1 for a negative regulator
  or depletion; 1/k for each of a node's k OR in-edges (output, variant pool and split, specialization). Cofactor,
  drug and pathogen nodes (the bundle's tables) are held at 0. Pins: over-expression log2 80, knockout log2 0.01 (the
  bottom of the 100× scale). Readout value: the mean over its nodes, as the benches read it; call UP / DOWN beyond
  log2 1.15 / log2 0.85. Solved by damped iteration, δ ← clip(A δ + s, ±log2 100), to a fixed point or 500 sweeps.
- **Per-edge weights.** A → W ∘ A, W_e ≥ 0, initial 1. The solver's derived edges (pools, carriers, twins, RNA supply)
  are not in the surrogate; S1 measures what that costs.

## Stages and gates (each gate fixed now)

**S1 — does the surrogate agree with the solver?** On the FIT perturbation cases (experimental, korkut; frozen build
`20261026-1349_8b6527e`), surrogate calls at W = 1 against the frozen solver's calls (`cf_off`, specs/129); no measured
value is read. Gate: **call agreement ≥ 85% on each dataset, and ≥ 70% on the cases where the solver calls a change.**
If S1 fails, stop: gradients of a surrogate that disagrees with the solver would steer weights for a different model.
Reported: agreement, the confusion table, the pathways where it fails, and three traced disagreements.

**S2 — is there anything to learn, in the surrogate?** Fit W on FIT measured calls (experimental and korkut) with
5-fold CV by perturbation (as specs/129): a smooth loss (cross-entropy on the surrogate's readout value through a
logistic link at the call thresholds) plus λ Σ (log W_e)², λ chosen by inner CV from {0.01, 0.1, 1, 10}. Only edges on a
route from some FIT perturbation's pins to some FIT readout are free; the rest stay 1. Gate: **CV balanced accuracy of
the fitted surrogate exceeds the W = 1 surrogate's, interval (bootstrap over perturbations, 10,000, seed 131) above
0.**

**S3 — does it carry to the real solver?** Requires a per-edge weight table in the solver (g(f) = f^{W_e} on that slot,
the specs/128 machinery per edge; registered as an amendment when built). The S2 weights fitted per fold are run in the
solver on that fold's held-out perturbations. Gate as S2, in the solver. Only if S3 passes is the fitted model a
candidate for stage 5.

## Predictions

- **P1.** S1 agreement on all cases is between 75% and 90%: the derived edges and the caps matter, but most calls
  follow the curated wiring.
- **P2.** If S1 passes, S2's fitted weights move most on edges near the perturbed targets (within three hops) and
  leave long routes weaker than 1 on average.
- **P3.** S2's CV gain, if any, comes mostly from korkut (more perturbations sharing routes).

Nothing here changes a default or reads a SELECT or TEST set.
