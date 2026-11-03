# specs/131 research

## S1: surrogate against the solver (2026-11-03)

`bench/surrogate_131.py` run through the benches' own pins, readouts and calls (`bench/run_surrogate_131.py`: the
server parses each network, every solve is the surrogate's), frozen build, against `cf_off` (specs/129, solver
`d210ae3`). Predicted calls only; no measured value read.

| dataset | cases | agreement | on solver change calls | surrogate converged | agreement where converged / not |
|---|---|---|---|---|---|
| experimental | 845 | **71.1%** | **68.1%** | 146 of 244 solves (326 cases) | 88% / 61% |
| korkut | 11,062 | **90.7%** | **90.3%** | 1,208 of 1,493 solves (7,604 cases) | 96% / 79% |

**S1 fails** (the gate needed ≥ 85% / ≥ 70% on each dataset). Experimental: 188 cases flip direction against the
solver, 179 of them in surrogate solves that never settle; by pathway, agreement is lowest in RAF (0.59), TP53 (0.60),
Cell Cycle Checkpoints (0.66), PIP3 and WNT (0.69), and highest in S Phase (1.00), Mitotic Prophase (0.96) and HDR
(0.92).

Traced: APC and AMER1 knockouts in WNT read about 10× at cytosolic β-catenin in the solver (the destruction complex
is gone) and 0.01–0.07× in the surrogate. In the solver the destruction complex is a release pool (specs/086) solved
by the cycle balance; in the surrogate the same cycle has loop gain ≥ 1, so the knockout's fall circulates and rails
the whole cycle at the bottom. The same holds for the RAF switch and the TP53 pools: the solver's loop rules (pools,
carrier loops, twins), which the surrogate does not have, are what decide those cases, P1's "derived edges matter",
concentrated in loops.

P1 (agreement between 75% and 90%) holds on the two datasets together (89.4% of 11,907 cases) and fails on each alone
(71% and 91%).

## Amendment 1 (2026-11-03): S2 on korkut only, registered before any fitting

S1 failed as registered, and the registered consequence (stop) stands for experimental: weights fitted on a surrogate
that disagrees with the solver in 29% of its cases would steer a different model. On korkut the surrogate passes the
S1 gate (90.7% / 90.3%), so S2 runs there alone, unchanged otherwise:

- Data: korkut scored cases; perturbation = drug set (as specs/129); 5 folds, seed 131.
- The loss uses only surrogate solves that converge at the current weights (non-converged solves cannot be
  differentiated at a fixed point); CV scoring uses every case, the surrogate's call as it comes.
- Gradient by implicit differentiation at the fixed point, over the free (unpinned, unclipped) nodes.
- Free weights: edges on a route from a training fold's pins to a korkut readout node.
- Gate as registered: CV balanced accuracy of the fitted surrogate above the W = 1 surrogate's, interval (10,000
  bootstrap resamples of drug sets, seed 131) above 0. S3 follows only if S2 passes.

Caveat recorded now: korkut's drug sets share single drugs (`901+AK`, `901+NT`), so a held-out drug set usually shares
a drug with training sets. A pass is then evidence for "unseen combination", not "unseen drug".

## S2 result (2026-11-03): fails

Run as amendment 1 registered it: `bench/fit_131.py` (code `780d187`) on the korkut problems collected by
`bench/collect_131.py` from the frozen build (1,481 solves, 11,062 scored cases, 53 networks); output
`~/deltasignal-catalogs/analysis/131/s2_korkut.json` and `.log` (outside the repository). 40 Adam epochs per fit,
learning rate 0.05, 5 outer folds by drug set, λ by inner 3-fold CV.

| | balanced accuracy (3-class) |
|---|---|
| surrogate, W = 1 | 0.3986 |
| surrogate, fitted per fold, on held-out drug sets | 0.4035 |
| **CV gain** | **+0.0049, 95% interval [−0.0005, +0.0104]** |

**The gate fails**: the interval includes 0. S3 is not built, as registered.

| fold | λ chosen | inner-CV spread over λ | edges free / moved | log-weight of moved edges, 5% / 50% / 95% |
|---|---|---|---|---|
| 0 | 0.1 | 0.4036–0.4049 | 27,072 / 2,402 | −1.24 / −0.03 / 0.29 |
| 1 | 0.01 | 0.4060–0.4075 | 27,072 / 2,403 | −1.57 / −0.53 / 0.52 |
| 2 | 1 | 0.4001–0.4019 | 27,072 / 2,407 | −1.18 / −0.01 / 0.13 |
| 3 | 1 | 0.4096–0.4104 | 27,072 / 2,331 | −1.16 / −0.01 / 0.11 |
| 4 | 1 | 0.3984–0.4001 | 27,072 / 2,407 | −1.21 / −0.01 / 0.14 |

What it says, plainly:

- **Learning edge strengths changes almost nothing.** The fitted surrogate is half a point better on drug
  combinations it did not see, and that half point is within what reshuffling the drug sets produces. Four values of
  λ, from almost free (0.01) to pulled hard toward weight 1 (10), score within 0.002 of each other in every fold: how far the
  weights are allowed to move hardly matters to the calls.
- **Only about 2,400 of 27,072 eligible edges get any signal**: the rest lie on routes from the pins that reach no
  readout in a solve that converges. The data reach a small part of the network.
- **What moves, mostly weakens.** Of the moved edges, the 5% tail sits near log-weight −1.2 (weight ≈ 0.3) in every
  fold, the upper tail at +0.1 to +0.5. The fit turns some routes down, the specs/128–129 direction again (fewer change
  calls), and gains about as little from it.
- **Both the fixed and the fitted surrogate sit near the floor.** Always-NOT scores 1/3 by balanced accuracy; the
  surrogate scores 0.40.

Predictions: **P2** (weights move most within three hops of the targets; long routes weaker on average) was not
measured: the fitter records the distribution of moved weights, not their distance from the pins, and with the gate
failed, adding that measurement now would be a post-hoc look at a null fit. The median moved weight is ≈ 1 in four
folds and 0.59 in the one where λ = 0.01, so there is no general weakening of routes. **P3** (gain mostly from korkut)
cannot be judged: amendment 1 left korkut the only dataset.

A timing trial run before the full fit (3 epochs, `--quick`; not a result) gave +0.005 [−0.0009, 0.011], the same as
the registered run: 37 more epochs bought nothing.

### Conclusion

On the one FIT dataset where a differentiable copy of the model agrees with the solver, edge strengths fitted to
measured calls do not carry to unseen drug combinations, even though such a set usually shares a drug with the
training sets (amendment 1's caveat: a pass would have been the easier claim). With specs/128–130, this is the fourth
arm since specs/128 in which fitting parameters of the frozen wiring to FIT data gains nothing that holds out.
The limit is not the number of free parameters: per-edge weights are the most flexible form registered, and they
gain no more than one shared curve did. What the model is missing is not a strength on the curated edges.

The surrogate itself (S1) is recorded as a tool: it reproduces the solver on 91% of korkut cases and 71% of
experimental ones, and the difference is the loop rules.
