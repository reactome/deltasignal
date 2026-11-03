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
