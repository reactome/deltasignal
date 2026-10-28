# specs/116 research

Run 2026-10-27: `bench/analysis/train_116.py` (tested in `bench/test_train_116.py`) on the specs/114 top-level case
logs (build `20261027-1346_8b6527e_top114`, code `a47898d`, `s114_top_kb`). Hop tables from
`bench/analysis/korkut_hops.py` and `bench/analysis/decryptm_hops.py`, in `~/deltasignal-catalogs/analysis/116/`;
output `train_116.json` there. Balanced accuracy (bal) is the registered metric; intervals are 95% bootstrap over
held-out drugs for trained minus the arm.

| evaluation | cases | drugs | trained | frozen | reach | topology only | always-NOT | trained − frozen |
|---|---|---|---|---|---|---|---|---|
| korkut, leave one drug out | 33,488* | 12 | 0.389 | 0.396 | 0.387 | 0.384 | 0.333 | [−0.048, +0.036] |
| decryptM, leave one drug out | 5,096 | 18 | 0.431 | 0.365 | 0.382 | 0.422 | 0.333 | [−0.016, +0.130] |
| train korkut → decryptM | 5,096 | 18 | 0.400 | 0.365 | 0.382 | 0.396 | 0.333 | [−0.047, +0.091] |
| train decryptM → korkut | 17,836 | 12 | 0.397 | 0.396 | 0.387 | 0.376 | 0.333 | [−0.031, +0.025] |

\* a combination is tested once per drug it contains (registered).

Trained minus topology only: korkut [−0.004, +0.025], decryptM [−0.009, +0.025], korkut → decryptM [−0.021, +0.028],
decryptM → korkut [+0.006, +0.032].

Macro-F1 of the trained layer is 0.17-0.22 against the frozen model's 0.36: the balanced-accuracy rule calls DOWN in
many cells where NOT dominates, trading NOT precision for recall on the rare classes.

**Hypotheses.**
- **T1 fails.** The trained layer does not beat the frozen model in leave-one-drug-out on either dataset (both intervals
  include 0; on korkut it is lower, 0.389 vs 0.396).
- **T2 fails.** It beats the topology-only layer only in one of four evaluations (decryptM → korkut, +0.006 to +0.032),
  not in either registered leave-one-drug-out.
- **T3 fails.** Neither transfer interval excludes 0.

**Reading.** The frozen model's predicted fold and the hop distance do not hold the information needed to say which
readouts move. Re-weighting them cannot recover it, which agrees with specs/115: a change call is right about as often
as the base rate. Training has to add information the model does not have: which readouts respond to drugs in general
(learnable from other drugs without seeing the held-out one), how strongly the drug inhibits its target, and the
cell's context. Those are separate registrations.
