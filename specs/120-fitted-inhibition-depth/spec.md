# specs/120: how hard a drug inhibits its target, fitted from the other drugs (training round 3)

**Status:** registered 2026-10-28, before any arm at a depth other than 0 was run on either dataset.

Adam, 2026-10-28, on specs/119: "nothing until the paper ... continue". Still no Reactome Release 98 and no new
datasets.

**Why.** The frozen model pins a drug's target at 0: complete inhibition. specs/115 found that 66-78% of the phospho
errors are false change calls (the model predicts a change, nothing was measured), and on korkut2015 the model scores
below always-NOT. A drug at the doses measured may inhibit its target only partly. The depth is one number, and it
can be learned from the other drugs without seeing the held-out one. This is a fitted parameter of DeltaSignal
itself, not a layer on top of it (specs/116-118).

## Arms

The specs/114 top-level build `20261027-1346_8b6527e_top114`, the bench's `DS_PERTURB_UI_DOWN` (the UI value a drug
target is pinned at), every other setting the frozen default. Solver and bench code are unchanged since `a47898d`
(`git diff a47898d HEAD -- src bench/korkut_bench.py bench/decryptm_bench.py bench/benchmark_vs_mpbiopath.py` is
empty).

| depth (pinned UI) | target activity | arm |
|---|---|---|
| 0 | none (frozen) | `s114_top_kb` (code `a47898d`, already run) |
| 0.1 | 10% of normal | `d120_010` |
| 0.25 | 25% | `d120_025` |
| 0.5 | 50% | `d120_050` |

Axes korkut and decryptm (kinobeads map, the specs/114 case and phospho files). The grid is fixed here and not
extended after seeing results.

## Evaluation

Cases scored in all four arms, measured call not `conflict`; the counts dropped are reported. Leave one drug out
within each dataset, as specs/116-117 (korkut: test every set containing the drug, train on every set without it).
For each held-out drug, the depth with the highest balanced accuracy of the model's call on the training cases is
chosen (ties: the deeper inhibition, i.e. nearer the frozen model), and the held-out cases are predicted at that depth.
Balanced accuracy over down / not / up (primary), macro-F1 and the false-change share reported; 95% bootstrap over
held-out drugs (2,000 resamples, seed 120). Also reported for reference: each fixed depth on all drugs (in-sample,
not a result), reachability and always-NOT.

**Hypotheses.**
- H1: the fitted depth beats the frozen model on korkut2015 (interval for fitted − frozen excludes 0).
- H2: the same on decryptM.
- H3: with the fitted-depth call in place of the frozen call, specs/117's R2 layer (responsiveness + DeltaSignal)
  improves on decryptM, where the model added information (interval excludes 0).
Each is reported whichever way it falls.

**What it is not.** The frozen model's default stays 0 whatever the result; a change of default is Adam's decision.
Both datasets are FIT, so a pass is a development result. HPN-DREAM is not used.
