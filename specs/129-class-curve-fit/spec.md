# specs/129: one activating curve, chosen on measured FIT data (specs/125 stage 3a)

**Status:** registered 2026-11-03, before any arm runs. First fitting step of the learnable model (specs/125 §8,
stage 3). Nothing is adopted here: the result is a candidate for the one scoring at stage 5 (CPPA test portions,
Release 98, the Lun 2019 reserve, specs/127).

## Why, and what specs/128 already showed

specs/128 built response curves on edges by class (`DS_EDGE_CURVES`) and measured fixed settings. Re-read by
three-class balanced accuracy, the measure the Lun test uses (specs/127), rather than raw accuracy (post hoc, this
spec's research.md §0):

| specs/128 arm | experimental | korkut | curator held-out |
|---|---|---|---|
| `ec_off` (today) | 0.689 | 0.413 | 0.868 |
| `ec_sub90` (substrates, e = 0.9) | 0.702 | 0.413 | 0.857 |
| `ec_all75` (all, e = 0.75) | 0.667 | 0.419 | 0.800 |

So the raw korkut gain of `ec_all75` (+353) was mostly fewer change calls on a set that is 91% "no change"; by balanced
accuracy the measured data favour at most mild contraction. These numbers were seen before this registration; the grid
below contains settings near them, which is why the judge is cross-validated and not the best grid point.

## The judge (decided now; specs/125 §6, Adam's plan of 2026-11-02)

Learned parameters are judged on **measured FIT data by cross-validation**, then on the locked tests at stage 5.
Curator agreement (SELECT held-out, FIT tuning) is **reported, not a gate**: curators expect changes at full strength
at any distance, so a curator gate pulls every curve back to the identity whatever the measurements say (specs/128).

## What is fitted

One curve shared by the three activating classes (`input`, `catalyst`, `regulator_pos`), g(f) = f^h (1 + K^h) /
(f^h + K^h). Inhibiting classes stay at `divide` (specs/128 P2: today's inhibition is already the k = 9, m = 1 curve).
Per-class curves are a later step (3b), registered separately after this one.

**Grid** (12 points; elasticity at rest e = h K^h / (1 + K^h)):

| K \ h | 0.5 | 0.75 | 1 |
|---|---|---|---|
| inf | e 0.50 | 0.75 | **1.00 (identity, today)** |
| 30 | 0.42 | 0.70 | 0.97 |
| 9 | 0.38 | 0.63 | 0.90 |
| 3 | 0.32 | 0.52 | 0.75 |

Tables: `bench/analysis/129/cf_K<K>_h<h>.tsv`. No h > 1 (amplification): it raises positive loop gains above 1, which
specs/125 §3 rules out. A grid, not gradient fitting: with two shared parameters a grid is exact enough and needs no
gradient through the current non-smooth solver (specs/125 §4 stays the route for many parameters).

## Runs

Canonical frozen build `20261026-1349_8b6527e`; one arm per grid point plus control `cf_off` (no table), all at one
commit, `run_arm.sh --axis curator --axis experimental --axis korkut`, four at a time.

## Fitting data and objective

- **Experimental** (MP-BioPath, FIT; valid cases; perturbation = gene; RAF set aside for selection per specs/054 D1,
  reported both ways) and **korkut2015** (FIT; status `scored` cases; perturbation = drug set).
- Measured classes UP / NOT / DOWN; prediction as each bench writes it.
- **Objective:** the mean of the two datasets' three-class balanced accuracies.

## Cross-validation (`bench/analysis/score_129.py`, tested before use)

1. Perturbations are split into 5 folds per dataset (seeded, 129), each case going with its perturbation.
2. For each fold k: on the other four folds of both datasets, choose the grid point with the highest objective (ties: the
   larger e, then the larger K, so the identity wins a tie). That point predicts fold k's cases.
3. The CV prediction of each case is its fold's chosen setting; the CV objective is computed on the pooled held-out
   predictions.
4. **Gain** = CV objective − identity objective on the same cases; 95% interval from 10,000 bootstrap resamples of
   perturbations within each dataset (seed 129).
5. The **full-data choice** (step 2 on all folds) is the stage-3a candidate.

## Verdict (fixed now)

- **The curve earns its contraction** if the gain's interval excludes 0 above. The full-data choice is then the
  candidate the learned model carries to stage 5 (with stage 2's resting levels if they exist by then).
- Otherwise the identity stays the candidate: the learned model starts stage 5 equal to the frozen one in its curves.
- Either way nothing changes a default.

Also reported for every grid point: curator held-out and tuning balanced accuracy and right counts, McNemar against
control on each axis, converged solves, rails (specs/128's definition); for the full-data choice, the pathways moved.

## Predictions (before any arm runs)

- **P1.** `cf_off` reproduces `ec_off` (specs/128) case for case on all three axes.
- **P2.** The full-data choice is a mild contraction (e between 0.75 and 0.97), not the identity.
- **P3.** The CV gain is small (< 0.02) and its interval includes 0.
- **P4.** Every non-identity grid point loses curator held-out balanced accuracy against `cf_off`.
- **P5.** No grid point has fewer converged solves than `cf_off` (contraction lowers loop gains).
