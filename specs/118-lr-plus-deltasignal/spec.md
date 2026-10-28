# specs/118: does DeltaSignal improve the best published korkut predictor?

**Status:** registered 2026-10-28, before the baseline was rebuilt or any combined prediction was made.

Adam, 2026-10-28: "shouldn't we look at how the best solutions run and then see if deltasignal can help improve it?"
Then "ok" to starting with korkut and the linear regression.

## The best published method

Long et al. 2024 (BMC Bioinformatics, "Causal models and prediction in cell line perturbation experiments") compare
methods on the korkut2015 data as the CellBox authors supply it (the files we hold, CellBox commit 6bed6b94). With each
drug held out, plain linear regression (LR) scores a mean correlation of 0.784 and CellBox 0.780. Their LR:
- X = R·D: X is the 87 responses of a condition (82 proteins, 5 phenotypes), D its 12 drug values.
- R is fitted by unregularised least squares.
- The held-out drug's column is 0. LR therefore learns nothing about the held-out drug; in a combination it predicts the
  partner drug's effect, and for the drug alone it predicts no change.
- Metric: for each held-out drug, the Pearson correlation between measured and predicted values over all its test
  conditions and the 87 responses, averaged over the 12 drugs.

That gap, a drug the regression has never seen, is where prior knowledge can help, and DeltaSignal's prediction
is made without any data.

## Data

korkut2015 (FIT): `expr.csv` (89 conditions x 99 nodes; the 87 response columns are used, the 12 activity columns are
not), `pert.csv` (the last 12 columns are D, CellBox's drug encoding; Long et al. say "drug concentrations"),
`expr_index.txt` (drug set per condition, read with `korkut_bench.conditions`). DeltaSignal: the specs/114 top-level case
log (build `20261027-1346_8b6527e_top114`, code `a47898d`, `s114_top_kb`, `korkut_cases.tsv`).

## Features (per condition and response, from the drug set; doses share a drug set)

- **DS (primary):** the mean, over the pathways where the case is scored, of DeltaSignal's call as a number (down −1,
  not 0, up +1). 0 where no pathway scores it (phenotypes, excluded or unreached readouts, trivial cases).
- **DS-fold (secondary):** the same mean over log2(pred_ui), each clipped to [−6.64, 6.64] (100x).
- **Reach (control):** the same mean over reachability's call (down −1, not 0).

## Models (fitted on the training conditions of each fold)

- **LR:** X ≈ R·D, as above.
- **LR + F:** X ≈ R·D + c·F, F one of the features and c one scalar shared by all responses. The fit is exact least
  squares: c = ⟨M F, M X⟩ / ⟨M F, M F⟩ with M the projection off D's column space, then R = D⁺(X − cF).
  c is learned from the training drugs only; the held-out drug's F is used at test time.
- **LR + F per response (secondary):** one c per response.

## Evaluation

Leave one drug out as Long et al.: test = every condition containing drug j, train = every condition without it.
Primary metric: Long et al.'s mean per-drug correlation. Secondary: balanced accuracy over down / not / up on the 82
proteins, each test condition called at log2(1.5) (bench threshold) for both measured and predicted values. 95%
bootstrap over the 12 held-out drugs (2,000 resamples, seed 118) for differences of the mean per-drug correlation.

**Reproduction check (not a hypothesis).** Our LR must give 0.784 ± 0.03. If it does not, D is replaced by the doses in
`expr_index.txt` (the drug's concentration, 0 if absent) and the check repeated; the version that reproduces is used
for everything. If neither reproduces, the comparison is reported but flagged as against a non-reproduced baseline.

**Hypotheses.**
- H1: LR + DS beats LR (interval excludes 0).
- H2: LR + DS beats LR + Reach (DeltaSignal's arithmetic beats the wiring alone).
- H3: on the conditions with the held-out drug alone (where LR predicts no change), the correlation of LR + DS with the
  measurement is above 0 (bootstrap over drugs).
Each is reported whichever way it falls. CellBox itself is not rebuilt here (LR matches it); that can follow.

**What it is not.** korkut2015 is FIT: a pass is a development result. No setting of DeltaSignal changes.
