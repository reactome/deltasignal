# specs/112 results

Frozen model (`freeze-2026-10-26`; solver identical to 05dc867), variant catalog `20261026-1643_8b6527e_hpn112`
(13 of 13 pathways, frozen generator), bench at 4264c58. Predictions at 69de560 (`hpn_predictions.tsv`), written
before any truth vector was scored; scored once (`score.json`). Solves: 110; every rule at its frozen default.

## Coverage (before scoring)

11 of the 13 networks hold a root form of MTOR. 20 of the 48 readouts occur in at least one of them; the rest
score 0. Readouts with a real effect under the MTOR pin: p70S6K T389, S6 S235/S236 and S240/S244, 4EBP1 T37/T46,
AKT S473 and T308 (mTOR signalling and the PECAM1/KDR network; full suppression, |ln| 13.8), PRAS40 T246 (6.25),
GSK3 S9/S21 and p27 T157 (PIP3 activates AKT; 6.22 under EGF, FGF1, HGF, NRG1 and Serum, 1.79 under IGF1,
Insulin and PBS). AMPK T172 is reachable but does not move.

## Result: mean AUROC over the 30 official contexts

| method | mean AUROC (30 official contexts) | all 32 contexts |
|---|---|---|
| **DeltaSignal (frozen)** | **0.693** | 0.689 |
| reachability from the MTOR pin (same networks) | 0.667 | 0.664 |
| constant score | 0.500 | 0.500 |

DeltaSignal is above reachability in 29 of the 30 official contexts and below it in one (BT20 NRG1, 0.660 vs
0.684). Per context it ranges from 0.541 (BT549 EGF) to 0.859 (MCF7 FGF1); by cell line, MCF7 is best (0.66–0.86)
and BT549 worst (0.54–0.74).

## Predictions

- H1 fails (recorded before scoring): 20 of 48 readouts in a network, not half.
- H2 holds: 0.693 > 0.5.
- H3 holds: 0.693 ≥ 0.667, in 29 of 30 contexts.
- H4 fails, in the model's favour: 0.693 is above the leaderboard median (0.532); see below.

## Reading

The frozen model, fitted to nothing in this dataset and given no cell-line context, ranks the readouts the mTOR
inhibitor changed above those it did not in every official context. Most of its signal is structural: which
readouts the curated networks connect to mTOR (reachability alone gives 0.667). The model adds to that by
grading effects (full suppression of S6K, S6 and 4EBP1 above a partial change in GSK3 and p27), which places the
canonical mTORC1 substrates first. Its ceiling is coverage: 28 of the 48 readouts are not in any network that acts on
mTOR and tie at 0.

## Leaderboard (H4)

`HPN-DREAM-Leaderboards-Official_v3.ods` (syn2289138; downloaded by Adam 2026-10-26 in the challenge's file bundle,
every MD5 matching its Synapse manifest; kept outside both repositories), sheet SC1A: 74 teams, each with a mean
AUROC over the same 30 official contexts.

| | mean AUROC |
|---|---|
| best team | 0.782 |
| 10th team | 0.677 |
| team median | 0.532 |
| team mean | 0.568 |
| worst team | 0.394 |

DeltaSignal's 0.693 is above 65 of the 74 teams and below 9, so it would place **10th of 75** by mean AUROC.
Reachability on the same networks (0.667) would place 11th: most of the standing comes from the curated wiring
itself, and the model's grading adds one place.

Caveats:
- The challenge's official rank is by mean rank over the contexts, not mean AUROC. The leaderboard does not hold
  per-context AUROCs, so our mean rank was not computed. The teams' submissions (networks, `SC1Asubmissions.zip`)
  are in the bundle, and re-scoring them would give it, but that step was not pre-registered.
- The teams fitted their networks to the challenge's training data, in the same cell lines; the frozen model was
  fitted to nothing here and has no cell-line context. The comparison is not like for like in either direction.
- The paper reports a mean AUROC of 0.82 for a weighted combination of the top methods and an aggregate prior
  network (Hill et al. 2016, Supplementary Fig. 9). That combination is not a single team.

## Amendment 2: validation failed (no rank reported)

Run at bench 6afc0dc (`hpn_dream_bench.py rank`, committed before the run) on DREAMTools 1.2.4's
`SC1A_aucs.json` (74 teams × 32 contexts) and the official SC1A sheet. Under the pre-registered criterion, 33 of
the 74 teams fail and one name is unmatched, so **no rank is reported**. DeltaSignal's place was not computed:
the run stopped at validation. Output (team-level, kept outside the repo):
`~/deltasignal-catalogs/external/hpn_dream/leaderboard_20261026/analysis/rank_validation_failed.json`.

Diagnosis (from the teams' figures only):
- **Mean AUROC.** The sheet truncates to four decimals rather than rounding. Truncated, all 74 teams match exactly;
  rounded, 44 do.
- **Mean rank.** The sheet's value minus ours is always a whole number of thirtieths (0, 1, 2 or 3 contexts) and
  never negative. That is what you get when the official ranks were computed with a few more entries present, which
  beat some teams in some contexts, and those entries were then left off the sheet. DREAMTools' `remove_users` drops
  duplicate and withdrawn accounts, which fits. The order of all 74 teams is identical to the sheet's final ranks.
- **Name.** `MLBHKKG` on the sheet is `ML_BHK_KG` in the per-context file.

So the per-context AUROCs are the leaderboard's own, and the failure is in how the check compared them. Changing the
check now would be a change after the gate failed, so it waits for Adam.

## Amendment 3: the rank

Adam approved changing the check ("do what you think is best"). Amendment 3 was committed (50d8a41) before any rank
was computed. Three runs followed, all on the same inputs; the per-team output stays outside the repo
(`builds/20261026-1643_8b6527e_hpn112/results/03e602e/hpn112/rank.json`):

1. b42d60f: 49 teams failed, because the ±1e-6 rank tolerance as written is finer than the sheet's printed
   precision. No rank was computed. The correction (±0.001, amendment 2's precision) was committed (7a5b1ff)
   before the next run.
2. 7a5b1ff: one failure, within the limit of two, so the rank was reported. The failure was the order check: two
   teams near the bottom (places 53 and 54) tie on our mean rank, and the sheet breaks the tie through the removed
   entries.
3. This commit: the order check accepts teams tied on our mean rank in either order. This was **changed after the
   rank was seen**. It changes no number and does not change whether the gate passes; it only removes the spurious
   entry. Validation failures: 0.

| entry (75 with the 74 teams) | mean AUROC | mean rank | place by mean rank | place by mean AUROC | 95% bootstrap interval, both measures |
|---|---|---|---|---|---|
| DeltaSignal | 0.693 | 16.2 | **10th** | 10th | 7th–13th |
| reachability | 0.667 | 19.6 | 11th | 11th | 9th–19th |

The bootstrap uses 10,000 resamples of the 30 contexts, numpy seed 112. It does not resample readouts.

- **Prediction:** place by mean rank between 8th and 16th: **holds** (10th). The AUROC-place interval spans at
  least five places: **holds** (7th–13th, seven places).
- **Removed entries.** These could raise DeltaSignal's official mean rank by up to 3/30, to at most 16.3. The next
  entry below is at 18.7, so the place is unchanged.
- **Reading.** By the challenge's own ranking rule the frozen model places 10th of 75. Its interval stays inside the
  top 13 and does not reach the top 6. Reachability alone, the wiring without the solver, places 11th. It has a wider
  interval, reaching 19th. So the solver's contribution is real but modest. The caveats in the Leaderboard section
  still apply.
