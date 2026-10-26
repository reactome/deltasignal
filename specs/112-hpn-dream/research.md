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
- H4 is open: it needs the official leaderboard (syn2289138), which is not yet downloaded.

## Reading

The frozen model, fitted to nothing in this dataset and given no cell-line context, ranks the readouts the mTOR
inhibitor changed above those it did not in every official context. Most of its signal is structural: which
readouts the curated networks connect to mTOR (reachability alone gives 0.667). The model adds to that by
grading effects (full suppression of S6K, S6 and 4EBP1 above a partial change in GSK3 and p27), which places the
canonical mTORC1 substrates first. Its ceiling is coverage: 28 of the 48 readouts are not in any network that acts on
mTOR and tie at 0.

The comparison with the 74 teams waits for the leaderboard. The paper reports a mean AUROC of 0.82 for a weighted
combination of the top methods and an aggregate prior network (Hill et al. 2016, Supplementary Fig. 9). That
combination is not a single team.
