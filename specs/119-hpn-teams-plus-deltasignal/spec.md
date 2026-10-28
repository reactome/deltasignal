# specs/119: does DeltaSignal improve the HPN-DREAM teams' own predictions?

**Status:** registered 2026-10-28, before any team's submission was converted to scores or combined with DeltaSignal.

Adam, 2026-10-28: look at how the best solutions run and see whether DeltaSignal can improve them; "ok" to the
HPN-DREAM step after specs/118 (korkut).

## Data (all already downloaded by Adam, specs/112; HPN-DREAM is FIT since specs/112)

- The 74 teams' sub-challenge 1A submissions (`SC1Asubmissions.zip`, syn bundle of 2026-10-26) and the aggregate
  ("wisdom of crowds") submission (`SC1AaggregateSubmission.zip`): one edge-score file per (cell line, stimulus).
- DREAMTools 1.2.4's per-team, per-context AUROCs (`dream8/D8C1/data/SC1A_aucs.json`), for validation only.
- The truth vectors (`SC1A_TrueDescVectors`), already scored in specs/112.
- DeltaSignal's frozen predictions (`specs/112-hpn-dream/hpn_predictions.tsv`, `ds_score`, NOISE rule of amendment 1)
  and reachability's (`reach_score`). They do not depend on the cell line.
Nothing from these files is copied into the repository: only aggregates are reported, and no team is named.

## A team's readout scores (DREAMTools' own rule, reimplemented)

As `HPNScoringNetwork`: parse each edge file (name fixes and ignored species as DREAMTools), |score| into an N x N
matrix, divide by the team's maximum over all 32 matrices, close it with the max-min Warshall recursion
path[i][j] = max(path[i][j], min(path[i][k], path[k][j])), and read mTOR's row. AUROC with ties counted half
(`hpn_dream_bench.auroc`).

**Validation, before anything is combined.** Each team's 30 official-context AUROCs must equal DREAMTools'
`SC1A_aucs.json` to 1e-6. Teams that do not are listed and left out; if more than 5 fail, the study stops.

## The combination (fixed here, not tuned)

Per context, each source's scores are turned into average ranks scaled to [0, 1] (ties share a rank). Combined score
= (team rank + DeltaSignal rank) / 2. Secondary: weight 0.25 on DeltaSignal. Controls, with the same rule:
- **Reach:** reachability's scores instead of DeltaSignal's.
- **Shuffled:** DeltaSignal's scores permuted among the readouts of each context (1,000 permutations, seed 119), to
  show what combining with any predictor of that shape gives.

## Measures

For every team and the aggregate: mean AUROC over the 30 official contexts, alone and combined. Reported: the best
team's change, the mean change over the top 10 teams, the number of the 74 teams improved, and the aggregate's change.
95% bootstrap over the 30 contexts (10,000 resamples, seed 119) for the best team's change.

**Hypotheses.**
- H1: combining with DeltaSignal raises the best team's mean AUROC (interval excludes 0).
- H2: DeltaSignal raises the top 10 teams' mean AUROC by more than reachability does.
- H3: the best team's gain from DeltaSignal is above the 95th percentile of its gain from shuffled DeltaSignal.
- H4: combining with DeltaSignal raises the aggregate submission's mean AUROC.
Each is reported whichever way it falls.

**What it is not.** The truth vectors were scored in specs/112; this is a development result on FIT data. No setting
of DeltaSignal changes, and the combination rule is fixed above, not chosen on the result.
