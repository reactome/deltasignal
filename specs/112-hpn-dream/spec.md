# specs/112: the frozen model on HPN-DREAM (the external test)

**Status:** pre-registered 2026-10-26, before any measured value was scored.

Adam, 2026-10-26: approved the HPN-DREAM breast cancer challenge as the external test ("if something is from a
competition and we can show that we can do better than them that would be perfect"), downloaded the data from
Synapse (syn1720047) after accepting its terms, and approved building the Reactome pathways that act on mTOR with
the frozen generator. The model is the freeze of specs/111 (tags `freeze-2026-10-26`); nothing in it changes for
this test.

## The data

Hill et al., Nature Methods 13:310 (2016). Reverse-phase protein arrays in four breast cancer lines (BT20, BT549,
MCF7, UACC812) under eight stimuli (EGF, FGF1, HGF, IGF1, Insulin, NRG1, PBS, Serum), 41–48 phosphoprotein readouts
per line. Sub-challenge 1A held out the mTOR inhibitor AZD8055: for each (cell line, stimulus) context, a team scored
every readout for being a descendant of mTOR, and was scored by AUROC against binary vectors of the readouts the
inhibitor changed (`SC1A_TrueDescVectors`), averaged over 30 official contexts. 74 teams took part.

The download (MD5s checked against Synapse's manifest) lives in `~/deltasignal-catalogs/external/hpn_dream/`,
outside both repositories, and is registered TEST in `bench/datasets.tsv`.

**Disclosed exposure.** While reading file headers, one truth vector (MCF7, EGF: 41 values, unlabelled) was printed
in full by mistake. The other 31 were not read. The model is frozen, so this cannot change a prediction; the
choices below that it could have influenced (networks, site map, scoring) are fixed by stated rules.

## Protocol (`bench/hpn_dream_bench.py`)

**Networks.** Every non-disease Reactome pathway with a reaction on an MTOR-containing entity (input, output,
catalyst or regulator, through complexes and sets), found in Neo4j by structure (`mtor_pathways.tsv`, 13
pathways; the disease pathways R-HSA-5674400 and R-HSA-9920951 are left out), built with the frozen generator as
the variant catalog `hpn112`. A network with no root form of MTOR is not used.

**Conditions.** In each network and for each stimulus:
- A: the stimulus ligands' root forms pinned at the bench's overexpression level (`DS_PERTURB_UI_UP`, 80×):
  EGF, FGF1, HGF, IGF1, INS (Insulin), NRG1; Serum all six at once; PBS no pin.
- B: A plus every root form of MTOR pinned at the knockout level (`DS_PERTURB_UI_DOWN`, 0), the drug bench's rule
  for an inhibitor (specs/021).
Pins follow the frozen bench's `pin_occurrences` (`DS_PIN_SCOPE=root_cycle`, `DS_PIN_ROOT_IGNORE=depletion`).

**Readouts.** `site_map.tsv` gives each readout its (gene, site) pairs, from the challenge's own antibody table
(`AbNames_CellLines.csv`) with four documented corrections: EGFR mature numbering +24 (Y1068 → Y1092), SRC chicken
numbering +3, MEK1 S217/S221 → S218/S222, "PDK1_pS241" is PDPK1, and the ERK antibody also reads ERK2 at T185/Y187.
A readout's value in a network is the mean over the unpinned nodes that are, or contain, a phosphorylated EWAS at
one of its sites (Reactome 97). Its score is |ln(B / A)|, the largest over the networks that hold it, and 0 where
none does. The model has no cell-line context: the four lines get the same scores.

**Scoring.** AUROC per context, ties counted half (equal to the challenge's ROC, which adds a point only where the
score changes, `test_hpn_dream_bench.py`); readouts whose truth is NaN are dropped, as the challenge did. The
headline is the mean over the 30 official contexts (BT20 Insulin and BT549 NRG1 left out, as in the final scoring).

**Baselines,** on the same contexts:
- constant score (AUROC 0.5);
- reachability: 1 if a readout node is reachable from an MTOR pin in any network, else 0;
- the 74 teams of the official leaderboard (`HPN-DREAM-Leaderboards-Official_v3.ods`, syn2289138), once Adam
  downloads it: our mean AUROC placed among theirs, and our mean rank over the 30 contexts if the file holds
  per-context AUROCs.

## Predictions

- H1. Before scoring: at least half of the readouts occur in at least one network.
- H2. DeltaSignal's mean AUROC over the 30 contexts is above 0.5.
- H3. DeltaSignal's mean AUROC is at least reachability's.
- H4. DeltaSignal places below the leaderboard median. Many readouts are absent from Reactome or from these
  networks and score 0, and the teams fitted their networks to the training data; the frozen model was fitted
  to nothing here.

## After scoring

The result is reported whatever it shows. The dataset then becomes FIT, and errors may be traced, but no change
made from them can be scored on it again.

## Amendment 1 (2026-10-26, after `predict`, before any truth vector was scored)

The predictions (`hpn_predictions.tsv`, committed at 69de560) give several readouts effects of 1e-15 to 1e-14
where the MTOR pin does not reach them: floating-point noise. Under AUROC such values would rank above true zeros,
breaking ties arbitrarily. An effect below 1e-6 (the solver's tolerance) now scores 0, in `score`; the
predictions file is unchanged. Readouts with a real effect are from 1.79 up.

Also recorded before scoring: H1 fails. 20 of the 48 readouts occur in at least one network (predicted: at least
half).

## Amendment 2 (2026-10-26, after H4, before any rank was computed): official-style rank and its interval

Adam, 2026-10-26: approved firming up the leaderboard placement. The official leaderboard file holds only each team's
mean AUROC, but DREAMTools 1.2.4 (the challenge's own scoring package) ships every team's AUROC per context
(`dream8/D8C1/data/SC1A_aucs.json`) and the ranking rule (`SC1A_ranking.get_ranking`): per official context, rank the
entries by AUROC (descending, ties at the lowest rank, pandas `method="min"`), then average the 30 ranks. No truth
vector is read in this step; it uses only `score.json` and the teams' AUROCs.

- **Validation, before any rank is reported.** The teams' AUROCs must reproduce the leaderboard: for every team, the
  mean over the 30 official contexts equals the sheet's `mean AUROC` to its printed precision (±0.0001 after
  rounding), and the mean rank computed among the 74 equals the sheet's `mean Rank` to ±0.001. Teams are matched by
  name; any unmatched or mismatched team is listed. If more than two teams fail, no rank is reported.
- **Rank.** DeltaSignal (and, separately, the reachability baseline) is added as a 75th entry with its 30 per-context
  AUROCs from `score.json`; reported: its mean rank, its place by mean rank, and its place by mean AUROC.
- **Interval.** A bootstrap over contexts: 10,000 resamples of the 30 official contexts with replacement (numpy
  seed 112), each recomputing every entry's mean AUROC and mean rank; reported: the 2.5–97.5 percentile interval of
  DeltaSignal's place on both measures. It does not resample readouts (the teams' per-readout scores are not in the
  file), so it understates the uncertainty.
- Prediction, written now: DeltaSignal's place by mean rank is between 8th and 16th, and the bootstrap interval of
  its place by mean AUROC spans at least five places.

## Amendment 3 (2026-10-27, after amendment 2's validation failed, before DeltaSignal's rank was computed)

Adam, 2026-10-27: "do what you think is best", on the proposal below. This changes a gate after it failed, so it is
recorded as a post-hoc change; it is written before DeltaSignal's or reachability's rank was computed or seen, and it
uses only the teams' figures and the sheet (research.md, amendment 2 diagnosis).

Amendment 2's check failed for reasons the diagnosis traced to how the sheet was produced, not to the per-context
AUROCs. The validation now reads:

- **AUROC:** the sheet truncates to four decimals: each team's computed mean must satisfy sheet ≤ mean < sheet + 0.0001
  (1e-9 slack).
- **Rank:** the sheet's mean rank may exceed ours by a whole number of thirtieths, 0 to 3 (±1e-6), never less; the
  official ranks were computed with entries later removed from the sheet. And the order must be identical: sorted by
  our mean rank, the sheet's `Final rank` never decreases.
- **Name:** `MLBHKKG` on the sheet is `ML_BHK_KG` in the per-context file.

The two-failure limit stands. Reported as in amendment 2. Because the removed entries may also have outranked
DeltaSignal in some contexts, its official-style mean rank may be worse by up to about 3/30 than computed; this is
reported beside the place. Amendment 2's prediction stands unchanged.

**Correction to amendment 3 (same day, before any rank was computed).** The first run under amendment 3 failed 49
teams, because the ±1e-6 written above is finer than the sheet prints mean rank (three or four decimals, e.g. 7.7333
for 7 + 22/30). The rank check is now: the sheet's mean rank equals ours plus k/30 for some k in 0–3, within ±0.001
(amendment 2's precision). DeltaSignal's and reachability's ranks were not computed by the failed run.

## Amendment 4 (2026-10-28, after scoring): two input errors found by code review, corrected and re-scored

A review of everything since the freeze found two errors in this test's inputs, not in the frozen model.
- **Site map.** Several antibodies were mapped to the nominal residue on every paralog, though paralogs carry the
  homologous site at a different number, so those entries could never match a curated phospho-form: AKT2 and AKT3
  (S474/T309, S472/T305), GSK3A (S21, the residue matching GSK3B S9), PRKAA1 (T183), ACACB (S222), and the
  pan-PKC hydrophobic motif (PRKCA S657, PRKCB S660/S661, PRKCD S664, PRKCE S729, PRKCH S675, PRKCQ S695). The
  corrected residues were checked in the Release 97 database (each except PRKCH S675 and ACACA S79 is curated). The
  JNK note is corrected (residue 185 is a tyrosine; nothing was added).
- **Unstimulated baseline.** Where no stimulus ligand is in a network (always for PBS), condition A was not solved and
  B was compared with a constant 0.01. The 1e-15 to 3.6e-14 "effects" amendment 1 recorded come from that, not from
  unreachable readouts. A is now always solved.

Registered before the re-run: `predict`, `score` and `rank` are re-run on the same build at the commit that makes these
changes, with every other setting unchanged. The original 0.693 (10th) stays the scored-once result of the frozen model
and is reported first; the corrected figures are reported beside it as a correction. HPN-DREAM is FIT since the first
scoring, so the re-run is not a second test, and whichever way it moves, it changes no rule.
