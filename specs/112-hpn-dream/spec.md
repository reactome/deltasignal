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
