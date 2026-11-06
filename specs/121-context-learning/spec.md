# specs/121: cell context and learned weights (project plan)

**Status:** registered 2026-10-28, before any CPPA measurement was read. Only the treatment design (drugs, cell lines,
doses, times) and the antibody names have been read. Adam approved the plan on 2026-10-28 ("I approve this plan").

**Why.** Three training rounds (specs/116-118, 120) and the error anatomy (specs/115) agree: the frozen model gets the
direction of a response right near the drug's target, but cannot tell which readouts move in a given cell. Reactome's
wiring does not hold that information. Two things could supply it: what the cell expresses, and how strong each curated
link is in practice. Adam, 2026-10-28: "I am fine with turning this into a larger project with expression data ... I
am fine with the idea of us needing to learn weights on edges. what ever is the proper way of doing this."

## Data

| dataset | role | use |
|---|---|---|
| korkut2015, decryptM, HPN-DREAM, LINCS P100 | FIT (already scored) | training and development |
| curator benchmark | SELECT | no-harm guard on every model change, as now |
| **CPPA v1.0** (Zhao et al. 2020, MD Anderson) | **split below** | training portion FIT; test portions TEST, scored once |
| **DepMap 24Q2 expression** (CC BY 4.0) | covariate | per-cell-line protein-coding expression; no response values |

CPPA is used under its terms of use, accepted by Adam for academic use; its data and anything derived from it stay
outside this repository. Only the split script, its seed and the hash of its output are committed. From DepMap's model
table only the identifier and name columns are read.

## The split (fixed now, before any value is read)

`bench/analysis/split_121.py` reads only the CPPA treatment summaries and writes the split; seed 121.
- **Units.** Samples with exactly one compound. Drug and cell-line names are compared case-insensitively.
- **Seen drugs.** A CPPA drug that is also a drug of korkut2015, decryptM, HPN-DREAM or LINCS P100 is never a test
  drug (it is training only).
- **Test cell lines.** 25% of CPPA cell lines, drawn at random.
- **Test drugs.** 40% of the drugs that are not seen drugs, drawn at random.
- **Training** = training drugs (seen drugs included) in training cell lines. **Test A, unseen drug** = test drugs in
  training cell lines. **Test B, unseen cell line** = training drugs in test cell lines. **Test C, both** = test drugs in
  test cell lines. Combination samples are kept out of every test set.

Every development step (antibody mapping, solver work, fitting) may read the training portion only. The three test
portions are read once, at the end, for all registered arms together.

## Phases

1. **Mapping and coverage** (no values read beyond training): CPPA antibodies to Reactome readouts, as
   `bench/analysis/062/readouts.tsv` does for korkut; drugs to targets, as for decryptM; cell lines to DepMap.
2. **Expression gating, no fitting.** A protein whose gene a cell line does not express (log2(TPM+1) below a threshold
   fixed in phase 2's own registration) holds no signal: its nodes are held at baseline like a drug. One rule, measured
   on the training portion and the curator guard.
3. **Learned weights, small first.** One elasticity per edge type (activation, inhibition, catalysis, regulation,
   depletion), then per-edge weights pulled toward Reactome's 1 by a penalty, fitted on the training portion (plus the
   FIT datasets). Needs the solver's sensitivity to each weight (implicit differentiation of the fixed point); the
   label-invariant sweep (specs/013) and the minimise mode (specs/003) are the starting points.
4. **Comparators, same splits:** reachability; specs/117's responsiveness (learned from the training drugs); the
   linear regression of specs/118 where it applies; CellBox if feasible.

## Endpoints (fixed now)

- **Primary:** per-sample Pearson correlation between predicted log fold (the model's readout UI on log scale) and
  CPPA delta-p over the mapped antibodies, averaged per drug, then over drugs (specs/118's measure); 95% bootstrap over
  test drugs (Test A, C) or test cell lines (Test B), 2,000 resamples, seed 121.
- **Secondary:** balanced accuracy over down / not / up with a change threshold fixed in phase 1 from training data.
- **The frozen model is scored on the test portions too**, with the other arms, so the narrow claim (right direction for
  unseen drugs) gets its untouched test whatever phases 2-3 find.

Each phase registers its own hypotheses and gate before it runs, as a numbered amendment here.

**What it is not.** No Reactome Release 98 (still held back). No change to a default without Adam's decision.

## The split, as drawn (2026-10-28)

`split_121.py` at this commit, seed 121; output kept at `~/deltasignal-catalogs/analysis/121/split.tsv`, SHA-256
`7dec6c069db69bd4c3099fbcf2bf1f028cca8892a6824a0abf44676a6f2a1234`. Drugs: 79 training, 52 test, 23 seen (training
only), 10 stimuli (training only). Cell lines: 238 training, 79 test. Single-compound samples: training 3,740; Test A
(unseen drug) 873; Test B (unseen cell line) 2,029; Test C (both) 348. Two name checks were made before the split was
used: bafilomycin is LINCS P100's bafilomycin A1, so it is a seen drug; I-BET726 and I-BET762 are different compounds.
Anonymised CPPA compounds ("MEK inhibitor-01") are drugs like any other; their targets come from the class named.

## Amendment 1 (2026-10-28): phase 1 mapping, before any CPPA value is read

- **Antibodies** (`~/deltasignal-catalogs/analysis/121/readouts_cppa.tsv`, private; korkut map format, made from the
  549 antibody names only, every phospho coordinate checked in Release 97): 375 total-protein, 106 phospho (each with
  at least one curated site), 68 excluded (uncurated sites, cleaved forms, histone and ubiquitin marks, loading
  controls, names that do not resolve). 44 ambiguous names are marked "check" with the reading used. 470 of the 481
  usable readouts are in at least one top-level network (`20261027-1346_8b6527e_top114`).
- **Drugs** (`bench/analysis/cppa_targets_121.py`; specs/113's rule: kinobeads at 1 µM, then ChEMBL mechanism, then
  ChEMBL binding; anonymised compounds take their named class; a few CPPA misspellings corrected; ChEMBL's hyphenated
  code numbers tried): 114 drugs have a target Reactome knows: 58 training, 35 of the 52 test drugs, 21 seen. The other
  17 test drugs (cytotoxics and tool compounds without a signalling target) cannot be predicted and are not scored.
- **Accession to gene**: each human protein's primary gene name in Release 97 (11,469 accessions, none ambiguous), not
  the kinobeads subset and not synonyms (a synonym another protein shares would pin the wrong protein).
- **Predictions** (`bench/cppa_bench.py`, values-blind): the frozen model, one solve per drug per top-level network, as
  korkut_bench. Stimuli are not solved in this phase.

## Amendment 2 (2026-10-28): phase 2, expression gating, registered before it is built or run

**Cell lines to expression.** A CPPA cell line takes the DepMap 24Q2 model whose stripped or full name equals its own
(normalised as the split normalises names); otherwise the DepMap line whose name is the longest prefix of it, at least
4 characters (a derived line: resistant, engineered, labelled; e.g. A375 for A375BR); otherwise none. Measured before
registration from names only: 97 lines exact, 56 by parent, 164 none; 5,291 + 247 of the 6,990 single-drug samples have
expression. Samples without expression keep the frozen prediction and are reported separately.

**The rule.** In a cell line, a gene with log2(TPM + 1) < 1 (TPM below 1) is absent. Every network node that is, or
contains, a protein of an absent gene (resolved as the bench resolves a target gene, `pin_occurrences` at scope `all`)
is pinned at baseline (UI 1) for that cell line: a complex that cannot form carries no change, while the cell's own
baseline is the reference, so absence is not a knockout. The drug's target pins win where they overlap. One solve per
(drug, cell line, network); the readout and scoring are as in phase 1. No parameter is fitted; the threshold is the
conventional TPM 1 and is not tuned.

**Hypotheses (training portion only; development).**
- G1: on training samples with expression, the gated model's per-drug correlation (fold arm) beats the frozen model's
  on the same samples (bootstrap over drugs, interval for gated − frozen excludes 0).
- G2: the gated model's false-change share (training samples, change threshold fixed in phase 1) is lower than the
  frozen model's.
Whatever G1-G2 show, the gated model is one of the arms scored once at the end; G1-G2 decide only whether phase 3
builds on it. The curator benchmark has no cell line, so this rule cannot touch it.

## Amendment 3 (2026-10-28): from the training portion, before any test portion is read

Phase 1 scored the frozen model on the training portion (development; `score_121.py`; numbers kept with the private
analysis, `~/deltasignal-catalogs/analysis/121/research_private.md`). Two decisions follow, both registered before any
test value is read:
- **Change threshold** for the secondary balanced accuracy: |delta-p| ≥ log2(1.5), korkut's threshold (specs/062),
  for consistency rather than tuned.
- **A secondary endpoint by time.** On the training portion the model's correlation is concentrated in samples treated
  for 12 hours or less and falls in samples treated for 16 hours or more, as expected of a model of signalling without
  feedback or gene expression. The primary endpoint is unchanged (all samples). Added as a secondary: the same
  correlation on test samples with a parsed treatment time of 12 hours or less (`5min`, `30min`, `4hr`, `0.5hr`,
  `3d`, ... parsed as minutes, hours or days; an unparsed time is left out of this endpoint).

## Amendment 4 (2026-10-29): phase 2 result (training portion)

Expression gating was run for the 282 training (drug, cell line) pairs with DepMap expression (`cppa_gated_bench.py`,
code `ff81b6a`; every solve converged). **G1 fails**: on the same training samples, the gated model's per-drug
correlation is not higher than the frozen model's (the interval includes 0; the point estimate is slightly lower).
**G2 holds descriptively**: gating halves the model's change calls, and the share of them that are false falls. The
calls it removes are mostly on readouts that barely move, which the correlation hardly weighs. Numbers are in the
private analysis under the CPPA terms. As registered, the gated model remains one of the arms scored once at the end;
phase 3 does not build on gating. Adam, 2026-10-29: expression "really isn't a measure of protein level at different
parts of a pathway"; phase 3 will use CPPA's own baseline protein measurements (p0) as cell features instead.

## Amendment 5 (2026-11-05): the gated arm reads totals as cppa_bench does

Found by code review. specs/134 changed how `cppa_bench.py` reads a total-protein antibody (abundance copies, the
default since 2026-11-04), but `cppa_gated_bench.py` kept reading every form, so at the one test scoring the readout
difference would have been credited to gating. The gated bench now resolves totals exactly as `cppa_bench.py` and
prints `DS_TOTAL_READOUT` in its Protocol line. The phase-2 training run (amendment 4, code `ff81b6a`) read forms, and
its verdict stands as recorded. At the final scoring, the frozen and gated arms are run with the same
`DS_TOTAL_READOUT`, and the logs must show it.

## Amendment 6 (2026-11-05): the close, registered before any test value is read

Adam ("yes", to closing development and scoring the test portions once). Development on the training portion is over:
the readout rule is settled (specs/134, amendments 1–3 adopted), and the three cell-context rules each failed their
registered gate (phase 2 expression gating, amendment 4; the readout's resting level, specs/130; the target's resting
activity, specs/139). Nothing from specs/130 or specs/139 is carried.

**Arms, fixed now.**
- **Frozen model**: one fresh prediction run at HEAD (`scripts/run_cppa.sh`, values-blind) under the adopted readout
  (`DS_TOTAL_READOUT=abundance` = `abundance_a3`; the log prints the version, amendment 5). As a reproduction check,
  its case rows must match the R9 run with the three amendment-3 pathways spliced in (`analysis/134/a3_splice`) before
  any scoring; a mismatch stops the close.
- **Expression-gated model** (phase 2, the registered arm): `cppa_gated_bench.py` on portions A, B and C
  (`scripts/run_cppa_gated.sh`, values-blind), same readout. A sample without a gated prediction keeps the frozen one.
- **Reachability**, from the same case logs.
- Phase 4 named responsiveness (specs/117), the specs/118 regression and CellBox as possible comparators. None was ever
  built on CPPA training; building one now would be development after the close, so they are not scored.

**Scoring** (`bench/analysis/score_121_final.py`, committed with this amendment; the only code that passes
`final=True`): for each portion A (unseen drug), B (unseen cell line), C (both), the registered endpoints —
- primary: per-sample Pearson over mapped antibodies (at least 3), averaged per drug, then over drugs; 95% bootstrap,
  2,000 resamples, seed 121, resampling drugs (A, C) or cell lines (B);
- secondary: three-class balanced accuracy at |delta-p| >= log2 1.5 (amendment 3);
- secondary: the primary on samples with a parsed treatment time of 12 hours or less (amendment 3).

The test portions are read once, by that script, after this amendment and the code are merged. Every measured number
stays in the private notes (the CPPA terms); the spec records qualitative verdicts. CPPA then becomes FIT
(`bench/datasets.tsv`).

**Expectations, stated before the read** (no gate; the scoring happens whatever they show):
- E1: the frozen model's per-drug correlation on Test A is positive and near its training value.
- E2: the frozen model is above reachability on every portion (as in training).
- E3: the gated model is not above the frozen model (its training gate failed).

## The final scoring (2026-11-06): done, once

Run at `94e8f79` exactly as amendment 6 registered: the fresh frozen run reproduced the R9+a3 splice on every row
before scoring; the gated runs covered portions A, B and C with every solve converged; `score_121_final.py` read the
test portions once. Every measured number, and the full reading, is in the private notes (the CPPA terms); the repo
records only whether the registered expectations held:

- **E1 fails.**
- **E2 holds on Test B only.**
- **E3 holds** on the registered secondary; no paired correlation endpoint was registered for the gated arm, so
  nothing more is claimed for it.

Both CPPA portions are now FIT (`bench/datasets.tsv`); there is no remaining untouched RPPA test. Lun 2019 (CyTOF,
overexpression) remains the sealed reserve.
