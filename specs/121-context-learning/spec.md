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
