# specs/139: does the drug's target work in this cell? (a cell-context gate on CPPA)

**Status:** registered 2026-11-05, before any CPPA value is read for this purpose. Adam, 2026-11-05: "lets go with what
you recommend" (cell context on CPPA, after specs/138). Data: CPPA v1.0 (specs/121; terms accepted by Adam for
academic use), training portion only. Every CPPA-derived number stays in `~/deltasignal-catalogs/analysis/139/` and the
private notes; this file records the rule, the verdict and nothing measured.

## Why this form

Two sources of cell context have failed the CPPA gate: expression gating of every node (specs/121 phase 2) and the
readout's resting level as a damping of its predicted change (specs/130). specs/130's reading: damping trades correct
change calls for correct no-change calls at about even rates; what the model lacks is which samples respond. specs/138
found the largest single class of korkut errors is drugs that did not move their own target's marker in that cell line.
The frozen model gives every cell line the same prediction for a drug. This spec asks whether the cell's own resting
activity of the drug's target, measured before treatment, says whether the drug will do anything there: a kinase that
is barely active at rest leaves little for an inhibitor to remove.

## The markers (values-blind, from Reactome Release 97 and the drug targets)

- **A drug's primary targets** (from `~/deltasignal-catalogs/analysis/121/drug_targets.tsv`, specs/121 amendment 1):
  kinobeads targets with Kd_app within 10-fold of the drug's lowest; otherwise every ChEMBL mechanism target, or every
  gene of an anonymised compound's class. A drug whose targets come only from ChEMBL binding assays has no primary
  target.
- **A marker** of drug d is a CPPA phospho antibody (`readouts_cppa.tsv`) at least one of whose sites Reactome curates
  as made by a reaction (the phosphorylated protein in an output, possibly inside a complex or set to depth 3, and not
  in any input) whose catalyst is one of d's primary targets: the catalyst activity's `activeUnit` if it has one, else
  its physical entity, counted only if that unit is not a complex (a protein, or a set of proteins, to depth 3). A
  catalyst complex without an active unit does not count, because its kinase cannot be told from its partners (RAF:MEK
  complexes made MEK inhibitors "mark" RAF's own sites).
- Built by `bench/analysis/139/markers.py` (Neo4j queries bounded server-side by `apoc.cypher.runTimeboxed`). Measured
  before registration from names only: 34 of 58 training drugs and 12 of 21 seen drugs have at least one marker.
  Examples: MEK inhibitors mark p-ERK (T202/Y204); mTOR inhibitors 4E-BP1, S6K T389, PRAS40, ULK1 S757 and AKT S473;
  AKT inhibitors GSK3, PRAS40, TSC2, MDM2 S166 and p27 T157; EGFR inhibitors EGFR's own tyrosines. Lost to the complex
  rule: CDK4/6 inhibitors (cyclin D:CDK4/6 has no active unit), PI3K inhibitors (a lipid kinase).

## The rule

For a training sample s (set σ, cell line c, drug d):

- **p0(a, c, σ)**: the mean of antibody a over the control samples of c in σ (specs/130's `control_baselines`).
- **Eligible set**: σ has controls for at least 4 training cell lines. CPPA sets are normalised separately, so a resting
  level is only compared with other cell lines of the same set; with fewer than 4 lines a rank says little. (Measured
  from treatment summaries only: 36 of 62 sets have at least 4.)
- **Rank** r(a, c, σ) ∈ [0, 1]: the mid-rank of p0(a, c, σ) among the training cell lines' p0(a, ·, σ), scaled
  (rank − 1) / (n − 1).
- **Target activity at rest** A(s): the mean rank over d's markers with p0 in (c, σ).
- **Gate at θ**: if A(s) < θ, every prediction for s is 0 (log2 fold; no change). Otherwise the frozen prediction.
- A sample of a drug without markers, in an ineligible set, or without marker p0 keeps the frozen prediction.
- **Grid** θ ∈ {0.2, 0.4, 0.6} and the identity (no gate). No other parameter.
- **Scoring leaves out each drug's own markers** in every arm, the identity included (as specs/138 amendment 1): the
  gate reads the marker's resting level, and its delta-p shares that control, so scoring it would score the gate on
  regression to the mean.

## Judge (fixed now)

**Primary for this spec: three-class balanced accuracy** at |delta-p| ≥ log2 1.5 (specs/121 amendment 3), pooled over
(sample, antibody) pairs, as specs/130's secondary. The reason is structural, not measured: the registered CPPA
correlation is per sample and does not change when a prediction is scaled, and a constant prediction scores 0
(`score_121.pearson`). A gate can therefore move the correlation only by replacing a gated sample's correlation with 0;
the question this spec asks (does the sample respond at all) is a question of calls.

1. Training cell lines are split into 5 folds (seeded, 139); the unseen thing a cell-context rule must generalise to
   is a cell line.
2. For each fold, θ is chosen on the other folds' samples (highest balanced accuracy; ties to the identity, then the
   larger θ); it predicts the held-out fold's samples.
3. **CV gain** = balanced accuracy of the CV predictions − the identity's, on the same pairs; 95% interval from 2,000
   bootstrap resamples of drugs (seed 139), balanced accuracy recomputed per resample.
4. **Random-gate null.** For each drug, the same number of its eligible samples as the CV gate removed, chosen at
   random among its eligible samples; 200 draws (seed 139). The gate must beat the 95th percentile of the null's gain:
   otherwise its gain comes from calling NOT more often, not from knowing which cell.
5. **No harm on the registered endpoint:** the per-drug correlation (`score_121.per_drug` measure, markers left out),
   CV gate − identity, 95% bootstrap interval over drugs, must not lie wholly below 0.
6. **Comparator, reported (not a gate):** the same rule with A(s) the mean rank over every phospho antibody instead of
   the drug's markers (a cell's general phospho-signal level). If it gains as much as the marker gate, the result is
   read as cell responsiveness, not target activity.

**Verdict.** The gate earns its place if 3 (interval above 0), 4 and 5 all hold. It is then an arm of the one CPPA test
scoring (specs/121; Test B and C, unseen cell lines, are where it is meant to work), with θ chosen on all training
samples, test cell lines ranked against the training lines of their set from their own controls at that scoring.
Otherwise it is not carried forward. It does not touch the curator benchmark (no cell line), korkut (one cell line, so
no rank), Lun 2019 (one cell line) or any default.

## Predictions (before any value is read for this spec)

- **P1.** The identity's balanced accuracy, markers left out, is within 0.01 of specs/121's frozen training figure.
- **P2.** At least half of the training samples of drugs with markers are eligible (eligible set, marker p0 present).
- **P3 (premise, descriptive).** Among eligible samples, the share of non-marker antibodies that change
  (|delta-p| ≥ log2 1.5) is higher in samples with A(s) ≥ 0.5 than in samples with A(s) < 0.5.
- **P4.** The gate fails. specs/129–130 found every way of calling fewer changes traded correct calls at even rates, and
  most training samples are not of marker drugs. It is registered because it is the one cell-context form that decides
  whether a sample responds rather than how much; a pass would be the first cell-specific signal on CPPA.

## Note on access

Before registration only names were read: antibody names, drug names and targets, the split, and which cell lines have
control samples in which set (treatment summaries). No p0, p1 or delta-p value was read for this spec.
