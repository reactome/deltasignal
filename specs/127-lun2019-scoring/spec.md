# specs/127: how Lun 2019 will be scored

**Status:** registered 2026-11-02, before any value of the dataset has been read and before any model has been
learned (specs/125 stage 1 has not started). The data stay sealed (specs/126). Adam, 2026-11-02: "stick to best
practices. make sure to be as thorough as possible. no rushing to conclusions"; and, the same day, "make the best
decision you can along the way". Every choice below was made by me on that authority, from the paper's methods, the
design tables (S1, S2) and the file names; none from the data.

**What this spec fixes:** which cases are scored, how the measured call is made from the cell-level data, how a model's
prediction is made, the one primary measure, the comparators, the uncertainty, and when the seal is broken. It does
not fix which learned model is scored: that is whatever specs/125 freezes at its stage 5.

## When it is scored: once, at the end

Scoring the frozen model now would show us its errors on this dataset, and those would steer the learned model. So
nothing is scored until specs/125 stage 5. Before the seal is broken:

1. the frozen model's predictions (tag `freeze-2026-10-26`) and the learned model's predictions are generated
   values-blind (`bench/lun_bench.py --predict-only`, which never opens the archive) and committed;
2. the scoring code is committed and tested on synthetic data;
3. the commit hash is recorded in this spec's research.md.

Then the archive's mode is restored, the format check (below) runs, and the scoring runs once. If specs/125 never
reaches stage 5, when to score the frozen model alone is Adam's decision.

## Cases

- **Perturbations:** the basal screen: each of the 649 proteins of interest (Table S1) over-expressed in HEK293T, no
  EGF, 18 h. Not in the primary score: the EGF condition, the 1 h time course, the kinase–phosphatase pairs, the
  kinase-dead mutants and the A375 subset (secondary analyses below).
- **Readouts:** the 32 scorable markers of `bench/analysis/127/readouts.tsv` (Table S2, checked against UniProt and
  Release 97 on 2026-11-02). GAPDH, cleaved PARP, p-HH3 and GFP are excluded. A readout covering several genes (p-ERK1/2,
  p-AKT, p-p38, ...) is one readout; its model value is the mean over its sites' nodes, as specs/062.
- **Self pairs** (the protein of interest is one of the readout's genes, e.g. MAP2K1 against p-MEK1/2) are not scored:
  the antibody sees the over-expressed protein itself. They are the measurement's positive control (below).
- **Primary population P:** the (protein, readout) pairs for which the frozen model gives a scored prediction on the
  whole-Reactome build `20261027-1346_8b6527e_top114`: the protein has a root form that is pinned (protocol of record,
  `DS_PIN_SCOPE=root_cycle`, `DS_PIN_ROOT_IGNORE=depletion`), and a readout node that is not pinned is in the same
  network. P is fixed by the values-blind run and committed. Every model and comparator is scored on P; a model with no
  prediction for a pair in P is counted as predicting no change.
- **Also reported:** the whole panel (every non-self pair; unpredicted pairs count as no change), and the frozen
  92-pathway build.

## The measured call (fixed before any value is read)

The measured call is computed by us from the per-cell screen data (`complete_dat.csv`), not taken from the authors'
binned or pre-analysed files or tables.

1. **Scale.** Every marker value is transformed by arcsinh(x / 5) as provided. Calls depend only on medians and on a
   threshold calibrated on the same scale (step 5), so they do not depend on whether the file is already transformed.
2. **Transfected cells.** A cell is transfected if its GFP exceeds the 99th percentile of GFP among the untransfected
   control cells (of the same batch, if the file identifies batches; otherwise pooled).
3. **Dose bins.** The transfected cells of the FLAG-GFP control samples (basal) are pooled; their GFP quintiles define
   five dose bins. An over-expressed fusion is compared with FLAG-GFP at the same GFP level, which removes the effect of
   transfection and of GFP itself.
4. **Effect.** For a protein-of-interest sample s and readout r, in each bin b with at least 20 transfected cells of s:
   d_b = median of r in s's cells in b − median of r in the FLAG-GFP control cells in b (same batch if batches exist).
   E(s, r) = the median of d_b over those bins; at least 2 bins are required, otherwise (s, r) is not measured. A
   protein with several samples (replicates, plates) gets the median of its E. A sample with fewer than 100 transfected
   cells is left out.
5. **Threshold, from the controls only.** For each FLAG-GFP control sample, E is computed against the other control
   samples (leave one out), per readout. δ_r = the 95th percentile of |E| over control samples. With fewer than 10
   FLAG-GFP control samples, δ_r = the largest |E| among them.
6. **Call.** UP if E > δ_r, DOWN if E < −δ_r, otherwise NOT.
7. **Positive control.** On self pairs, at least half must read UP. If they do not, the measurement pipeline is wrong:
   scoring stops, and the failure is recorded before anything else is looked at. The fix is registered as an
   amendment; it is not tuned against model agreement.

## Format check (the only look before scoring)

Once the archive is unsealed, the first and only step before scoring reads: the column header of `complete_dat.csv`,
the distinct values of its sample-identifying columns (protein, condition, time, plate or batch, replicate, control
type), and `Data S1 description.docx` if the header is not self-explanatory. Reading stops at the first sign of results
in the description. No marker values are read. What is found is committed to research.md. If the file does not carry
what steps 1–6 need (for example, no way to tell FLAG-GFP from untransfected controls), the steps are amended and
committed **before** any marker value is read.

## The model's prediction

- **Over-expression:** the protein's gene pinned at its root forms at the bench's over-expression pin (80× at the
  default resting level), as the curator bench pins it.
- **Value:** per network holding the protein's root forms and a readout node, the mean over the readout's unpinned
  nodes, classified with the curator cutoffs (UP ≥ 1.15, DOWN < 0.85).
- **Several networks:** the pair's prediction is the call on which every network that predicts a change agrees; if
  networks disagree it is `conflict` and counted as NOT (and the number is reported); if none predicts a change, NOT.
- **Comparators.** Always-NOT. Reachability as specs/021, mirrored for over-expression: UP if any readout node is
  reachable from the pinned roots in a network holding both, otherwise NOT. The frozen model. Any learned model.

## The primary measure and the verdict

**Three-class balanced accuracy** (the mean, over measured UP, DOWN and NOT, of the share predicted correctly) on P.
Always-NOT scores 1/3 by construction. Uncertainty: 10,000 bootstrap resamples **of proteins of interest** (a protein's
readouts are not independent), 95% percentile intervals, for

- frozen model − reachability, frozen model − always-NOT;
- learned model − frozen model, learned model − reachability.

A difference is called a gain if its interval lies above 0, a loss if below, and unresolved otherwise. All four are
reported whatever they show. This is the test of record for specs/125: adoption remains Adam's decision.

## Secondary analyses (reported, never the verdict)

- The whole panel; the frozen 92-pathway build; a 90th-percentile threshold.
- Change versus no change: AUROC of |log fold| (models) against measured regulated or not; reachability as a binary
  score.
- Magnitude: Spearman correlation of model log fold with E, on pairs measured and predicted regulated.
- Per pathway and per protein family (kinase, phosphatase).
- EGF 10 min, if a time model exists by then (specs/125 stage 4); otherwise not scored.
- A375 with vemurafenib and/or CI1040, if the archive holds it (scored like a drug bench: BRAF, MAP2K1/2 inhibited).
- The kinase-dead mutants, as the model sees them: the model cannot represent a catalytically dead protein, so they are
  listed, not scored.

## What this test can and cannot say

It tests over-expression of a GFP fusion, transiently and at doses well above normal, in one non-cancer line. A tag can
change a protein's function; HEK293T carries SV40 large T. It cannot test knockouts or drugs (except the A375
subset). The resting levels of HEK293T are not used by the frozen model; a learned model that uses them takes them from
a source chosen before scoring (specs/126 records the candidates; any new download needs Adam's approval).
