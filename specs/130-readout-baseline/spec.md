# specs/130: a readout's resting level in the cell (specs/125 stage 2a)

**Status:** registered 2026-11-03, before any CPPA baseline value is read for this purpose. Stage 2 of the learnable
model (specs/125 §8), in its smallest form. Data: CPPA v1.0 (specs/121; terms accepted by Adam for academic use), the
training portion only. Every CPPA-derived number stays in `~/deltasignal-catalogs/analysis/130/` and the private notes;
this file records the rule, the verdict and nothing measured.

## Why this form

specs/125 stage 2 proposed per-cell resting levels, with expression as the source. That was already measured: specs/121
phase 2 gated nodes by DepMap expression on the CPPA training portion and its gate failed (G1, amendment 4). Adam then
pointed at CPPA's own baseline protein measurements (p0) instead: expression "really isn't a measure of protein level".
CPPA measures p0 only for its ~480 antibodies, not for every node, so the first use of it is at the readout: **a readout
that rests near its background in a cell cannot fall much there, and a rise is diluted by the background**. This is a
measurement model on top of the frozen prediction, not a change to the network solve: it is the readout's resting
level, the one place p0 reaches.

## The rule

For a sample s (set σ, cell line c, drug d) and antibody a with frozen predicted fold F (specs/121's `fold` arm: the
predicted log2 fold, clipped at ±log2 100, mean over networks):

- **p0(a, c, σ)**: the mean of antibody a over the control samples of cell line c in set σ (compound named CONTROL,
  DMSO, PBS, MOCK or VEHICLE, any time point), from the p0_p1 RPPA tables (log2, centred per set).
- **q_θ(a)**: the θ-quantile of p0(a, ·, ·) over every training cell line and set that has controls.
- **Share above background:** w = 1 / (1 + 2^(q_θ(a) − p0(a, c, σ))).
- **Predicted measured change:** m = log2(w · 2^F + 1 − w). As w → 1 it is F (today); a knockout (F = −∞) reads
  log2(1 − w); a rise is compressed by the background.
- A sample whose (set, cell line) has no controls, or an antibody without p0 there, keeps F.

**Grid** θ ∈ {0.05, 0.25, 0.5, 0.75} and the identity (w = 1, the frozen model). No other parameter.

## Judge (fixed now)

The registered CPPA endpoint (specs/121): per-sample Pearson correlation between prediction and delta-p over mapped
antibodies, averaged per drug, then over drugs (`score_121.per_drug`), training portion only.

1. Training drugs are split into 5 folds (seeded, 130).
2. For each fold, θ is chosen on the other four folds' drugs (highest mean per-drug correlation; ties to the identity,
   then the smaller θ); it predicts the held-out fold's samples.
3. **CV gain** = mean per-drug correlation of the CV predictions − the identity's, over the same drugs; 95% interval
   from 2,000 bootstrap resamples of drugs (seed 130).
4. **Secondary:** three-class balanced accuracy at |delta-p| ≥ log2 1.5 (specs/121 amendment 3), same CV.

**Verdict.** The readout baseline earns its place if the CV gain's interval excludes 0 above. It is then an arm of the
one CPPA test scoring (stage 5, with the frozen model), with θ from all training drugs. Otherwise it is not carried
forward. It does not touch the curator benchmark (no cell line) or any default.

## Predictions (before any value is read for this spec)

- **P1.** The identity reproduces specs/121's frozen training score exactly.
- **P2.** At least 90% of training samples have controls for their (set, cell line).
- **P3.** The CV gain is positive, and smaller than 0.03.
- **P4.** The secondary balanced accuracy rises under the chosen θ (fewer false changes on readouts at background).

## Note on access

While checking the p0_p1 table layout for this registration, three control rows' first values were printed (set101,
cell line SKMEL5, a training line). No response (delta-p or treated p1) value was read.

**Lun 2019.** HEK293T is not in CPPA; this rule cannot reach the reserve test unless its own control cells serve as
p0 at the one scoring. That would be a separate registration before the seal is broken; it is not made here.
