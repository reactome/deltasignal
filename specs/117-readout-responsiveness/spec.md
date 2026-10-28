# specs/117: readout responsiveness learned from other drugs (training round 2)

**Status:** registered 2026-10-27, before any responsiveness feature was computed or any prediction scored.

Adam, 2026-10-27, on specs/116: "continue" (go-ahead for round 2). Still no Reactome Release 98 and no new datasets, so
cell-line expression (round 2, item 3) is not used here.

**Why.** specs/116 found that the frozen model's predicted fold and the hop distance cannot say which readouts move
(T1-T3 fail). The missing information is plausibly a property of the readout: some sites and antibodies respond to
many drugs, most to none. That can be learned from the *other* drugs without seeing the held-out one. The question that
matters for DeltaSignal is whether its output adds anything once responsiveness is known.

## Data

As specs/116: the specs/114 top-level case logs (build `20261027-1346_8b6527e_top114`, code `a47898d`, `s114_top_kb`),
korkut2015 and decryptM (kinobeads map), both FIT (`dataset_roles.require_role`). Scored cases, measured call not
`conflict`. A readout is the korkut antibody (`readout`) or the decryptM site (`uniprot`, `position`).

## Features, from training drugs only

For a readout r, over the distinct (drug group, r) measurements among the training cases (a measurement scored in
several pathways counts once):
- **Responsiveness bin:** no training measurement; changed share 0; (0, 0.1]; (0.1, 0.25]; (0.25, 0.5]; (0.5, 1].
- **Direction:** down if more training measurements are down than up, up if the reverse, else tie.
- **Frozen call:** the frozen model's predicted class for the case (down, not, up).

## Layers (lookup tables, rule fixed as in specs/116)

Each cell predicts the class with the largest n(class in cell) / n(class in training); MIN_CELL = 10.
- **R1, responsiveness only** (the control): cell (bin, direction); fallback NOT.
- **R2, responsiveness + DeltaSignal:** cell (bin, direction, frozen call); fallback the R1 cell, then NOT.
- **R3, responsiveness + reachability** (the second control): cell (bin, direction, reachability call); same fallbacks.

## Evaluation

Leave one drug out within each dataset, exactly as specs/116 (korkut: test every set containing X, train on every set
not containing X; predictions pooled). No transfer across datasets: korkut antibodies and decryptM sites are different
readouts, so responsiveness does not carry over. Balanced accuracy over down / not / up (primary), macro-F1 reported.
95% bootstrap over held-out drugs (2,000 resamples, seed 117). Also reported: frozen, reachability, always-NOT, and
specs/116's trained layer figures for reference.

**Hypotheses.**
- H1: R2 beats R1 on both datasets (interval for R2 − R1 excludes 0): DeltaSignal adds information beyond how readily a
  readout responds.
- H2: R2 beats R3 on both datasets: that information is more than whether a route exists.
- H3: R1 beats the frozen model on both datasets (responsiveness alone is a stronger predictor than the model).
Each is reported whichever way it falls.

**What it is not.** No setting of the frozen model changes. Both datasets are FIT; a pass is a development result whose
test is unseen data, which Adam has asked to hold back.
