# specs/116: a trained readout layer on the frozen model (korkut2015, decryptM)

**Status:** registered 2026-10-27, before any training was run or any trained prediction was scored.

Adam, 2026-10-27: no Reactome Release 98 and no new datasets yet; "I would like to see what you come up with with regards
to training". Earlier: train on public datasets, but show accuracy with tests that are not circular.

**Why this layer.** specs/115: 66-78% of errors are change calls for readouts that did not move; within 20 hops the
direction is right 76-87% where a change is called on both sides; beyond 20 hops the calls are worse than chance. The
frozen model's output (a predicted fold) and the readout's hop distance from the drug's target are the two
quantities that carry that pattern. A layer trained on them can learn **which** predicted changes to believe, without
changing the solver or the networks.

## Data (both FIT, `bench/dataset_roles.require_role(..., "fit")`)

The specs/114 case logs on the top-level build `20261027-1346_8b6527e_top114`, code `a47898d`, run `s114_top_kb`:
korkut2015 (`korkut_cases.tsv`) and decryptM with the kinobeads map (`decryptm/decryptm_cases.tsv`). Scored cases
whose measured call is not `conflict`. Hop distance: shortest path in the bench's adjacency from any pinned target node
to any unpinned readout node of the same case (`bench/analysis/korkut_hops.py` and a decryptM counterpart using the
bench's own `readout_uuids`); none if no route. No solve is re-run.

## The layer (fixed here, not tuned)

A lookup table, so every learned decision can be read:
- **Fold bin** of the frozen model's signed log2 predicted fold, `log2(pred_ui)`: exactly 0 (`pred_ui` = 0);
  ≤ −3.32; (−3.32, −1]; (−1, −0.585]; (−0.585, 0.585); [0.585, 1); [1, 3.32); ≥ 3.32.
- **Hop band:** 1-5, 6-10, 11-20, 21-40, > 40, no route.
- **Rule:** for each (fold bin, hop band) cell, predict the class c (down, not, up) with the largest
  n(c in cell) / n(c in training), i.e. the call that maximises balanced accuracy on the training cases. A cell with
  fewer than 10 training cases uses its hop band alone; a band with fewer than 10 uses NOT.
- **Topology-only layer** (control): the same rule over hop band alone, no DeltaSignal output.

## Evaluation (no case is ever predicted by a table that saw its drug)

1. **Leave one drug out, within each dataset.** decryptM: one fold per drug. korkut: one fold per single drug X;
   test = every drug set containing X, train = every set not containing X (as Korkut et al. 2015's leave-one-drug-out).
   A combination is tested once per drug it contains; all test predictions are pooled.
2. **Across datasets.** Train on all of korkut, predict all of decryptM; and the reverse.

**Metric:** balanced accuracy over down / not / up (primary; always-NOT = 0.333), macro-F1 reported too.
**Compared with:** the frozen model's own call, reachability, the topology-only layer, always-NOT.
**Uncertainty:** 95% bootstrap interval over held-out drugs (2,000 resamples of drugs, seed 116) for the difference
trained-layer minus each comparator.

**Hypotheses.**
- T1: the trained layer beats the frozen model (interval excludes 0) in leave-one-drug-out on both datasets.
- T2: the trained layer beats the topology-only layer on both datasets (DeltaSignal's output adds information beyond
  graph distance).
- T3: trained on one dataset, it beats the frozen model on the other, in both directions.
Each is reported whichever way it falls.

**What it is not.** No setting of the frozen model changes; this is a layer after it. A pass here is a development
result: both datasets are FIT, and the layer's real test is unseen data, which Adam has asked to hold back for now.
