# specs/118 research

Run 2026-10-28: `bench/analysis/lr_118.py` (tested in `bench/test_lr_118.py`) on korkut2015 (CellBox commit 6bed6b94
files) and the specs/114 top-level korkut case log (build `20261027-1346_8b6527e_top114`, code `a47898d`,
`s114_top_kb`). Output `~/deltasignal-catalogs/analysis/118/lr_118.json`.

**Reproduction check passes.** Our LR, with D from `pert.csv`, gives a mean per-drug correlation of 0.770 against
Long et al.'s 0.784 (registered tolerance ±0.03). With D as the doses in `expr_index.txt` it is 0.772. `pert.csv` is
used throughout, as registered.

| arm | mean per-drug corr | drug-alone corr | balanced accuracy (82 proteins) |
|---|---|---|---|
| LR (published best, rebuilt) | 0.7697 | 0 (predicts no change) | 0.522 |
| LR + DeltaSignal's call | 0.7710 | **0.157** | 0.522 |
| LR + DeltaSignal's log2 fold | 0.7719 | 0.161 | 0.522 |
| LR + reachability's call | 0.7684 | −0.125 | 0.521 |
| LR + DeltaSignal, one c per response | 0.665 | 0.037 | 0.519 |
| LR + reachability, one c per response | 0.764 | 0.034 | 0.522 |

95% bootstrap intervals over the 12 held-out drugs:
- LR + DS − LR: [+0.0005, +0.0023]
- LR + DS-fold − LR: [+0.0002, +0.0049]
- LR + Reach − LR: [−0.0023, −0.0005]
- LR + DS − LR + Reach: [+0.0010, +0.0045]
- LR + DS-fold − LR + Reach: [+0.0008, +0.0069]
- LR + DS on the drug alone: [+0.089, +0.235]

**Hypotheses.**
- **H1 passes, but the gain is very small** (+0.0013 in mean correlation).
- **H2 passes.** Reachability's call makes LR slightly worse.
- **H3 passes.** For the drug given alone, LR has nothing to say. DeltaSignal's prediction correlates 0.16 with the
  measurement, positive for all 12 drugs. Reachability correlates −0.13: the scale it is fitted with on combinations
  points the wrong way for single drugs.

**Why the overall gain is small.** Each held-out drug has two single-drug conditions and 8-21 combinations.
In combinations, the partner drug's learned effect already explains most of the response. So the fitted scale c is
small (about 0.017 per unit call), and DeltaSignal's term barely moves the combination predictions. Fitting one c per
response overfits (0.665).

**Traced, the drug given alone:**
- **Right where the target's own signalling is curated.** AKT inhibitor: p-AKT S473/T308 and p-GSK3 down
  (correlation 0.45). Temsirolimus: p-S6 S235/S240 and p-p70S6K down (0.38). PI3K inhibitor: p-AKT and p-GSK3 down
  (0.18).
- **Wrong where feedback is missing.** Temsirolimus raises p-AKT S473 (measured +1.47); DeltaSignal says down. This
  is the known relief of mTORC1 → IRS1 negative feedback, which our networks do not carry. PLX4720: DeltaSignal
  calls p-ERK up and measured is down (the RAF rail, specs/054).
- **Near zero** for HDAC, MDM2 and STAT3 inhibitors, whose effects are transcriptional and not curated as signalling.

**Reading.** DeltaSignal adds real information about a drug the best published model has never seen (H3). It
barely improves that model's overall score, because the score is dominated by combinations whose partner drug is known.
The useful role is predicting new drugs, not refining combinations of known ones. korkut2015 is FIT; development
result.
