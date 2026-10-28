# specs/120 research: drug inhibition depth fitted from the other drugs

Arms on `20261027-1346_8b6527e_top114`: `s114_top_kb` (depth 0, code `a47898d`) and `d120_010`, `d120_025`,
`d120_050` (code `9e0694f`; solver and bench code identical to `a47898d`). Each arm's `ARM.json` records its
`DS_PERTURB_UI_DOWN`, and `depth_120.py` checks it. Every arm converged 793 of 814 korkut solves, as at depth 0. Every
case is scored in all four arms (none dropped): korkut 33,488 leave-one-drug-out predictions over 12 drugs (17,836
cases, a combination counted once per drug in it), decryptM 5,096 over 18 drugs. Outputs:
`~/deltasignal-catalogs/analysis/120/{korkut,decryptm}.json`.

## Result

Balanced accuracy over down / not / up, leave one drug out (the fixed depths are all drugs, in-sample, for reference):

| | korkut | decryptM |
|---|---|---|
| frozen (depth 0) | 0.396 | 0.365 |
| **fitted depth** | **0.392** | **0.365** |
| fixed 0.1 / 0.25 / 0.5 | 0.397 / 0.396 / 0.397 | 0.364 / 0.363 / 0.370 |
| fitted − frozen, 95% interval | −0.005 [−0.010, +0.000] | −0.000 [−0.009, +0.007] |
| specs/117 R2 with the frozen call | 0.676 | 0.573 |
| R2 with the fitted call | 0.677 | 0.567 |
| R2 fitted − R2 frozen | +0.001 [−0.000, +0.004] | −0.006 [−0.010, −0.003] |
| reachability / always-NOT | 0.387 / 0.333 | 0.382 / 0.333 |

- **H1 fails** (korkut), **H2 fails** (decryptM), **H3 fails** (on decryptM the fitted call makes R2 slightly worse).
- Depths chosen: decryptM picks 0.5 for 17 of 18 held-out drugs (lapatinib 0). Korkut is split: 0 for five drugs,
  0.1 for three, 0.5 for four. The choice is not stable on korkut, and 0.5 is the edge of the registered grid.
- The false-change share hardly moves: korkut 0.869 → 0.862 at depth 0.5, decryptM 0.797 → 0.793.

## Trace: what depth 0.5 changes (against depth 0, all cases)

- **It removes false alarms, as intended.** On korkut, 483 change calls on unchanged readouts become NOT (266 UP, 217
  DOWN); on decryptM, 159. Folds roughly halve: the median |log2 fold| of the frozen model's change calls goes 1.77 →
  0.99 on korkut, 2.1 → 1.0 on decryptM.
- **But most false alarms stay.** A target pinned at 50% is still a large change against the 0.85 / 1.15 cutoffs, so
  86% of korkut change calls keep their call. The false alarms come from routes, not from how deep the inhibition is.
- **It also breaks right calls, and those weigh more.** Measured changes are rare (91% of korkut cases are NOT), so
  under balanced accuracy one lost right DOWN call costs about as much as fifteen removed false alarms. Korkut loses 23
  right calls to NOT and 33 right DOWN calls flip to UP.
- **Partial inhibition flips direction.** 125 korkut calls go DOWN → UP and 63 on decryptM, mostly in combinations with
  the AKT inhibitor 901 (MK2206) or PLX4720. A half-inhibited target lets a loop that the full knockout silenced respond
  the other way. This is the loop sensitivity specs/054 and 108 describe (rails), not a depth effect to fit.

**Reading.** One global depth is not a lever: the model's false change calls come from the routes it takes, and a
partial pin trades them against right calls almost one for one under balanced accuracy. Per-drug depths would need a
measured drug-specific inhibition (for example, the kinobeads potency at the dose used), which is new data.
