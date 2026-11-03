# specs/129 research

## §0 specs/128 re-read by balanced accuracy (post hoc, 2026-11-03, before this spec was registered)

specs/128 compared arms by cases right. The Lun test (specs/127) and this stage use three-class balanced accuracy, so
the specs/128 arms (solver `110bae4`, frozen build) were re-read that way. Experimental excludes nothing here (845 valid
cases); korkut is its 11,062 scored cases (always-NOT 0.333).

| arm | curator held-out | curator tuning | experimental | korkut |
|---|---|---|---|---|
| `ec_off` | 0.8681 | 0.8681 | 0.6885 | 0.4133 |
| `ec_inh` | 0.8681 | 0.8669 | 0.6925 | 0.4117 |
| `ec_sub90` | 0.8565 | 0.8454 | 0.7020 | 0.4129 |
| `ec_all90` | 0.8524 | 0.8101 | 0.6977 | 0.4132 |
| `ec_all75` | 0.7996 | 0.6801 | 0.6671 | 0.4194 |

The +353 korkut cases of `ec_all75` are +0.006 balanced accuracy: they came from calling fewer changes on a set that is
91% "no change". specs/128's reading 1 ("measured drug data reward attenuation") holds for raw accuracy only; by
balanced accuracy the measured sets favour at most mild contraction (experimental best at `ec_sub90`).

## Runs

Solver `d210ae3` (no `src/` change since `110bae4`), frozen build `20261026-1349_8b6527e`, arms `cf_off` and eleven
grid points (`bench/analysis/129/`), curator, experimental and korkut axes, four at a time. Scored with
`bench/analysis/score_129.py` (`~/deltasignal-catalogs/analysis/129/score_129.json`); side reports with specs/128's
scorer against `cf_off` (`analysis/129/cf_*.json`). Cases in the objective: experimental 796 (74 genes, RAF set aside),
korkut 11,062 (78 drug sets).

## Cross-validation

| | |
|---|---|
| chosen per fold | `cf_Kinf_h075` in four folds, `cf_off` in one |
| full-data choice | `cf_Kinf_h075`: g(f) = f^0.75 (K = inf, h = 0.75, e = 0.75) |
| **CV gain** (objective, CV − identity) | **−0.0004**, 95% interval [−0.028, 0.030] |
| full choice, in sample | +0.018 [−0.023, 0.058] |

**The curve does not earn its contraction.** Under the registered verdict the identity stays the candidate: the learned
model starts stage 5 with today's curves.

## The grid (balanced accuracy; objective = mean of experimental and korkut)

| arm | e | objective | experimental | (with RAF) | korkut | curator held-out | tuning |
|---|---|---|---|---|---|---|---|
| `cf_off` | 1.00 | 0.5491 | 0.6848 | 0.6885 | 0.4133 | 0.8681 | 0.8681 |
| `cf_K30_h1` | 0.97 | 0.5478 | 0.6839 | 0.6875 | 0.4118 | 0.8640 | 0.8592 |
| `cf_K9_h1` | 0.90 | 0.5517 | 0.6900 | 0.6947 | 0.4134 | 0.8527 | 0.8105 |
| `cf_K3_h1` | 0.75 | 0.5464 | 0.6724 | 0.6694 | 0.4204 | 0.8001 | 0.6840 |
| `cf_Kinf_h075` | 0.75 | **0.5670** | 0.7194 | 0.7250 | 0.4146 | 0.8461 | 0.7611 |
| `cf_K30_h075` | 0.70 | 0.5491 | 0.6798 | 0.6761 | 0.4184 | 0.8054 | 0.6808 |
| `cf_K9_h075` | 0.63 | 0.5154 | 0.6084 | 0.6074 | 0.4224 | 0.7462 | 0.6109 |
| `cf_K3_h075` | 0.52 | 0.4829 | 0.5417 | 0.5328 | 0.4240 | 0.6705 | 0.5173 |
| `cf_Kinf_h05` | 0.50 | 0.4949 | 0.5686 | 0.5608 | 0.4212 | 0.6980 | 0.5382 |
| `cf_K30_h05` | 0.42 | 0.4705 | 0.5185 | 0.5103 | 0.4226 | 0.6449 | 0.4960 |
| `cf_K9_h05` | 0.38 | 0.4570 | 0.4899 | 0.4838 | 0.4241 | 0.6087 | 0.4667 |
| `cf_K3_h05` | 0.32 | 0.4518 | 0.4796 | 0.4743 | 0.4240 | 0.5790 | 0.4605 |

Cases right against `cf_off` (fixed / broke, McNemar p); rails held-out curator and experimental; pathways with a net
curator change:

| arm | held-out (16,778) | tuning (4,415) | experimental (643 / 845) | korkut (7,117) | rails | pathways |
|---|---|---|---|---|---|---|
| `cf_K30_h1` | 16,755 (30 / 53, 0.02) | 4,380 | 642 (1 / 2) | 7,106 (15 / 26) | 3,009→2,071; 216→176 | 16 |
| `cf_K9_h1` | 16,657 (47 / 168) | 4,193 | 629 (7 / 21, p 0.013) | 7,165 (75 / 27) | →1,706; →137 | 32 |
| `cf_Kinf_h075` | 16,593 (50 / 235) | 4,006 | 612 (22 / 53, p 0.0005) | 7,337 (259 / 39) | →1,844; →111 | 37 |
| `cf_K3_h1` | 16,184 (103 / 697) | 3,697 | 535 (25 / 133) | 7,431 (337 / 23) | →1,647; →124 | 61 |
| `cf_K30_h075` | 16,223 (86 / 641) | 3,684 | 537 (27 / 133) | 7,440 (362 / 39) | →1,579; →97 | 57 |
| stronger (h 0.5, or K ≤ 9 at h 0.75) | 14,192–15,704 | 2,810–3,403 | 279–443 | 7,571–7,895 | →1,203–1,524; →66–96 | 70–80 |

Every contraction loses curator calls; korkut cases right rise monotonically with contraction (fewer change calls),
while korkut balanced accuracy moves only 0.413 → 0.424.

## Where the chosen curve's in-sample gain comes from (traced)

`cf_Kinf_h075` against `cf_off`: experimental 88 cases change class, **net −31 cases right**, yet balanced accuracy
rises 0.035, because 21 of the fixes are measured "no change" (a small class, weighted up) and the 64 losses are correct
changes weakened to NORMAL (29 DOWN, 23 UP and 12 direction-crossed cases). The gains concentrate in RBL1 and RBL2 in
Mitotic G1 (+17 of the positive net); the losses are spread (ERBB3 −9, KIT −7, WNT5A −5, MDM2 −4, WNT1 −4, CCND1 −4).
korkut: net +220, 199 of them DOWN → NOT on measured no-change, concentrated in drug sets containing SR. A gain this
concentrated in a few perturbations is what the grouped folds remove: the fold without RBL1/RBL2 chose `cf_off`.

## Predictions

| prediction | result | verdict |
|---|---|---|
| **P1.** `cf_off` reproduces `ec_off` case for case | byte-identical on all three case tables | **holds** |
| **P2.** the full-data choice has e between 0.75 and 0.97 | e = 0.75 (`cf_Kinf_h075`) | **holds** (at the edge) |
| **P3.** CV gain < 0.02, interval includes 0 | −0.0004 [−0.028, 0.030] | **holds** |
| **P4.** every non-identity point loses curator held-out balanced accuracy | all eleven below 0.8681 | **holds** |
| **P5.** no point converges less than `cf_off` | `cf_K3_h1`: curator 1,723 vs 1,724, experimental 243 vs 244; all others ≥ control | **fails** (one point, one solve per axis) |

## Reading

1. **One shared curve does not help on measured data out of sample.** The best in-sample setting's gain is a handful of
   perturbations; cross-validated it is zero. Stage 5 carries the identity curves unless a later stage earns otherwise.
2. **Balanced accuracy and cases right disagree** on the chosen curve's experimental effect (+0.035 against −31 cases).
   The objective was registered as balanced accuracy, matching the Lun scoring; both are reported.
3. **Contraction is not where the measured errors are.** It trades correct change calls for correct no-change calls at
   roughly the rate both classes are weighted, so balanced accuracy barely moves (korkut 0.413 → 0.424 across the whole
   range). What separates a responding readout from a non-responding one is not distance-uniform attenuation.
4. Per-class curves (3b) would split the same trade across classes; with a single curve's CV gain at zero they are
   unlikely to earn more and would cost three times the runs. Not run unless a reason appears. The next lever with a
   different mechanism is stage 2 (resting levels per cell: a gene that is off cannot respond).
