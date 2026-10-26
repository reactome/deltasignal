# specs/110 results

Solver e00d64b (no solver change; `DS_DEDUP_ACTIVATORS` predates this spec), canonical build `20261024-2211_7f4b28a`.
Control `dflt102`; arm `dedup110` (`DS_DEDUP_ACTIVATORS=1`).

## Behaviour checks (G1): fail

| catalog | arm | converged | rails | monotone | relabel |
|---|---|---|---|---|---|
| benchmark (21 pathways) | defaults | 655 / 669 (97.9%) | 36,270 / 97,108 (37.4%) | 46,588 / 48,554 (95.95%) | 1.3e-10 |
| benchmark | dedup | 669 / 669 (100%) | 23,238 (23.9%) | 44,121 (90.87%) | 1.0e-13 |
| sample (16 pathways) | defaults | 436 / 436 | 9,561 / 29,560 (32.3%) | 14,735 / 14,780 (99.70%) | 7.8e-14 |
| sample | dedup | 436 / 436 | 8,572 (29.0%) | 14,619 (98.91%) | 2.2e-14 |

Convergence improves and rails fall on both catalogs: TGF-beta 151 → 42, R-HSA-5693606 4,094 → 0, DNA
Damage Bypass 960 → 88, insulin receptor 266 → 6, EGFR 343 → 7, the MyD88 cascades to 0. Dose monotonicity
falls on both, so **G1 fails** as registered. Most new reversals are small. Read after the fact, not as a
registered check: pairs whose reversal exceeds 10% fall (benchmark 3,776 → 1,475, sample 124 → 76), while those
between 1e-6 and 1% rise. Several are railed plateaus that were exactly at the cap and now sit just below it.
The largest losses are Class I MHC (1,914 pairs) and RAF/MAP (R-HSA-5673001, 807; e.g. 2× reads 94.6×, 10×
reads 81.5×).

## Benchmark (G2–G5): fail

| axis | baseline | arm | net | fixed / broke | p |
|---|---|---|---|---|---|
| curator held-out | 0.9079 | 0.9046 | **−60** | 3 / 63 | < 1e-4 |
| curator tuning | 0.8765 | 0.8714 | −26 | 2 / 28 | < 1e-4 |
| curator tuning, RAF set aside | 0.8779 | 0.8726 | −26 | 2 / 28 | < 1e-4 |
| experimental (tuning), RAF set aside | 0.7538 | 0.7550 | +1 | 1 / 0 | 1 |

G3/G4: of the 63 held-out losses, 2 were rail-dependent at baseline (`rail_audit.py` against `mild108`), so at
least 61 are unjustified. Every loss is a change call that becomes NORMAL (43 DOWN, 20 UP). Losses by pathway:
Interferon gamma 34 (all SOCS1, both directions), Base Excision Repair 13 (MBD4, SMUG1, OGG1 knockouts), FGFR2
6 (FGF10), EGFR 5, PTK6 3, DSB repair 2. Tuning: Mitotic G1 25 (SRC, ABL1, JAK2, PTK6 knockouts).

The values show why: a moderate effect that the squared step carried past a threshold now falls inside NORMAL
(SOCS1 knockout 1.36× → 1.08×; FGF10 knockout 0.78× → 0.88×; SRC knockout 0.84× → 0.93×). For the BER
glycosylases the change is larger (MBD4 knockout 0.001× → 0.89×).

## Predictions

- D1 fails: rails fall on both catalogs, but G1 fails on dose monotonicity.
- D2 fails as a range (−60 is outside −20 to +5); its second half holds: the gate fails on G4, and the losses are
  not rail-dependent.
- D3 holds (experimental +1).

## Reading

specs/012's result replicates on today's networks, larger: counting a self-modifying complex once removes
rails, and it also removes a doubling that many correct curator calls lean on to clear the 0.85 / 1.15
thresholds. The losses are not artefacts by the specs/108 test. Not adopted; the double count stays. The loops
it closes therefore need a rule that acts on the loop, not on the step: a pool where the generator's
regeneration test drops the path (TGF-beta, insulin receptor, macroautophagy), a carrier for Ku-type cycles.
