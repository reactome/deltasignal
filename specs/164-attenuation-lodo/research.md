# specs/164 research: the measurement (2026-11-16)

Twelve arms (the identity plus the eleven specs/129 curve tables) on the
served 98 build; every server's curve rule was verified positively from a
solve before its benches ran; all three benches finished in every arm. Scorer
`bench/analysis/164/score_164.py` (fixture-checked, committed before any arm
finished). CPPA numbers stay private under CPPA's terms (verdict only here);
korkut and decryptM are public FIT data.

## The verdict: does NOT generalize (2 of 3 held-out folds pass)

**Selection is perfectly consistent**: with each dataset held out in turn,
the other two choose the SAME grid point every time, `cf_K3_h1` (activating
slots read g(f) = 4f / (f + 3): unchanged at rest, near-linear for drops,
saturating rises at 4×).

| held out | moved: fixed / broken | sensitivity (identity → chosen) | specificity | precision | balanced | gate |
|---|---|---|---|---|---|---|
| korkut (9,476 units) | 122 / 15 (p ≈ 3e−22) | 0.645 → 0.624 | 0.183 → 0.282 | 0.109 → 0.118 | 0.4031 → 0.4038 | **PASS** |
| decryptM (2,232 units) | 83 / 11 (p ≈ 4e−15) | 0.890 → 0.838 | 0.167 → 0.301 | 0.262 → 0.285 | 0.3839 → 0.3909 | **FAIL**: guard (−0.052 vs the 0.05 limit) |
| CPPA (private) | fixed far more than broken | drop within the guard | up | up | up | **PASS** |

The registration requires all three. decryptM's paired test passes strongly
and its precision and balanced accuracy rise, but its sensitivity falls by
0.052, just past the registered 0.05 guard. **That fold fails, so attenuation
does not generalize as registered.** This is a fail by the margin it was
registered with, recorded without relabelling.

## The context that decides adoption anyway: the curator and genetic axes

specs/129 measured this exact setting (`cf_K3_h1`) on the frozen build:
curator held-out balanced accuracy **0.8681 → 0.8001** (103 fixed / 697
broken), tuning 0.8681 → 0.6840, and experimental 643 → 535 right (25 / 133).
The curve that the three DRUG datasets agree on badly damages the GENE
perturbation axes (curator expectations, MP-BioPath knockouts and
over-expressions).

## The reading

Drug responses in cancer cells at 24 h and genetic perturbations are best fit
by different signal strengths: drug datasets consistently prefer saturated,
weaker propagation, while gene-level axes need today's full-strength
propagation. A single global attenuation cannot serve both, which is why
every attempt to "fix magnitude" globally has traded one axis for another
(specs/128–131, and now here). If attenuation is pursued again, it should be
a perturbation-type-specific option (a "drug response" mode), registered and
tested as such. It is not a default.
