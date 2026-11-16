# specs/161 research: the measurement (2026-11-16)

The engagement-pin arm (one replay of every recorded korkut solve on the
served 98 build; env probe passed, every registry rule canonical) is compared
with the specs/157 canonical full-knockout run on the same rows. 9,476 scored
rows; 993 marker readouts are excluded from BOTH arms (the markers set the pin
values). Scorer `bench/analysis/161/score_161.py` (fixture-checked, committed
before the arm finished).

| | accuracy | balanced | Δ balanced [95% boot, 78 drug sets] |
|---|---|---|---|
| full-knockout pins | 0.8580 | 0.4043 | — |
| engagement pins | 0.8842 | 0.4070 | +0.0028 [−0.0060, +0.0089] |

**Registered gate: FAIL** (the interval includes zero). RAF set aside: +0.0031.

## The anatomy: a large, nearly one-sided correction that balanced accuracy barely registers

264 calls move: **256 fixed, 7 broken** (sign test p ≈ 10⁻⁶⁴).
- Weak or no engagement sets: 152 moved, 151 fixed, 1 broken.
- Strongly engaged sets: 112 moved, 105 fixed, 6 broken.

Reachable-case scorecard (reachable means the full-knockout fold is not 1; markers excluded):

| | sensitivity | specificity | precision | AUROC |
|---|---|---|---|---|
| full-knockout | 0.430 | 0.697 | 0.113 | 0.650 |
| engagement | 0.405 | 0.785 | 0.141 | 0.672 |

Why balanced accuracy is unmoved: the fixes are false change calls turned into
"no change". That raises the "not" class's recall, which is already high
across the whole record because unreachable rows are almost all correct
"not"s. Meanwhile "down" recall dips (sensitivity −0.025). Balanced accuracy
averages the three class recalls, so a large one-sided cleanup of false
alarms registers as small. This is the specs/144 pattern again: removing
wrong change calls cannot create correct ones.

## A correction this measurement forced

The reachable-case scorecard given to Adam earlier the same day (sensitivity
0.846, specificity 0.133, precision 0.213 on 1,999 reachable cases) INCLUDED
each drug's own direct markers. Those are 263 of the 449 "real changes", caught
98.9% of the time, because pinning a target trivially moves its own marker.
Excluding them, the record's reachable-case sensitivity is **0.645**,
specificity 0.183, and precision 0.109 (1,307 cases, 186 real changes). The
scorecard must exclude a drug's own markers, and from here on it does.

## Serving

Per the registration, a fail records the negative. The engine already accepts
any pin value, so a user who knows their drug's measured engagement can pin
at it today. What the gate failure withholds is a documented, validated
protocol claim. Whether the 256-vs-7 raw-accuracy cleanup justifies
documenting it anyway is Adam's call.
