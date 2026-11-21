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

## Amendment, 2026-11-20 (code review): the reachable scorecard counted untouched rows

`score_161.py` took "reachable" as `base_fold != 1.0` on replayed folds, which carry float noise, so
2,000-odd untouched korkut rows entered the scorecard. Re-scored with the record-precision rule
(`bench/analysis/fold_rules.py`), same runs, n = 1,241:

| | sensitivity | specificity | precision | AUROC |
|---|---|---|---|---|
| full-knockout | 0.659 | 0.166 | 0.113 | 0.640 |
| engagement | 0.621 | 0.408 | 0.141 | 0.733 |

The specificities above (0.697, 0.785) were mostly untouched rows read as "no change". The gate does
not use reach: it stands (+0.0028 [−0.0060, +0.0089], FAIL), as do the 256 fixed / 7 broken.

The anatomy grouped drug sets with a no-marker drug (HN, RO: engagement placeholder 0) as "strongly
engaged". Regrouped: weak/none 132 moved (131 fixed, 1 broken), strongly engaged 89 (82 / 6),
engagement unknown 43 (43 / 0). n = 1,241 here against 1,307 in the record: the replay valued korkut
phospho readouts differently (specs/157 amendment of the same date).

## Amendment, 2026-11-21 (issue #278): korkut phospho readouts valued as the record values them

Both arms above (the specs/157 canonical run and the engagement arm) sent korkut phospho readouts to
`/api/readouts` as `forms` queries: the first stable id that resolves, every copy including pinned ones.
`korkut_bench.py` reads the union of the readout's sites' ids over unpinned copies, and 442 of the 3,565
phospho folds differed from the record. The comparison above was internally consistent (same valuation
both sides) but not the record's statistic. `/api/readouts` now has a `sites` mode that reproduces the
record exactly (10,469 of 10,469 korkut folds and calls at print precision; specs/157 amendment of the
same date), and `arm_161.py` sends phospho readouts in it (`--phospho-mode forms` reproduces the arm
above).

Re-run because it is cheap (one replay, 1,469 solves, about 12 minutes): the engagement arm at `16ae4061`
(`runs/calls_engagement_sites.tsv.gz`) against the corrected specs/157 canonical run
(`runs_278/calls_canonical.tsv.gz` in specs/157), on the served 98 build, same rows (9,476 scored, 993
marker readouts excluded from both):

| | forms (2026-11-16, re-scored 11-20) | sites (2026-11-21) |
|---|---|---|
| full-knockout acc / balanced | 0.8580 / 0.4043 | 0.8545 / 0.4031 |
| engagement acc / balanced | 0.8842 / 0.4070 | 0.8832 / 0.4067 |
| gate Δ balanced [95% boot, 78 sets] | +0.0028 [−0.0060, +0.0089] FAIL | **+0.0036 [−0.0052, +0.0099] FAIL** |
| RAF set aside | +0.0031 | +0.0040 |
| calls moved: fixed / broken | 264: 256 / 7 | **288: 279 / 7** |
| weak/none · strong · unknown (fixed/broken) | 131/1 · 82/6 · 43/0 | 145/1 · 87/6 · 47/0 |
| reachable n | 1,241 | 1,307 (the record's own count) |
| full-knockout sens / spec / prec / AUROC | 0.659 / 0.166 / 0.113 / 0.640 | 0.645 / 0.183 / 0.109 / 0.667 |
| engagement sens / spec / prec / AUROC | 0.621 / 0.408 / 0.141 / 0.733 | 0.608 / 0.432 / 0.140 / 0.744 |

The verdict does not change: the registered gate fails. The one-sided cleanup is slightly larger (279
fixed / 7 broken). The reachable scorecard now counts the record's 1,307 reachable rows, and its
full-knockout row (sensitivity 0.645, specificity 0.183, precision 0.109) equals the record figures in
the "correction" section above, as it should. The 11-20 note's n = 1,241 discrepancy is explained: it
came from the `forms` valuation.
