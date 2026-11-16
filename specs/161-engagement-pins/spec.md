# Feature Specification: Engagement-scaled drug pins

**Feature Branch**: `specs/161-engagement-pins`

**Created**: 2026-11-16

**Status**: Draft, registered BEFORE any solve

**Input**: Adam, 2026-11-16 ("sure try it out"), after the false-alarm decomposition: of the 1,344 false change calls on reachable korkut cases, 703 (52%) are in drug sets where no drug was shown to engage its target in these cells (specs/138), and 561 (42%) are in mixed sets. Only 80 (6%) are where every drug engaged.

## The problem

The recorded protocol models "drug added" as "target gone": it pins the target's starting forms to 0. In this melanoma line, five of the twelve korkut drugs never measurably moved their own direct marker (specs/138), so the model propagates a confident effect that did not happen.

## The registered change (the ONLY difference from the recorded protocol)

Each drug's pins are set to its MEASURED engagement e (UI scale, 1 = baseline) instead of 0. e is computed by `bench/analysis/161/engagement_161.py` from the drug's own SINGLE-drug conditions (both doses) and its specs/138 direct markers, which were fixed before specs/138 read most drugs' values:

- inhibitor markers: e = min(1, 2^mean log2 marker fold);
- positive-sign markers (Nutlin: MDM2 blocked, so p53/p21 rise): e = min(1, 1/fold);
- no marker on the array (HDAC, PKC): e = 0, the recorded full knockout;
- a pin shared by two drugs' targets takes the smaller e.

The values are committed as `engagement.tsv` (901 0.51, AK 0.20, NT 0.89, P6 0.88, PLX 0.60, RY 0.68, SR 0.83, ST 0.90, Tm 0.28, ZS 0.39; HN/RO 0). Everything else is the recorded protocol: the same pin resolution, doses, solve, and /api/readouts valuation.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - The registered measurement (Priority: P1)

Solve every recorded korkut (pathway, drug set) with engagement-scaled pins on the served 98 build (replay server; nothing else changed). Compare against the specs/157 canonical run (same build, same valuation, the full-knockout pins) on the SAME rows.

**No leakage**: a readout that is a direct marker of any drug in the drug set is EXCLUDED from scoring in both arms. The markers determined the pin values, so they must never grade them.

**Gate**: balanced accuracy beats the canonical run with a 95% bootstrap over drug sets (10k, seed 161) excluding zero. Reported whatever it shows: raw accuracy, the RAF-set-aside figure, the reachable-case scorecard (sensitivity, specificity, precision, AUROC on predicted |log2 fold|), and the per-engagement-group anatomy (fixed and broken calls).

**What this test is and is not**: it measures the PRODUCT protocol "the user tells the model how strongly the drug actually hit its target". On korkut that engagement comes from the record's own marker data, so this is not an independent prediction of engagement. korkut is FIT, and the claim is scoped accordingly.

### User Story 2 - Nothing else moves (Priority: P1)

No solver or engine change. The curator and experimental axes are untouched.

### User Story 3 - Serving, only per the gate (Priority: P2)

The engine already accepts any pin value (observations take 0–100), so a pass needs no engine work. A pass would add a documented protocol, "pin a drug's target at its measured engagement", to docs/API.md and the scope statement, plus the WA option of entering a drug's measured engagement; adoption of that protocol is Adam's call. A fail records the negative.

## Requirements *(mandatory)*

- **FR-001**: the engagement values and the arm runner are committed before the run; role guard korkut2015 FIT.
- **FR-002**: the arm changes only the pin values; the env is verified by the specs/157 positive probe (every registry rule canonical).
- **FR-003**: the gate, exclusions, and reports as US1; numbers committed whatever they show.

## Success Criteria *(mandatory)*

- **SC-001**: the measurement is on the record with intervals and the engagement-group anatomy.
- **SC-002**: any served claim is scoped to "engagement supplied by the user".

## Non-goals

- No engagement PREDICTION (that would need cell-context data the model lacks; specs/139 failed on CPPA).
- No change to the recorded protocol's defaults, the solver, or the cutoffs.
