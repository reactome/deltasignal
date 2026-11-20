# Feature Specification: Is a change call better than the base rate among reached readouts?

**Feature Branch**: `specs/167-reached-base-rate`

**Created**: 2026-11-20

**Status**: Analysis, POST HOC. No rule, no gate, nothing served.

**Input**: a review of the record after specs/161–166 (Adam, 2026-11-20: "look through everything we have been doing and come up with the next thing"). specs/165 concluded that "the low precision is specific to drug readouts": about 88% of change calls are right on the gene axes, about 11% on korkut. That reading compares axes that differ in two ways at once. This spec measures how much of the difference is the model and how much is the axis.

## What was already seen

Before this file was written, an ad-hoc count on four axes of the served Release 98 record (curator held-out, curator tuning, experimental, korkut with own markers included) showed that the model calls a change on 89–97% of the readouts a route reaches, and that the share of those calls landing on a measured change is within about 0.02 of the share of ALL reached readouts that changed. Not yet seen when this was committed: korkut with own markers excluded (the scorecard convention), decryptM, any interval, and the readouts the model reaches but declines to call.

## The measurement

Read from recorded case files; no solve. Per axis:

- **reached**: the predicted fold is not 1.
- **base rate**: the share of labels that are a change, among reached readouts and among unreached ones.
- **call rate**: change calls / reached.
- **precision**: among change calls, the share whose label is a change (`change`) and the share in the called direction (`exact`, the specs/165 figure).
- **lift** = precision(change) − base rate among reached, with a 95% bootstrap over perturbation units (10,000 resamples, seed 167; units as specs/165: korkut drug sets, decryptM drugs, otherwise (pathway, gene, direction)).
- **declined**: reached readouts called no change, and the share of those whose label is a change.

**Axes** (the served Release 98 record, which equals the frozen record case for case, specs/152):

- curator held-out, curator tuning, experimental: `20261113-1311_baf608b_r98fix` / `12b5db3` / `r98fix_frozen`, valid cases;
- korkut, own markers excluded (and, beside it, included): same build / `cff961d` / `r98fix_korkut`, scored cases;
- decryptM: same build / `552ba89` / `r98fix_dm153`, scored cases.

CPPA is not read (its terms keep derived numbers out of this repository).

## How the result is read (fixed before the committed script ran)

- **An axis where the model's call adds to reachability** is one whose lift interval excludes 0.
- **If precision tracks the reached base rate on every axis** (call rate ≥ 0.85 and |lift| ≤ 0.05 everywhere), the difference in precision between the gene axes and the drug axes is a difference in how often a reached readout changed, not in how well the model picks among reached readouts. It then cannot be attributed to the perturbation type from this record, because the recorded axes confound two things: every gene-perturbation axis has readouts chosen because a change was expected or published, and every unselected readout panel is a drug axis.
- **Selection is visible in the unreached base rate**: on an unselected panel, readouts the network cannot reach rarely change; on a set assembled from reported effects they often do.

## Requirements

- **FR-001**: the reader (`bench/analysis/167/base_rate.py`) is fixture-tested (`bench/test_base_rate_167.py`) and committed before it is run on the record.
- **FR-002**: nothing is fitted. The curator held-out split (SELECT) is read descriptively, as in specs/165.
- **FR-003**: the result is recorded in `research.md` whatever it says, with the correction to specs/165's reading if one is due.
- **FR-004**: `research.md` states what data would separate perturbation type from readout selection, and how much of it the served networks could score (a structural count, no external data read).

## Non-goals

- No change to any call, cutoff, served table or scope statement; no new dataset is downloaded or scored here.
