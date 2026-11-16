# Feature Specification: Rank, don't call

**Feature Branch**: `specs/163-rank-not-call`

**Created**: 2026-11-16

**Status**: Draft, registered BEFORE any scoring

**Input**: Adam, 2026-11-16 ("ok to your plan"). Every attempt to turn the model's predicted magnitudes into better yes/no calls has failed (specs/137, 144, 161, 162). Yet the predicted |log2 fold| ranks real changes above non-changes consistently: AUROC 0.65–0.67 on korkut and 0.63 on CPPA (non-marker readouts). The question: is a per-condition RANKING of readouts by predicted effect a measurably useful product output where calls are not?

## User Scenarios & Testing *(mandatory)*

### User Story 1 - The registered measurement (Priority: P1)

Offline, on RECORDED predictions only (nothing solved, nothing fitted):
- **korkut**: the 98 record (`results/cff961d/r98fix_korkut/korkut_cases.tsv`); a condition is a drug set.
- **CPPA**: the specs/121 final frozen run (`builds/20261027-1346_8b6527e_top114/results/94e8f79/cppa_final_frozen/`, four shards); a condition is a single-compound sample (drug, cell line, time, dose). Numbers stay outside the repository (CPPA terms).

Per condition, each readout's prediction is the mean over the pathways that score it of the clipped log2 fold (as specs/121). The drug's OWN direct markers are excluded (specs/138 markers on korkut, specs/139 on CPPA). Readouts are ranked by predicted |log2 fold|, and tied readouts take the tied block's mean outcome. A readout is a HIT if it measurably changed in the PREDICTED direction (korkut's primary measured call; CPPA delta-p beyond ±log2 1.5).

**Primary metric**: precision@10, the hits among a condition's top 10 readouts, averaged per drug, then over drugs. **Comparators**: reachability (rank 1 for readouts a route reaches, else 0, ties at the block mean, so "any reachable readout") and random (the condition's hit base rate).

**Gate, per dataset**: the model's precision@10 exceeds reachability's with a 95% bootstrap over drugs (10k resamples, seed 163) excluding zero. Reported whatever they show: precision@20, recall@10, and the random baseline. A claim is made only for the dataset(s) where the gate passes.

### User Story 2 - Nothing else moves (Priority: P1)

Analysis only: no engine, call, or scoring change.

### User Story 3 - Serving, only per the gate (Priority: P2)

A pass on a dataset licenses a ranked view in the product ("readouts ranked by how likely they are to change"), with that dataset's measured precision@10 as its stated record. That would be a later, separate spec (engine field and WA view). A fail on both records the negative.

## Requirements *(mandatory)*

- **FR-001**: the scorer is committed and fixture-checked before the real read; role guards for korkut2015 and CPPA (FIT).
- **FR-002**: the metric, comparators, tie rule, exclusions, and gate are fixed as above.
- **FR-003**: CPPA numbers stay outside the repository (verdict only); korkut numbers go in research.md.

## Non-goals

- No fitting, no new cutoff, no new solves; no change to the served calls.
