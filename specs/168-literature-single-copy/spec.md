# Feature Specification: Literature edges attached once, and literature from new papers

**Feature Branch**: `specs/168-literature-single-copy`

**Created**: 2026-11-20

**Status**: RESULT: FAIL (2026-11-20, research.md). Registered before any arm was built or run.

**Input**: Adam, 2026-11-20: the FullTextMining session "left something in the mailbox for you to run". Their 014 is
a single-copy rewiring of their v3 bundle (`composite-h12b-v3-single`), answering our mailbox 014. That note found
the specs/166 bundle fed 431 AND slots with several copies of one curated entity (up to 20 in one slot), so a
perturbation counted once per copy. Their 015 is three pilot bundles from 441 papers the pipeline had never seen,
for PIP3/AKT, TP53 and Mitotic G1/S (`composite-pilot`, `-strict`, `-single`).

## Questions

1. **Gate.** Do literature edges attached once per end improve the primary pool against the curated network?
2. **Copy.** How much of the specs/166 harm was copy multiplication? Compare the same edges, attached to every
   copy (v3) and once (v3-single).
3. **Description.** What do new papers add, and what does the strict selection do?

## Arms

Every arm is built with `bench/analysis/168/build_arms_168.py` (the specs/166 machinery: `without` removes
literature rows and nodes; shuffles are degree- and sign-preserving). Every arm is run with `scripts/run_arm.sh`
under run name `lit168`, on curator, experimental and korkut, with `--bench DS_SOLVE_TIMEOUT=3600`. Shuffles run
on experimental and korkut only.

| build prefix | source bundle (external/, manifests verified 2026-11-20) | arms |
|---|---|---|
| `20261120-lit-v3` | literature-h12b-v3 | with, without |
| `20261120-lit-v3-single` | literature-h12b-v3-single | with; shuf01–19 only if the gate's net > 0 |
| `20261120-lit-v3-strict` | literature-h12b-v3-strict | with |
| `20261120-lit-pilot` | literature-pilot | with, without |
| `20261120-lit-pilot-single` | literature-pilot-single | with |
| `20261120-lit-pilot-strict` | literature-pilot-strict | with |

Each single and strict bundle carries the same curated rows as its parent. That is checked when the arms are built,
so each parent's `without` arm serves as the base for all its variants.

## The registered gate (question 1, the specs/166 gate)

On the primary pool (korkut with each drug set's own markers excluded, plus MP-BioPath experimental), over the
cases scored in both arms:
net(X) = #(right in X, wrong in v3-without) − #(wrong in X, right in v3-without).

PASS if and only if both hold:
- net(v3-single) > 0;
- net(v3-single) is above every finished shuffle (rank p ≤ 0.05).

Shuffles run only if net(v3-single) > 0, because they cannot rescue a net at or below 0. A shuffle still running
after 21,600 s is stopped and recorded as not finished.

## Described, not gated

- **Copy:** v3-single against v3, per axis: fixed, broken and net, plus change calls in each arm and how many are
  right.
- **v3 and v3-strict against v3-without.** Strict was chosen from the specs/166 slice results ("sign stated"), so
  this is in-sample.
- **pilot, pilot-strict and pilot-single against pilot-without**, on three pathways. The papers were searched using
  korkut drug-target/readout gene pairs. Korkut reach is therefore targeted, though no label was used. These arms
  are descriptive only.
- For every pair: the specs/166 three numbers (coverage gained, accuracy on newly covered readouts, fixed/broken
  on previously covered ones), plus change-call precision on newly covered readouts. The last is the touched
  over-call rate that the korkut route split asked to be reported beside any reach gain.
- Curator held-out (SELECT) is reported, never gated.

## Requirements

- **FR-001:** this file, the builder and the scorer (`bench/analysis/168/score_168.py`, tests in
  `bench/test_score_168.py`) are committed before any arm is built.
- **FR-002:** nothing is fitted or selected on these results inside this spec. Any edge set chosen from them is
  judged on new data.
- **FR-003:** results are recorded in `research.md` whatever they say. The verdict goes to the FullTextMining
  mailbox. Per-edge attribution is not sent (an open question for Adam).

## Non-goals

- No change to the served catalog, rules or defaults.
