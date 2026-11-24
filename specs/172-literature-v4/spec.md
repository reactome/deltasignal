# Feature Specification: literature bundle v4 under curated pins (a specs/170 re-run)

**Feature Branch**: `specs/172-literature-v4`
**Created**: 2026-11-23
**Status**: Registered before any arm is run
**Input**: Adam, 2026-11-23: continue with what the literature source asked for, which is to score its two new
bundles under the specs/170 pin rule.

## Question

The literature source rebuilt its bundles from the same extracted facts as v3, in two forms:
- **v4-single:** every edge, attached once. It is the same 10,978 edges as v3-single, with attachment ties now broken
  by stable identity instead of uuid.
- **v4-strict-single:** only edges with a stated sign and a named mechanism (2,708). This selection has no root filter.

Under curated pins (specs/170), does either bundle improve the primary pool against the curated network?

## What is already known (checked before registering)

- **The curated network is the same as before.** Each bundle's curated part (`without`) is byte-identical to the
  specs/170 control (`20261120-lit-v3-without`) in every file the solver or the bench reads. The only differences are
  files neither reads (`single_copy_flagged.csv`, and extra literature rows in `leaf_refs.csv`, which only the
  off-by-default `DS_CARRIER_MULTISTEP` reads).
- **The code is the same.** `src/`, `bench/` and `scripts/` are unchanged from `be1fbe08`, the specs/170 run.
- **Pin check (the specs/170 FR-003).** Without the rule, literature edges unroot 266 curated roots in v4-single and 64
  in v4-strict-single. Under the rule there are 0 differences from curated root status in either.
- **The PIP3 edges are where they were.** The four literature edges specs/170 traced (PI3K ⊣ MTOR) are all in
  v4-single. v4-strict-single has none of them, and no literature edge into MTOR at all.

## Arms (run name `lit172`; every arm `--bench DS_SOLVE_TIMEOUT=3600 --bench DS_PIN_CURATED_ROOT=<builds>/20261123-lit-v4-single-without`)

| arm | build | role |
|---|---|---|
| control | `20261123-lit-v4-single-without` | base; must reproduce the specs/170 control (`v3-without`, `lit170`) class for class (FR-002) |
| gate | `20261123-lit-v4-strict-single-with` | the registered gate |
| described | `20261123-lit-v4-single-with` | described, not gated |
| shuffles | `20261123-lit-v4-strict-single-shuf01`–`19` (specs/166 degree- and sign-preserving rewiring) | only if the gate's net is > 0 |

## The registered gate (the specs/168 and specs/170 gate, unchanged)

The primary pool is korkut (with each drug set's own markers excluded) plus MP-BioPath experimental, over the cases
both arms score. Define net(X) = #(right in X, wrong in control) − #(wrong in X, right in control).

The gate PASSES if and only if both hold:
- net(gate arm) > 0;
- net(gate arm) is above every finished shuffle (rank p ≤ 0.05).

A shuffle still running after 21,600 s is stopped and recorded as not finished.

**Prediction, written before running.**
- **v4-strict-single fails the gate.** Removing the MTOR edges should leave experimental near 0 (PIP3 no longer
  collapses). korkut should still lose: on newly reached readouts the change-call precision was 0.015 under
  specs/170, and a sign-and-mechanism filter does not change what the edges reach. I expect a korkut loss somewhere
  between v3-strict's specs/168 figure and a quarter of v4-single's, because there are about a quarter as many edges.
- **v4-single reproduces specs/170 within a few cases** (gate net about −554), because it has the same facts and the
  same four edges. Changes in tie-breaking can move only a few attachments.

## Described, not gated

- For each pair: the specs/166 three numbers, change-call precision on newly reached readouts, and the curator
  held-out split required by the 2026-11-20 rule (curator results are reported, never gated).
- v4-single against `20261120-lit-v3-single-with` (`lit170`): paired, per axis. This measures what the tie-breaking
  change moved.
- The PIP3 over-expression cases in v4-strict-single: do the 78 cases that specs/170 traced read UP again?

## Requirements

- **FR-001:** this file, the scheduler and the scorer (`bench/analysis/172/`) are committed and pushed before any
  arm runs.
- **FR-002:** the control under `lit172` must give the same class on every case as `v3-without` under `lit170`. If it
  does not, the difference is traced before anything is scored.
- **FR-003:** nothing is fitted or selected. Results go in `research.md` whatever they say. The verdict goes to the
  literature source. This public repository gets verdicts only: no edge identifiers, bundle provenance or paper
  counts.

## Non-goals

- No change to the solver, the served catalog, any default or the protocol of record.
- No edge is chosen here. Which edges to keep is the literature source's decision, on its own data.
