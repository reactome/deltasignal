# DeltaSignal Constitution

DeltaSignal is the **engine + HTTP API** that solves logic networks produced by
the logic-network-generator. The UI lives in the WebsiteAngular workspace and
talks to this API over the contract in `docs/API.md`. These principles govern
what belongs here and what does not.

## Core Principles

### I. Processing is DeltaSignal's job (NON-NEGOTIABLE)

The generator represents pathways **as curators intended them to be designed**;
DeltaSignal figures out how best to process that representation.

The corollary is the one that bites: when a faithful upstream representation
regresses the benchmark, the fix belongs **here**, as a better aggregator or
solver capability — not upstream as an unfaithful encoding. The open worked
example is redundant alternatives: a set-valued catalyst marked OR is
biologically correct and costs 51 of 223 cases, because `max(and, or)`
discards OR-cluster loss entirely (fold stays 1.0 even at total kill). That is
a missing capability in this repo, not an upstream mistake.

### II. Measure, then claim

No behavioural change is described as an improvement without a benchmark
before/after and a paired significance test. Accuracy alone is not the metric —
optimise **macro-F1**, so that "predict no change" is not rewarded.

This principle exists because it has been violated. A "122 of 223 solves are
non-converged" finding was reported that turned out to be the measuring
change's own λ-scaling bug; the true figure is 0. Retract loudly and correct
the record when it happens.

### III. Negative results are results

Record the numbers for what did not work, in the spec, permanently. Most
things tried here have been neutral or negative, and the written record is
what stops the next attempt from re-running them blind. A spec that only
documents wins is a trap for the next reader.

### IV. New behaviour ships default-OFF

Anything that can change predictions arrives behind a `DS_*` flag, defaulting
to current behaviour, until it is measured. The **code defaults are the
validated winning configuration** — environment variables exist to override
for sweeps, not to carry the real config. A config that only lives in the
environment is a config that gets lost.

### V. Honest solver reporting

Convergence, residuals, and errors must reflect what actually happened.
Do not report a residual that was never computed, do not call a non-finite
result converged, and do not report success for a branch that was skipped.
A solver that lies about convergence is worse than one that fails loudly:
`last_change` was once only assigned in one branch, so the common case
returned `converged=true, residual=0.0` unconditionally.

### VI. The API boundary is a trust boundary

Validate shape, type, and range on every input at the boundary. Do not echo
user-supplied values back in errors; log them instead. Map malformed input to
4xx and reserve 5xx for genuine faults. Julia's coercions are permissive in
ways that silently corrupt scientific input — `Float64("80"[1])` is `56.0`.

## Additional Constraints

- **Scales are asymmetric and load-bearing.** Observations arrive on 0–100;
  `node_activities` are returned on 0–1. Benchmark clients depend on both.
  Baseline is `x₀ = 0.01` universally.
- **The API contract in `docs/API.md` is the interface.** Breaking it breaks a
  separately-deployed UI; change it deliberately and version it.
- Comparable benchmark numbers require a stated ground truth (curator vs
  experimental), pathway set, and thresholds. A number without those is not
  comparable to another number.

## Development Workflow

- Paired comparisons use McNemar's exact test; report p, not just the delta.
- Verify the container environment actually carries the flags before quoting a
  benchmark number — Docker env passthrough has silently dropped them before.
- Findings outside the change at hand become tracked tasks in `specs/`, not
  inline scope creep.

## Governance

This constitution supersedes convention and habit. Specs and plans under
`specs/` are checked against it; a plan that conflicts with a principle must
either change or state the justification explicitly. Amendments require a
stated rationale and a version bump.

**Version**: 1.0.0 | **Ratified**: 2026-09-09 | **Last Amended**: 2026-09-09
