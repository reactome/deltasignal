# Feature Specification: Literature edges with curated pins (the specs/168 re-run)

**Feature Branch**: `specs/170-literature-curated-pins`

**Created**: 2026-11-23

**Status**: REGISTERED 2026-11-23, before the rule was coded or any arm was run.

**Input**: Adam, 2026-11-23, confirmed in conversation: "Approved: run the stepping
stone (the new clamp rule only for genes a literature edge feeds) as a re-run of specs/168. Then register full option B
as its own spec." The background:
- Adam's rule: roots and terminals are computed on the combined network, curated plus literature.
- A perturbation acts on the gene's own copies whatever feeds them.
- In specs/168, 49 of the 90 broken experimental cases had the perturbed gene's starting form unrooted by a
  literature edge, so the protocol stopped pinning it. The traced case is PHF8 → ERBB2 under ERBB2 over-expression.

## Question

Under a pin rule where a literature edge cannot stop a curated gene being perturbed, do the literature edges, attached
once (v3-single), improve the primary pool against the curated network?

## The rule (bench `DS_PIN_CURATED_ROOT=<catalog dir>`)

Given the curated network, the `without` arm of the same bundle:
- a node that the curated network contains is a root exactly when it is a root there;
- every other node (a literature node) is a root when it is a root in the network as served, curated plus literature.
  PHF8, fed by nothing, is a root.

Implemented as one filter on the edge rows the bench reads when choosing pins. An edge into a curated node is dropped
unless the curated network has that edge (source, target, edge type, sign). The specs/024 fallback reads the same
filtered rows. Literature edges into literature nodes, and literature edges out of curated nodes, still count. The
solve itself uses every edge: the filter only decides what is pinned.

So for a curated gene, the pins are those of the curated network. A pin clamps the node and cuts its in-edges for that
solve, so under ERBB2 over-expression PHF8 → ERBB2 is cut, and under a PHF8 perturbation it carries.

This is the stepping stone, not option B. It keeps the curated pins rather than clamping every stock copy, so curated
networks are untouched and the specs/168 arms are a direct comparison. Option B is registered separately.

The one expected exception: a curated node that the curated network does not pin, because it has no out-edge (specs/024
and the 2026-10-15 week review), but that gains a literature out-edge, can now be pinned. The pin check (FR-003) counts
these cases.

## Arms

The catalogs are the specs/168 builds; nothing is rebuilt. Arms are run with `scripts/run_arm.sh`, run name `lit170`, at
one commit, on curator, experimental and korkut, with `--bench DS_SOLVE_TIMEOUT=3600` and
`--bench DS_PIN_CURATED_ROOT=<builds>/20261120-lit-v3-without`.

| arm | catalog | role |
|---|---|---|
| control | `20261120-lit-v3-without` | base. The filter is a no-op; it must reproduce `lit168` exactly (FR-002) |
| gate | `20261120-lit-v3-single-with` | the gate |
| v3 | `20261120-lit-v3-with` | described: every-copy attachment under the rule |
| shuffles | `20261120-lit-v3-single-shuf01`–`19`, built with `build_arms_168.py` | only if the gate's net > 0; gate axes only |

v3-strict and the pilots are not re-run. Strict has no literature edge into a curated root, so its pins cannot change.
The pilots did not finish in specs/168.

## The registered gate (the specs/168 gate, unchanged)

On the primary pool (korkut with each drug set's own markers excluded, plus MP-BioPath experimental), over the cases
scored in both arms:
net(X) = #(right in X, wrong in control) − #(wrong in X, right in control).

PASS if and only if both hold:
- net(gate arm) > 0;
- net(gate arm) is above every finished shuffle (rank p ≤ 0.05).

Shuffles run only if the net is > 0. A shuffle still running after 21,600 s is stopped and recorded as not finished.

**Prediction, written before running.** I expect the gate to fail on korkut. In specs/168, korkut's −305 was mostly
readouts that only the literature reaches: 401 newly covered, with change-call precision 0.009. The previously covered
korkut cases were 72 fixed / 51 broken. A pin rule does not change reach. I expect the experimental loss (−69, 72 of
90 breaks being PIP3 over-expressions) to shrink substantially, and ERBB2 over-expression in PIP3 to read UP again.

## Described, not gated

- **Protocol effect:** the gate arm against `20261120-lit-v3-single-with` under `lit168`, per axis (fixed, broken,
  net). This measures how much of the specs/168 loss was the pins.
- **v3 against the control**, per axis.
- For each pair: the specs/166 three numbers, change-call precision on newly covered readouts, and the curator held-out
  split required by the 2026-11-20 rule. Readouts the curated network already reaches are the harm check; readouts
  only the literature reaches are reported as uninformative. Curator results are reported, never gated.
- **The traced case:** ERBB2 over-expression in PIP3 (R-HSA-1257604). Its pins under each arm, and its p-AKT readouts.

## Requirements

- **FR-001:** this file, the rule (with unit tests that run real code on a fixture network) and the scorer
  (`bench/analysis/170/score_170.py`, reusing `score_168`) are committed before any arm is run.
- **FR-002, control parity:** the control under `lit170` must give the same class for every case as `v3-without` under
  `lit168`. If it does not, the difference is traced before anything is scored.
- **FR-003, pin check:** before running, every gene in every pathway of the gate and v3 arms is compared. Its pins on
  curated nodes under the rule are compared with its pins in the control. Every difference is listed with its reason.
  An unexplained difference stops the run.
- **FR-004:** nothing is fitted or selected here. Results are recorded in `research.md` whatever they say, and the
  verdict goes to the literature source.

## Non-goals

- No change to the served catalog, the solver, any default, or the protocol of record. The rule is a bench flag for
  this measurement, never a served protocol. The structural rules that also key on starting copies are
  untouched here; option B takes them up.
