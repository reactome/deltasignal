# AND multiplication fidelity, and the default AND mode

**Created**: 2026-09-17
**Status**: Landed
**Supersedes**: the `DS_AND_MODE=hill_log` choice in
`specs/002-upregulation-propagation` and its re-measurement in
`specs/009-solver-defaults`

## Why this exists

The specification for AND is multiplication of fold-changes capped at 100 —
Adam, 2026-09-17: *"making sure the math comes out when numbers are below 1 is
as important as those above one. 1/2 * 1/2 should be close to 1/4"* and
*"100*100 should be 100"*.

Measured against that specification, the shipped default was accurate above
baseline and badly wrong below it, and a knockout was not a knockout.

| inputs | specification | `hill_log` (was default) | error |
|---|---|---|---|
| 0.5 x 2.0 | 1.0 | 1.00005 | 0.0% |
| 2.0 x 2.0 | 4.0 | 3.9645 | −0.9% |
| 0.5 x 0.5 | 0.25 | 0.2523 | +0.9% |
| 0.25 x 0.25 | 0.0625 | 0.0670 | **+7.2%** |
| 0.1 x 0.1 | 0.01 | 0.0135 | **+35.2%** |
| 0.001 x 1.0 | 0.001 | 0.0027 | **+167.6%** |
| **0 x 100** | **0** | **0.0135** | never zero |
| 10 x 10 | 100 | 74.06 | −25.9% |

Two properties of that table matter more than the individual numbers.

**The error below baseline is systematically UPWARD**, so every
down-regulated value is lifted toward baseline. That is a directional bias
against detecting DOWN, not a rail artifact.

**`0 x anything` was never 0, and it ROSE with the co-input** (0.0007 at co=1,
0.0135 at co=100). An abundant partner could partially "rescue" a knockout,
which is the opposite of AND semantics.

## Root causes

**In `hill_log`: the sigmoid is symmetric in log space.**
`log_out = z_max * tanh(log_fold / z_max)` compresses downward exactly as hard
as upward. But there is nothing to saturate against downward. The internal
domain is `[0, 1]` with baseline `0.01`, so an upward fold is genuinely capped
at 100x — that ceiling is real and worth modelling — while a downward fold of
0.01 is perfectly representable. Compressing it is error, not a floor.

**In `hill_log`: a hardcoded `eps = 1e-6` in the fold ratio**, against a
baseline of `0.01`. A knockout contributed `log(1e-6/0.010001) = −9.21` instead
of `−∞`. Third instance of the epsilon sizing rule already applied to
`DS_HILL_SAT_EPS` and `DS_INHIBITOR_EPS`; this one was not even configurable.

**In `hill_sat`: the same hardcoded ratio epsilon**, flooring every fold at
`1e-4`, so `0 x 100` read `0.0100`.

**In `hill_sat`: `DS_HILL_SAT_EPS = 1e-5`**, the smooth-max floor against zero.
It bounds how far below baseline a value may travel and contributes nothing
above it. A fold of `0.001` read 20.7% high; a fold of `0.0001` read **452%**
high.

**In `hill_sat`: catastrophic cancellation in the 100x cap.** The smooth-min
was written `(raw + max − sqrt(d² + eps²)) / 2` with `d = raw − max`. Once
`raw` is large this loses `max` entirely: with 100 AND inputs at fold 2,
`raw = 1.3e28`, `sqrt(d²) == d` to machine precision, and the expression
evaluates `(raw + 1 − (raw − 1))/2` as **0**. The smooth-max then returned
`eps/2`, so a maximally **elevated** wide reaction read as ~0 — the saturation
inverted into a collapse. `Class_I_MHC_mediated_antigen_processing_presentation`,
which has reactions carrying hundreds of nodes, did not solve at all.

Fixed with the algebraically identical, numerically stable form
`max − eps² / (2·(d + sqrt(d² + eps²)))`.

## Why the previous comparison was void

`specs/002` and `CLAUDE.md` recorded **`hill_sat` scoring 90 cases worse than
`hill_log`**, with the note *"why compression helps is unexplained"*. That
comparison was run against a `hill_sat` carrying two defects unknown at the
time: it could not represent a knockout, and it inverted its own saturation on
wide reactions, losing an entire pathway to `solve_failed`.

Correctly implemented, the gap is **−16 held-out cases, not −90**. Roughly 74
of the 90 were bugs, not model behaviour.

## Method

Release97, one shared catalog build regenerated from LNG `main` for this
measurement. The pathway list has 93 entries; **81 are actually scored** — 11
tuning (5,100 cases) and 70 held out (18,808 cases). Of the rest, 8 have no
curator file and `R-HSA-9025112 / _NEW_ROCK_signaling_regulates_MRLC_phosphorylation`
does not exist in Release97 at all (the `_NEW_` prefix marks it as a
placeholder in MP-BioPath's list). Quote 81, not 93 — 93 is the input list. Guards applied before any number was taken,
because two earlier A/B runs in this project were invalidated by exactly these:

- **Catalog freshness gated automatically.** The catalog's own
  `cache/fingerprint.json` carries `src_sha256`; the harness recomputes it from
  the live LNG source and refuses to run on a mismatch. The pre-existing
  catalogs failed this gate (`d2dd9166…` vs `9d3a4391…`) and predated
  `containment.csv` entirely.
- **Per-arm env verified by `docker exec`** — `DS_AND_MODE`,
  `DS_HILL_SAT_EPS` and the mounted catalog's pathway count — and the arm
  refused if any disagreed with what was requested.
- **Pairing conditioned** on `n_gene_uuids`/`n_ko_uuids`, so both arms answer
  the same question. 0 of 23,908 cases dropped, all arms.
- **Tuning/held-out split** per the MP-BioPath protocol: the paper's ten are
  tuning, the other 83 are held out and are what is reported.
- **Denominators compared across arms before metrics.** This is how the
  cancellation bug was found: the broken `hill_sat` arm scored 23,788 cases
  against the others' 23,908, and a 0.003 macro-F1 gap on a smaller set is
  indistinguishable from a real result.

## Results

All arms: 23,908 scored cases, same catalog, conditioned, 0 dropped.

| arm | held-out net | held-out macro-F1 | fixed / broke | all-pathway macro-F1 |
|---|---|---|---|---|
| `hill_log` (was default) | — | **0.8152** | — | 0.7846 |
| `hill_sat` (specification) | **−16** | 0.8138 (−0.0014) | 11 / 30 | 0.7835 |
| `hill_log_asym` (hybrid) | −12 | 0.8143 (−0.0010) | 2 / 14 | 0.7840 |

`hill_log_asym` keeps `hill_log`'s upward curve and multiplies exactly below
baseline. It is −12, so **exact sub-baseline multiplication does not improve
classification accuracy on its own.** Whatever the downward compression is
doing empirically, it is not an artifact of the zero handling.

`DS_HILL_SAT_EPS` 1e-5 → 1e-9 was isolated as its own arm on the same catalog
and is **exactly neutral**: 0 of 23,908 predictions and 0 raw `pred_ui` values
differ.

## The experimental axis

Added 2026-09-17 after an adversarial review found the decision had been
justified on curator ground truth alone. `docs/RESULTS.md` calls the
experimental axis "the most important caveat in this document" — DeltaSignal
reproduces curator *reasoning* far better than its predecessor and is no
better at predicting experimental *outcomes* — so a default change measured
only against curators is a gap in the justification, not just in the
reporting.

Same fresh catalog, same guards, conditioned pairing. 10 pathways, 849 of
849 cases comparable:

| arm | macro-F1 | change-F1 | accuracy |
|---|---|---|---|
| `hill_log` | 0.6506 | 0.7919 | 617/849 = 72.67% |
| `hill_sat` | 0.6487 | 0.7891 | 615/849 = 72.44% |

**Net −2 cases, 0 fixed and 2 broke, both in
`Transcriptional_Regulation_by_TP53`. That is NOT a detectable difference:**
McNemar exact on the discordant pairs gives **p = 0.50**.

(An earlier draft of this table reported 11 fixed / 30 broke and p = 0.0043 for
the held-out split. Those were the POOLED discordant pairs across both
splits, mislabelled. The held-out figures are 9 and 25. The conclusion is
unchanged — the cost is real — but the numbers were wrong.)

The two axes are therefore not the same result, and an earlier draft of this
document wrongly described them as "marginally worse on both axes, in the same
direction, at comparable magnitude". They are not comparable:

| axis | fixed / broke | net | McNemar p | verdict |
|---|---|---|---|---|
| curator held-out | 9 / 25 | −16 of 18,808 | **0.0090** | small but REAL cost |
| curator tuning | 2 / 5 | −3 of 5,100 | 0.4531 | no detectable effect |
| experimental | 0 / 2 | −2 of 849 | **0.50** | no detectable effect |

As percentages (−0.085pp and −0.24pp) they look similar, which is exactly how
the over-read happened. On the paired test they are a real signal and a coin
flip.

So the honest statement is: correct AND arithmetic costs a small, measurable
amount of agreement with curator reasoning, and has **no measurable effect on
predicting experimental outcomes**. There is no empirical upside, but there is
no measured empirical cost either.

**Limitation, unchanged from `RESULTS.md`:** every pathway with experimental
evidence is inside the paper's tuning ten, so there is no held-out empirical
test available in this dataset — for either tool. The 849 cases cannot be
split tuning/held-out the way the curator set can.

## The decision, and what it costs

**Default changed to `DS_AND_MODE=hill_sat`.** This trades 16 held-out cases
(0.085pp, macro-F1 −0.0014) for correct arithmetic and for magnitude fidelity.

The magnitude argument is the decisive one and is the part that belongs in a
paper. A **single-input** relay has nothing to combine, so it should be the
identity. Fold at the end of a chain of single-input reactions:

| source fold | `hill_log` d1 | d5 | d10 | `hill_sat` d1 / d5 / d10 |
|---|---|---|---|---|
| 100x | 74.1 | 33.5 | **19.1** | 100 / 100 / 100 |
| 10x | 9.6 | 8.3 | 7.2 | 10 / 10 / 10 |
| 0.1x | 0.104 | 0.120 | 0.139 | 0.1 / 0.1 / 0.1 |

`hill_log` loses **81% of a 100x signal over ten hops, and the loss depends on
path length**. Median path length in this catalog is ~10 hops, so a genuine
100x perturbation arrives at a typical readout as 19x while the same
perturbation one hop away arrives as 74x.

That undercuts the reason `hill_log` was adopted. It was chosen in `specs/002`
for continuous outputs, to make rank-correlation and GSEA-style claims
possible. Magnitudes that shrink with depth are not a usable scale: two
readouts with identical biology and different path lengths receive different
predicted folds. `hill_sat` is depth-invariant.

So the choice is: 0.085pp of three-class accuracy, against arithmetic that
matches the specification and fold-magnitudes that can be compared across
readouts. Recorded explicitly so the −16 is not discovered later as a surprise.

## Open

**Why downward compression helps classification at all is still unexplained.**
It is now a 16-case effect rather than a 90-case one, but it is a real paired
signal in a consistent direction across two independent formulations
(`hill_sat` −16, `hill_log_asym` −12). The most likely explanation is that
lifting suppressed values toward baseline suppresses false DOWN calls in
over-coupled regions, which would make it a symptom of the over-coupling
problem rather than a property of AND. Untested.

**A residual floor remains below baseline.** Suppressed wide reactions floor at
the smooth-max epsilon rather than reaching the true product (`0.5^50 =
8.9e-16` reads ~5e-8 in fold terms). That is seven orders below the 0.85 DOWN
cutoff and cannot change a call, so it is recorded rather than fixed.

## Pinned by

- `test/test_and_curves.jl` — 71 assertions, including the specification table
  above, the zero cases, and wide reactions at widths 10 to 400. Nothing
  covered wide reactions before, which is why the cancellation survived.
- `test/test_config_validation.jl` — the defaults, including
  `hill_sat_eps == 1e-9`.
