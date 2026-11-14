# specs/155 research: where the pairs live and how the table serves

Resolved from the code and the recorded results, 2026-11-14.

## R1. The pair sources, axis by axis

**Decision**: two candidate axes, resolved separately.

- **korkut (drug-RPPA)**: the predicted side is the recorded canonical scoring
  (`korkut_cases.tsv` of the adopted specs/134 amendment-3 protocol on the
  frozen record — 10,469 scored cases with `pred_ui`, 78 drug sets; the sound
  Release-98 record reproduces these predictions exactly, specs/152). The
  measured side is NUMERIC and already in the dataset:
  `external/korkut2015/expr.csv` holds per-condition log2 ratios, and
  `korkut_bench.measured_calls` shows the exact join (node_Index.csv names ×
  expr_index.txt conditions → values per (drug_set, readout)). A drug set
  measured under several conditions contributes the per-condition values; the
  pair's measured fold is 2^log2 per condition, one pair per (drug set,
  readout, condition) with the drug set as the LODO unit. The fitting script
  re-reads expr.csv directly (with the specs/062 alignment guards) rather than
  re-running the bench.
- **decryptM (phospho)**: the recorded `measured` is the paper's curve
  REGULATION CALL, not a number (`decryptm_coverage.py`: "the measured
  direction is the curve's regulation call"). A numeric measured fold would be
  the curve's fitted effect size (plateau fold) from the curves archives. R1a
  (to resolve in implementation, before fitting): if the curve files carry a
  fitted plateau/effect parameter per (drug, site), decryptM contributes pairs
  with that number and the drug as the LODO unit; if extracting one number
  requires choosing among doses or re-fitting curves, decryptM is EXCLUDED
  from this spec with the reason recorded here — a dose-free predicted fold
  has no canonical dose-bearing measured partner, and inventing one is a new
  protocol, not a lookup.

**Rationale**: korkut alone already gives every class far more than the
30-pair minimum and a 78-unit LODO; decryptM is additive if and only if its
own record supplies the number.

**Alternatives considered**: re-running benches to emit measured numerics
(changes recorded artifacts for no reason — the raw datasets are on disk and
registered FIT); P100/HPN as axes (P100 was scored as reachability only, HPN
is AUROC-over-contexts — neither recorded a per-case predicted fold under the
served protocol).

## R2. The class key: reuse the specs/146 stratum machinery

**Decision**: the table is keyed by the SAME stratum the served confidence
uses — `readout_confidence` (src/api/readouts.jl) already derives the level
string (call, magnitude bin, railed, resolution, convergence; hierarchical
fallback with `min_cell`) and returns it as `stratum`. Conditioning (a) =
quantiles per (axis, stratum level). Conditioning (b) = additionally per
predicted-fold band inside the level, bands chosen as recorded-fold quantile
bins (monotone by construction; no new binning code paths server-side — the
bands ship in the table as explicit [lo, hi) fold ranges).

**Rationale**: one stratum definition product-wide; a second, slightly
different class key would make the two honesty statements (confidence and
range) disagree about what "answers like this" means.

## R3. The shipped table and the lookup

**Decision**: `src/api/fold_ranges.json`, the exact specs/146 shape:
`version`, `spec`, `fit_record` (the korkut cases file's results path and the
dataset identity), `catalog_build`, `solver_commit`, `min_cell` (30),
`nominal_shares` ([0.5, 0.9]), `gates` (the measured verdicts, including the
losing arm's numbers), `table`: per axis → per stratum level → either the (a)
entry or the (b) band list, each `{lo, hi, share, n}` per nominal share.
Loaded by a `load_fold_ranges_table()` twin of `load_confidence_table()`;
refused (with the reason served in `protocol.fold_ranges`) when
`catalog_build`/`solver_commit` do not match the serving record — the
specs/143/146 provenance rule verbatim. Lookup: after `readout_confidence`
returns `(confidence, stratum)`, a readout with a call looks up
(axis, stratum); a (b) table picks the band containing the predicted fold and
falls back to the level's (a) entry when the fold is outside every band
(edge-case rule from the spec: refinement, never extrapolation).

**Response shape**: each entry with at least one covered axis gains
`measured_fold_ranges`: a list of `{axis, share, lo, hi, n}` (one element per
axis × nominal share). The `protocol` object gains `fold_ranges` ("ok:
<version>" or the refusal reason) and `fold_ranges_min_cell`. Nothing else
changes; absence is field absence (FR-004).

## R4. The fitting script and the gates

**Decision**: `bench/analysis/155/fit_155.py`, pure python (no scipy; the
specs/153 Clopper-Pearson precedent for hand-rolled stats), calling
`dataset_roles.require_role("korkut2015", "fit")` (and decryptM's role if R1a
admits it). Steps: build pairs (R1); derive each pair's stratum with the SAME
rules the server uses (ported constants, pinned against each other by a test,
as `fit_confidence.py` does for specs/146); LODO by drug set: for each
held-out unit, intervals fitted on the rest, held-out pairs checked for
coverage and width; gate 1 per axis: the drug-bootstrap (10,000 resamples of
units) 95% interval of held-out coverage contains the nominal share, for each
share; gate 2: (b) passes gate 1 and reduces median held-out interval width
by ≥10% against (a). The shipped table is then fitted on ALL pairs of every
axis that passed, arm per gate 2's verdict. Verdicts and numbers land in this
research file and in the table's `gates` block whatever they are.

**Rationale**: LODO by whole drug is the specs/051 discipline inside FIT data;
the bootstrap-by-unit respects the only real independence structure korkut
has (conditions of one drug set are not independent).

## R5. Tests

**Decision**: new `test/test_fold_ranges.jl` (bare name registered in
`test/asserting_suites.txt` — NEVER a path, NEVER sort the file — and the
CLAUDE.md count table): the lookup on a fixture table (level hit, band hit,
band miss falling back to (a), min-cell refusal, untouched/call-less absence),
the provenance refusal, list shape, and FR-001 on the existing readouts
fixture (response identical when the table is absent vs present except the
new fields). `test_api_errors.jl` gains the malformed-table cases, the
specs/146 pattern.

## R6. Parity and the serving flow

**Decision**: FR-001 is proven the specs/154 way — the specs/141 readout dump
and specs/142 rails dump replayed through the feature branch with the table
installed, compared field-by-field on existing fields (the new fields
ignored by the existing comparators by construction, asserted not to change
any compared byte). The serving decision then follows the gates, is recorded
here, and the PR ships table + lookup together. WebsiteAngular consumption is
its own spec afterwards (the specs/154 → WA specs/026 pattern).

## The measurement (P1, run 2026-11-14 — the registered gates)

**R1a resolved first**: the decryptM record's own `all_sites.tsv` carries the
paper's fitted `effect_size` per (site, experiment) — plateau ratio − 1 (down
calls bounded near −1, up calls unbounded) — so the measured fold is
max(0, 1 + effect_size) and decryptM IS admitted (LODO unit: drug). Its cases
file records no convergence flag, so its strata resolve at the
converged-free levels. korkut pairs: 11,413 over 78 drug sets (canonical
`a490cdc/k134_abund_a3` scoring × `expr.csv` per-condition folds); decryptM:
5,931 over 18 drugs (`a47898d/s114_f92_hand` × effect sizes).

**A construction bug caught by its own coverage, before any gate reading**:
the first band construction sorted pairs by (predicted, measured) and split
tied predicted folds across bands — stratifying bands by the MEASURED value,
which is circular. It covered 17% at nominal 50% IN-SAMPLE, which is how it
was caught. Fixed: bands never split a tied predicted fold (`cell_b`), and
the rerun is the measurement below.

**Gate 1 (coverage, LODO by whole drug, 95% bootstrap over units) — PASSES
for both arms on both axes**:

| axis | arm | 50% share | 90% share |
|---|---|---|---|
| korkut_rppa (11,413 pairs, 78 units) | (a) | 0.499 [0.466, 0.533] | 0.898 [0.881, 0.914] |
| korkut_rppa | (b) | 0.499 [0.467, 0.531] | 0.895 [0.878, 0.911] |
| decryptm_phospho (5,931 pairs, 18 units) | (a) | 0.479 [0.417, 0.551] | 0.911 [0.886, 0.927] |
| decryptm_phospho | (b) | 0.472 [0.424, 0.545] | 0.909 [0.883, 0.922] |

**Gate 2 (arm (b) ≥10% median width reduction) — FAILS decisively**: 0.17%
on korkut, 3.0% on decryptM. Conditioning on the predicted fold narrows the
intervals by essentially nothing: the predicted magnitude carries no interval
information beyond the answer's class. This confirms specs/128–131's finding
(the magnitude signal is too weak to fit) from an entirely independent angle.

**Verdict, per the registration: serve arm (a) on both axes.** The shipped
table (`src/api/fold_ranges.json`, 9 KB) is fitted on all pairs, stamped for
the served record (`20261113-1311_baf608b_r98fix` / `12b5db3`; the pairs
carry the frozen-97 predicted folds, which specs/152 proved that record
reproduces exactly). What it says is the feature's point: a confident strong
down call (magnitude bin gt3, converged, not railed) on the drug-RPPA record
measured between **0.64× and 0.95×** half the time (n 95) — the model's
magnitudes overshoot measured reality, and now the product says so with the
record's own numbers.

## FR-001 parity, end to end (P3, 2026-11-14)

The whole recorded benchmark (every (pathway, gene, direction) of the tracked
98-record dumps, `specs/152-release-98-redo/g4/`) replayed through the
feature branch's SERVED stack — parse → pins → solve → readouts — on the
serving catalog with the fold-ranges table installed and serving: **all
24,356 recorded folds reproduce to the dump's own rounding (max |diff|
5e-7), zero past 1e-5** (`report/parity_folds_155.tsv`). The live response
carries `measured_fold_ranges` and `fold_ranges: ok: 2026-11-14` beside
unchanged folds, calls, railed flags, confidences and switch fields.

Two traps hit and recorded on the way:
- A first replay mixed the Release-97 record's pin dumps with the 98 serving
  catalog: 5 folds differed, all RUNX1 over-expression — the 97 stable_ids
  resolve a different pin set on the re-curated 98 network. Not a code
  effect; the right dumps show zero. Lesson: a parity replay's query dumps
  must come from the SAME record as the catalog they drive.
- A second replay silently compared a stale output file: the runner had
  crashed on a missing input (untracked working copies of the g4 dumps from
  the specs/154 session no longer exist) and the `| tail` pipe ate the error
  while the old TSV sat where the result should be. Caught because the
  "corrected" numbers were impossibly identical. The tracked g4 dumps are the
  durable inputs; untracked copies are not a record.

## Amendment 1 (pre-serving, 2026-11-14): the fitter mirrors the served population

The adversarial review of the branch caught a population mismatch: the first
measurement's pairs included cases with predicted fold exactly 1.0 —
untouched readouts, which the server NEVER serves a range (the specs/146
amendment-1 rule). They were most of the korkut pairs (11,413 → 2,177 touched;
decryptM 5,931 → 2,761). The fitter now excludes them from fitting and
evaluation both, so the gates are measured on exactly the population the
product serves. The first measurement stands above as the mixed-population
record; THIS one is the serving decision:

| axis | arm | 50% share | 90% share | median width cut |
|---|---|---|---|---|
| korkut_rppa (2,177 pairs, 78 units) | (a) | 0.491 [0.439, 0.540] | 0.892 [0.865, 0.915] | — |
| korkut_rppa | (b) | 0.490 [0.447, 0.531] | **0.875 [0.851, 0.897] — fails** | 25.9% |
| decryptm_phospho (2,761 pairs, 18 units) | (a) | 0.455 [0.395, 0.538] | 0.924 [0.872, 0.939] | — |
| decryptm_phospho | (b) | 0.439 [0.400, 0.514] | 0.917 [0.851, 0.936] | 26.4% |

**Verdict, per the registered gates: korkut_rppa serves (a); decryptm_phospho
serves (b)** — on the touched population the predicted magnitude DOES carry
interval information on the phospho axis (26% narrower, still calibrated at
both shares), while korkut's banded intervals under-cover at the 90% share and
are refused. The mixed-population "narrows nothing" conclusion above was an
artifact of drowning the touched pairs in untouched ones. The shipped table
(`src/api/fold_ranges.json`) carries these verdicts in its `gates` block.
