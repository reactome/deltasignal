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
