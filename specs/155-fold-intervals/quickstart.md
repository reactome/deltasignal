# specs/155 quickstart: validating measured fold intervals

## 1. The fitting run (P1, before any serving claim)

```bash
python3 bench/analysis/155/fit_155.py \
  <canonical korkut results dir>/korkut_cases.tsv \
  ~/deltasignal-catalogs/external/korkut2015 \
  specs/155-fold-intervals/report/
```

Expected: the role guard admits korkut2015 (FIT); the R1a decision on
decryptM printed and recorded; per-axis LODO coverage for both arms and both
shares with the drug-bootstrap intervals; the gate verdicts; the table
written only for a passing arm.

## 2. Unit fixtures

```bash
docker compose -f docker-compose.dev.yml run --rm julia-api \
  julia --project=/app /app/test/test_fold_ranges.jl
```

Expected: lookup hit/miss/fallback, min-cell refusal, provenance refusal,
absence rules, and the existing readouts fixture unchanged field-for-field.

## 3. Parity replays (FR-001)

Replay the specs/141 readout dump and specs/142 rails dump through the
feature branch with the table installed (the specs/152 g4 pattern). Expected:
zero differences on existing fields, 24,356 folds.

## 4. The live scenario

Solve a benchmark perturbation on the served catalog; a readout with a call
in a covered class carries `measured_fold_ranges` matching the shipped table
exactly; an untouched readout carries nothing; `protocol.fold_ranges` names
the table version.
