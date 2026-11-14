# Implementation Plan: Measured fold intervals beside the point fold

**Branch**: `specs/155-fold-intervals` | **Date**: 2026-11-14 | **Spec**: [spec.md](spec.md)

## Summary

Serve, per readout with a call, the measured-fold interval(s) the record holds
for answers in its class — per axis, with pair counts — from a provenance-
stamped table fitted by a registered script under pre-registered out-of-sample
gates. Additive only; honest absence everywhere the record is silent.

## Technical Context

- **Engine**: Julia API; the lookup mirrors specs/146's confidence table
  (`src/api/confidence.json`, `load_confidence_table`, `readout_confidence`)
  — see research.md R2/R3.
- **Fitting**: python under `bench/analysis/155/`, pure stdlib, role-guarded
  (R4). Pairs: korkut recorded predictions × `external/korkut2015/expr.csv`
  numerics; decryptM admitted only if its curves carry a fitted effect size
  (R1a), else excluded with the reason recorded.
- **Registered constants** (constitution II, set before measuring): nominal
  shares 50% and 90%; min cell 30 pairs; gate 1 = drug-bootstrap 95% coverage
  interval contains the nominal share, per axis and share; gate 2 = (b) valid
  and ≥10% median held-out width reduction over (a); LODO unit = drug set
  (korkut) / drug (decryptM).
- **Unknowns**: R1a (decryptM numeric availability) — resolved by inspection
  before fitting, never by inventing a dose protocol.

## Constitution Check

- I/IV (additive, default-on honesty fields): FR-001 parity replay is the
  proof, the specs/142/146/154 precedent.
- II (measure then claim): gates registered in spec + research BEFORE the
  fitting script runs; all numbers recorded whatever they are.
- III (tests that can fail): fixture-driven suite, each new test proven to
  fail under a mutation before the PR.
- V (honest absence): FR-004; refusals served with reasons.

## Phases

1. **P1 — pairs + gates (bench)**: `fit_155.py` builds pairs, runs LODO +
   bootstrap, prints the gate verdicts, writes `fold_ranges.json` for the
   winning arm (or nothing); verdicts appended to research.md. R1a resolved
   first and recorded.
2. **P2 — serving (engine)**: `load_fold_ranges_table` + lookup in the
   readouts handler; `measured_fold_ranges` entries; `protocol.fold_ranges`;
   provenance refusal.
3. **P3 — tests + parity**: `test_fold_ranges.jl` (+ registration in
   asserting_suites.txt and CLAUDE.md), api-error cases, the 141/142 dump
   replays with the table installed.
4. **P4 — serving decision**: per the gates; research.md records the decision;
   PR with the adversarial self-review + green CI flow; WA consumption spec
   afterwards if served.

## Project Structure

- `bench/analysis/155/fit_155.py` (+ a small `report_155.py` if the verdict
  print outgrows the fitter)
- `src/api/fold_ranges.json` (shipped only if a gate passes)
- `src/api/readouts.jl` (lookup), `src/api/server.jl` (response fields)
- `test/test_fold_ranges.jl`, `test/asserting_suites.txt`, `CLAUDE.md` table
- `specs/155-fold-intervals/` (this record)
