# Implementation Plan: Soft observations

**Branch**: `specs/156-soft-observations` | **Date**: 2026-11-14 | **Spec**: [spec.md](spec.md)

## Summary

A soft observation (own request field, never the existing confidence) makes
its node settle at `w*y + (1-w)*F(x)` under every solve path, default OFF;
a registered korkut leave-one-readout-out arm decides whether solver-level
conditioning helps before anything serves.

## Constitution Check

- IV (default off): DS_SOFT_OBS=off; the flag on with no soft entries is
  bit-identical by construction; dumps replay unchanged (FR-001).
- II (measure then claim): the arm (R4) runs before any serving decision;
  weights registered, no fitting; numbers recorded whatever they are.
- III: fixtures compute the blend's fixed points; mutations proven caught.
- V: overridden soft entries are reported with reasons, never averaged
  silently; flag-off requests carrying soft entries are refused, not ignored.

## Phases

1. **P1 solver**: the soft map threaded to the R1 sites; root rule; seeding;
   diagnostics. 2. **P2 API/CLI**: request field, 400s, response fields.
3. **P3 tests + parity**: test_soft_observations.jl, registration, 141/142
   replays with the flag ON (and no soft entries). 4. **P4 the arm**:
   bench/analysis/156/arm_156.py via run_arm.sh; research records the
   verdict; serving decision (gate pass → serve; fail → Adam's call on
   shipping the enabler unmeasured, US3 scenario 2). 5. WA consumption only
   if served.

## Structure

- src/solvers/steady_state.jl (R1 sites), src/api/server.jl, cli
- test/test_soft_observations.jl, asserting_suites.txt, CLAUDE.md, docs/API.md
- bench/analysis/156/arm_156.py, specs/156-soft-observations/ (record)
