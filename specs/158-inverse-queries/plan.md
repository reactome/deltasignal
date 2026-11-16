# specs/158 plan

## Constitution check
- I/IV: a new endpoint only; the 141/142 dumps replay unchanged.
- II: cost and ordering MEASURED before defaults (research R1/R2); the
  influence suggestion is recorded as a measured negative.
- III: fixtures compute; inverse-vs-forward bit equality is the central test.
- V: caps, ordering, unreachable and unmoved counts all stated; refusals
  reuse the readouts endpoint's honest errors.

## Phases
1. **P1 research** — DONE (R1 cap, R2 ordering, R3 shape).
2. **P2 endpoint**: /api/inverse in server.jl reusing the pins/solve/readouts
   internals verbatim; docs/API.md.
3. **P3 tests**: test_inverse_queries.jl registered; mutation-proven.
4. **P4**: parity replay (dumps unchanged), adversarial review, PR, CI,
   merge, live smoke. WA consumption is its own spec once served.
