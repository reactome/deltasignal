# specs/157 plan

## Constitution check

- I/IV additive-only: serving (if any) adds fields beside specs/146's; the
  141/142 dumps replay unchanged — the parity replay is the proof.
- II measure-then-claim: the registry freezes before the run; the gate is
  registered (research R2/R4); numbers commit whatever they say.
- III fixtures compute: the scorer is fixture-checked before the run; served
  lookup tests follow the specs/146 suite pattern.
- V honest absence: non-contributing variants, excluded pairs, silent
  classes, provenance refusal — all reported, never substituted.

## Phases

1. **P1 registry**: registry.tsv + research R1 (DONE at registration; frozen
   on commit).
2. **P2 runner + scorer**: bench/analysis/157/ensemble_157.py (per-variant
   replay servers, positive diagnostic probes, recorded-case replay) and
   score_157.py (R2 labels, per-axis gate, seed 157), scorer fixture-checked.
3. **P3 the run**: up to 4 servers, variants sequenced; TSVs committed under
   the spec; research records the verdict, flip anatomy, informativeness.
4. **P4 serving (gate-dependent)**: the specs/146-pattern table
   (fit script role-guarded, provenance, min cell 30), /api/readouts fields,
   test_fragility.jl, docs, CLAUDE.md row, parity replay, PR.
   A full gate failure serves nothing without Adam's recorded call.
5. WA consumption only if served — its own spec.

## Structure

- specs/157-structural-fragility/{spec,research,plan}.md, registry.tsv,
  run TSVs + report under the spec.
- bench/analysis/157/{ensemble_157.py,score_157.py}.
- If served: src/api (table + readouts fields), test/test_fragility.jl.
