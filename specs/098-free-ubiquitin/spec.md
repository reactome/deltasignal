# specs/098: free ubiquitin is a pool, not a signal

**Status:** pre-registered 2026-10-23, before the implementation and the arm. Arm `ubi` void (branch lacked specs/095). Arm `ubi2`: gate passes narrowly, held-out +6 (USP1), tuning +3, experimental −4 (CCNE1, whose route runs through free ubiquitin in the p27 degradation cycle). **Not adopted** (Adam, 2026-10-23: "lets not adopt then"); `DS_UBIQUITIN_MODE` stays `propagate`.

From the loop audit (specs/097), which looks for unhandled loop types by structure, not by
errors. Adam, 2026-10-23: "continue".

## The finding (specs/097, canonical `20261022-2011_233665c`)

Of the 57 unhandled positive (gain-1) cyclic components in the 92 benchmark pathways, **13,
holding 1,007 nodes in 10 pathways, are closed by free ubiquitin**: removing catalyst and
regulator edges leaves a material cycle, and it runs through a free ubiquitin node (R-HSA-68524,
R-HSA-113595, R-HSA-9660007). Ubiquitination consumes free ubiquitin; deubiquitination and the
proteasome release it; so every ubiquitinated substrate in a pathway is wired into one positive
loop through the free pool.

## Why hold it

specs/007 left ubiquitin out of the cofactor list because "ubiquitin and SUMO … transfer IS the
regulatory event". That is true of the **conjugate**: Ub-substrate is a different entity and
carries the signal. It is not true of the **free pool**. Free ubiquitin is abundant and not
limiting, and a perturbation of one substrate does not change the cell's free ubiquitin. Routing
signal through it couples every ubiquitination in a pathway, and its loops have gain 1.

## The rule (`DS_UBIQUITIN_MODE`, default `propagate`; `inert`)

Under `inert`, every node whose stable id is a free ubiquitin (the three stIds above, exactly; a
complex or conjugate containing ubiquitin is not touched) is held at baseline, as drugs
(specs/055) and pathogens (specs/048) are. It stays an AND input of every ubiquitination, so
completeness is unchanged. An observation that pins wins, as for drugs. The ubiquitin-pass
depletion edges (specs/089) are untouched. The solve reports `ubiquitin_rule` and
`ubiquitin_held`.

## Census (after the implementation, before the arm)

`loop_audit_hold.jl` (specs/097's audit with a held free ubiquitin's in-edges removed) on canonical:
- **342 free ubiquitin nodes in 40 pathways**;
- ubiquitin-closed positive components **13 → 0**; unhandled components 97 → 85; pathways with an
  unhandled positive loop 30 → 24;
- one large component was closed by ubiquitin and also by material steps: it stays cyclic as a
  material loop (material-closed nodes 249 → 1,136).

## Arm

| arm | build | solver | against |
|---|---|---|---|
| `ubi` | canonical `20261022-2011_233665c` | `--server DS_UBIQUITIN_MODE=inert` | `results/b61b365/dflt095` |

**Gate (as specs/077):** curator held-out net ≥ 0, experimental net ≥ −5, no pathway below −10.

## Predictions

1. Only the 10 pathways with a ubiquitin-closed loop move a class (and possibly others whose
   free ubiquitin is not in a loop, but they should not move).
2. Of the curator cases whose route crosses a ubiquitin-closed component, more are fixed than
   broken.
3. No prediction for the held-out or experimental nets.
