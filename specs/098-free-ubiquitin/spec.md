# specs/098: free ubiquitin is a pool, not a signal

**Status:** pre-registered 2026-10-23, before the implementation and the arm.

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

Free ubiquitin nodes per pathway, and how many of specs/097's 13 ubiquitin-closed components
remain cyclic once they are held (re-run of `loop_audit.jl` with the held nodes' in-edges removed).

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
