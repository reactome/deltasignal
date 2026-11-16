# specs/158 research: cost, ordering, and the shape of the endpoint

Measured 2026-11-15 on the served 98 build (replay server, three pathways by
node-count percentile; 25 sampled knockout candidates each, every number an
ordinary forward solve + readout valuation).

## R1. Cost (sets the cap)

| pathway | candidates (distinct pin sets) | s/solve | full inverse |
|---|---|---|---|
| R-HSA-5685939 (small) | 22 | 1.63 (includes warm-up) | ~36 s |
| R-HSA-451927 (median) | 58 | 0.06 | ~3 s |
| R-HSA-73894 (large) | 499 | 2.07 | ~17 min |

A typical pathway's full inverse is seconds; the large tail is minutes.
**Registered: default cap 50 candidates, server maximum 500**; the endpoint
stays synchronous and the response always carries
`candidates_total / reachable / solved / moved / capped`, so a capped answer
can never read as exhaustive (FR-004).

## R2. Ordering — influence does NOT pre-rank (measured, negative)

The backlog note suggested influence_scores for ranking. Measured: take one
candidate solve's `influence_scores`, rank the other candidates by the max
|influence| over their pins, and ask whether top-half candidates move the
probe readout more often than bottom-half. On the only pathway with enough
movers to tell (R-HSA-73894, 16 of 25 sampled KOs move it): **top half 8,
bottom half 8** — no signal; the small pathways' single movers split 1/0 and
0/1. Influence explains the solve that produced it, not other solves.
**Recorded: no influence pre-ranking.**

The principled pre-filter is **ancestry**: a candidate whose pins have no
directed path to any of the target readout's resolved nodes over the solved
network's own edge set (derived edges included — twins, recycles, bridges all
live in `network.edges`) provably cannot move it through the network; it is
counted under `reachable`, never solved. Within the reachable set the order
is deterministic (ascending entity id) and the response SAYS the order is
arbitrary (`ordering: "entity_id"`) — no ranking claim is made before solving.
Returned rows are ranked by measured |log2 fold| AFTER solving, which claims
nothing the solves did not show.

## R3. The endpoint shape (resolved from the code)

`POST /api/inverse`: `{network_id, target: {stable_ids, members?, mode?},
direction: "down"|"up"|"any", dose?: {ko: 0, oe: 80}, cap?}`. The handler
reuses, with identical arguments, the exact internal paths the product
serves: pin resolution (the /api/pins machinery) for the candidate set, the
ordinary `solve_steady_state` call the solve handler makes, and the readouts
handler's valuation ladder for the target — FR-001's bit-equality test runs a
returned row's candidate through the public /api/solve + /api/readouts and
requires identical numbers. Rows carry the full honesty set (fold, call,
railed, confidence, measured ranges, fragility, switch state) exactly as the
readouts endpoint serves them; refusals reuse that endpoint's errors.

## R4. Tests

`test_inverse_queries.jl` (registered in the asserting list + CLAUDE.md): the
US1 fixture (A → B → C, D ⊣ C) with inverse-vs-forward bit equality, the
direction filter, self-exclusion ("what lowers C" never answers "pin C"),
pin-set dedup, the ancestry filter (an entity with no path to the target is
reachable=false and unsolved), cap reporting, and the honest zero ("0 of N
move this readout here").

**Additive parity (P4, 2026-11-16)**: the g4 replay through the feature branch
(the inverse endpoint registered and live): all 24,356 recorded folds
reproduce (max |diff| 5e-7, zero past 1e-5; `parity_folds_158.tsv`).
