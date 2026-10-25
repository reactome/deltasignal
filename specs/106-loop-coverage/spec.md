# specs/106: loop coverage (which loop types exist, which rule handles each)

**Status:** pre-registered analysis, 2026-10-25. Step 1 of the loop-type review (specs/096
section 2, continuing specs/097).

Adam, 2026-10-25: "I am starting to worry that we are over fitting our method to the tests. are
the things we are doing things that also improve other loops that are found in other pathways?
what percentage of the loops in other pathways are effected by our changes? ... Another thing we
could look into is working on other loop types in the other files. we would not have the
experimental evidence but I think we may find a way to say we have looked at every loop and feel
like we have methods that handle it well."

Decision recorded with it: no further rules are chosen by tracing experimental cases. New loop
rules are judged by a stated rationale and behaviour checks on networks without ground truth
(step 2), and accuracy is claimed only from a frozen external test (specs/096 step 3).

## Why specs/097 is not enough

097 counted loops on the network the solver runs, after the carrier-loop rules (specs/049, 077,
084) had already dropped their edges, so loops those rules dissolve were never counted, and it
did not say which rule handled which loop. It also predates specs/100–102.

## Method (`loop_coverage.jl`)

For each pathway, on the network AS GENERATED (every edge, recycle edges included):

1. the cyclic components (non-trivial strongly connected components): the **loops**;
2. each mechanism's handled edges, all as the default solver uses them:
   - `recycle`: recycle edges (specs/081), which the pool balance solves;
   - `carrier`: carrier-loop edges the default `closed_twin` + `fed_single_ub` rules drop
     (specs/049, 077, 088);
   - `pool`: an edge between two nodes of one interconversion, recycle or release pool
     (specs/039, 081, 086);
   - `self`: an inhibitor or activator slot built from its own reaction's input
     (specs/022, 040, 083, 095);
3. per loop, the mechanisms with at least one handled edge inside it (**touched**), and the
   cyclic components left inside it once every handled edge is removed and the solver's added
   edges (twins, RNA supply) are added (**remnants**), labelled as in 097 (`positive` by
   `catalyst` / `ubiquitin` / `material`, `self_depleting`, `negative`);
4. a loop is **dissolved** if no remnant is left, **partly handled** if touched with a remnant,
   **untouched** otherwise. A `negative` remnant is reported separately: negative feedback has a
   stable fixed point the default iteration finds, so it is a type the solver handles by design,
   and a `self_depleting` remnant is handled numerically (specs/066), not structurally.

Catalogs: the 92 benchmark pathways on canonical `20261024-2211_7f4b28a`, and the specs/097
sample of 101 other pathways rebuilt with the same generator commit (`--variant sample101`).

## Pre-registered expectations

- P1: the share of loops touched by a rule in the sample is within 10 points of the benchmark's.
- P2: carrier-loop rules dissolve more loops than any other single mechanism in both catalogs
  (097 could not see them).
- P3: material-closed positive remnants are the largest unhandled class in the sample by
  pathways, as in 097.

## Out of scope

No rule changes. Step 2 (behaviour checks per loop type on networks without ground truth) is a
separate spec.
