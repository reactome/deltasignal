# specs/036 — Over-cap set bundles as one input per curated component

**Status:** pre-registered 2026-09-27, before any code or arm.

## The defect (traced in specs/034 §8)

The variant cap (`LNG_MAX_VARIANTS`, 512; issue #40) exists because expanding
every set combination is intractable. Uncapped, RAF alone asks for a
52.6 GiB (84,000 × 84,000) array.

**The problem is what the cap does:** when a reaction's input or output
combinations exceed it, `get_broken_apart_ids` merges every alternative of every
set into ONE member list, a single variant that **requires all of them**.

**RAF, "MAPs and MAPKs bind to the activated RAF complex" (R-HSA-5672972):**
- Curated with 4 inputs: RAS:GTP:activated RAF; MAPKs (a set of 2); MAP2K dimers
  (a set of 3); RAF/MAPK scaffolds (a CandidateSet of 8).
- The generator log: "variant cap hit for R-HSA-5672972: 12000 combinations >
  512; bundling 4 components into one node".
- The network node requires 28 inputs at once.
- The other members of the RAF loop are capped the same way: R-HSA-5672720,
  5672723, 5672724 and 5672980.
- About 20 of those inputs are recycled from the downstream dissociation, so a
  0.6% step down becomes 0.994^20 ≈ 0.89 per pass. The loop collapses to 0 under
  every RAS/RAF perturbation (experimental 14 of 49, against MP-BioPath's 46).

**Scale:** the canonical build's `generate.log` records **96 cap hits in 21
pathways**: RUNX1 24, HOX 19, WNT 6, CD28 6, RAF 5, Rho GTPases 5, and others.

## The rule

**`LNG_CAP_POOLS=1`:** at emission, for a reaction or complex whose inputs come
from a capped bundle, the bundled leaves are grouped by the CURATED input
component they belong to (from Neo4j: the reaction's inputs and their set
membership).
- A component that is a set becomes a **set pool** (specs/033): its leaves feed
  the pool by `set_member` edges, and the pool feeds the reaction as ONE AND
  input.
- A component that is not a set is wired as before.
- The reaction then has one input per curated component, as Reactome curates it.
- The matching layer (the cap itself, and the content hashes that line inputs up
  with outputs) is **unchanged**.

The solver reads pools with `DS_SET_POOL_MODE` (default `product`). For a pool
of alternatives, `product` would re-create the all-required bundle, so this
change is **measured with `extreme` and `mean`** as well. Which one to use is
chosen on the tuning split, as specs/033 did.

## Pre-registration

**Arms,** all through `scripts/run_arm.sh`:
- `base`: the canonical build `20260926-1221_590301c`, code defaults;
- `cap_product`, `cap_extreme`, `cap_mean`: one build with `LNG_CAP_POOLS=1` and
  `LNG_SET_POOL=1`, with each pool mode.

**Choosing the mode:** by tuning curator macro-F1 plus experimental macro-F1;
ties go to product < extreme < mean.

**Adopt** the chosen arm only if all hold:
- curator held-out net ≥ −15 and not significantly negative;
- experimental net ≥ 0;
- no single pathway loses more than 10 held-out cases.

It is a faithfulness fix (one input per curated component), so it does not need
a held-out gain to be adopted.

**Predictions:**
- RAF experimental rises well above 17 of 49 under `extreme` or `mean`, and the
  RAF solve converges;
- `cap_product` is close to base;
- the 21 capped pathways are the only ones that move, together with the
  specs/033 set-pool pathways, which are measured separately there.

**Also reported:**
- per-pathway net on both axes;
- concentration;
- convergence counts;
- the RAF iteration trace re-run on the new build;
- cyclic node count excluding pools.
