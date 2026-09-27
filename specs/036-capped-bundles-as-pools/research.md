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

## Sanity trace on RAF alone (2026-09-27, before any arm)

**Setup:**
- RAF (R-HSA-5673001) regenerated at LNG feat/set-pool + WIP `LNG_CAP_POOLS=1`,
  in the scratchpad (a diagnostic, not a catalog).
- "MAP2Ks and MAPKs bind to the activated RAF complex": **28 inputs → 4**. The
  4 are RAS:GTP:activated RAF plus pools of MAPKs (2), MAP2K dimers (3) and
  RAF/MAPK scaffolds (22 leaves of the 8 candidates).
- **Limitation:** that pool is over leaves, not candidates, so larger candidates
  weigh more.

HRAS pinned 80x or 0. Readout: p-T,Y MAPKs in the nucleus (R-HSA-5674340).

| configuration | HRAS up | HRAS KO |
|---|---|---|
| canonical | 0 (collapse) | – |
| cap pools, any pool mode | 0 (collapse), not converged | 100x (inverted) |
| + drug-derived nodes held at baseline | 0 (collapse) | 100x (inverted) |
| + `DS_DEDUP_ACTIVATORS=1`, `mean` pools | **0.83x**, no collapse (0 zeros), not converged | **0.016x** (correct direction) |
| + `DS_DEDUP_ACTIVATORS=1`, `extreme` pools | 0 (collapse) | 0 |

**Amplifiers found, in the order they appeared:**
1. Capped bundles: all alternatives required at once.
2. Drug-bound complexes inhibiting "RAF phosphorylates MAP2K dimer" and "MAP2Ks
   phosphorylate MAPKs" (R-HSA-9653109 and others rise with RAS).
3. **Catalyst = substrate squaring.**
   - "RAF phosphorylates MAP2K dimer" reads the activated RAF:scaffold complex
     as both input and catalyst.
   - "MAP2Ks phosphorylate MAPKs" does the same with its complex.
   - Each squares its own fold (specs/012; de-duplicating was held-out −15
     catalog-wide there).
4. **Open:** with 1–3 addressed and `mean` pools, up-regulation still reads
   0.83x and the component does not converge. Not yet traced.

**Consequence for the pre-registration:** the RAF prediction ("rises well above
17/49") cannot be met by cap pools alone. Before a combined arm (cap pools +
inert drugs + dedup) is pre-registered, amplifier 4 is traced. Each component
was individually neutral or negative catalog-wide, so the combination is
measured as a combination, and the RAF effect is not generalised from this
trace.
