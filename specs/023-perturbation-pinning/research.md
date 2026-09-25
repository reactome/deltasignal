# 023 — Pin where the gene enters the network, and let the rest propagate

## What the benchmark does today (measured 2026-09-25, build `20260925-1039_d4f4f64`)

A perturbed gene resolves to the entities that reference it (Neo4j
`referenceEntity`, `neo4j_gene_to_stids`). **Since 2026-07-14 (`0bd4565`)** each
of those stids also matches every node whose `member_leaves` include it, so the
benchmark pins **every complex that contains the gene, wherever it sits**:

- **Scale:** 856 perturbations pin 9,378 nodes. Only **1,023 (11%) are roots**
  (in-degree 0).
- **Mid-pathway pins:** 8,355, made up of 6,371 complexes, 1,333 dissociation
  sinks and 651 entity forms produced by a reaction (e.g. a glycosylated
  secreted form).
- **Coverage:** 805 of 856 perturbations pin at least one mid-pathway node. 77
  have no root occurrence at all.
- **Example:** WNT5A in Signaling by WNT pins 1 root entity and 36 downstream
  nodes. These include generic set complexes (`WLS:WNT`, `WIF1:WNT`), and
  pinning them at 0 knocks out *every* WNT ligand. That is the specs/021
  reversal.

A hard pin overrides the node's own producing reaction. A downstream complex
is then not *computed* from the perturbed gene; it is *set*, including for
partners the knockdown should not affect.

## Intended protocol (Adam, 2026-09-25)

*"the idea is to pin the root input that is either a particular entity or a
complex containing that entity … breaking apart the complexes into root
entities where the root would look like A entity and B entity -> AB complex ->
reaction. In this way both the entity that is not part of the complex and the
entity inside the complex can be perturbed."*

## The rule (`DS_PIN_SCOPE=entry`; default `all` = today)

Let G be the nodes the gene resolves to today. Pin only the **entry
occurrences**: the members of G that cannot be reached from any other member
of G. Roots in G are always entries.
- A gene first produced mid-pathway (e.g. by transcription) gets that first
  occurrence pinned.
- When every member of G is reachable from another (G lies inside one cycle),
  all of G is pinned, as today.
- Everything else propagates.

## Pre-registration (committed before the arm runs)

- **Arms:** control = production scoring (`ae84de9`, byte-identical to the
  `ffb3aa1` control in specs/022); arm = `DS_PIN_SCOPE=entry`, code defaults
  otherwise.
- **Decides:** curator held-out net and macro-F1, McNemar p, pathways moved,
  distinct genes. **Reported:** tuning, experimental axis, per-pathway net,
  pinned-node counts, and the WNT5A dose ladder.
- **This is a protocol correction, not a model change.** The primary question
  is whether the benchmark measures what it claims to, so it is adopted if the
  held-out result is not worse than the noise floor (−15). Adam's statement of
  intent is the justification, not the accuracy. A gain is reported, not
  required. A loss beyond the noise floor is traced before any decision.
- **Prediction:** large churn, because 94% of perturbations change their pin
  set. Direction unknown. The false-change rate should fall, because
  downstream complexes stop being set directly.

## Result, arm 1: `DS_PIN_SCOPE=entry` alone — a large loss, traced

Build `20260925-1039_d4f4f64`, bench `ec81afd`, solver unchanged; control =
production `ae84de9`. Same 23,268 valid cases in both arms.

| | control | entry pinning | net | fixed / broke | p |
|---|---|---|---|---|---|
| curator held-out | 0.8687 / mF1 0.8299 | 0.8281 / mF1 0.7617 | | | |
| curator all | 0.8477 / 0.8145 | 0.7994 / 0.7401 | **−1,164** | 179 / 1,343 | 1e-220 |
| experimental | | | **−137** | 3 / 140 | 9e-38 |

- **Breadth:** 71 pathways got worse and 2 got better.
- **Direction:** 1,149 of the 1,343 breaks (86%) are overexpression cases,
  and 1,148 are UP → NORMAL.

**Traced (KMT2C OE, Chromatin modifying enzymes; readout 100 → 1.0).**
- KMT2C, ASH2L, RBBP5, WDR5 and DPY30 feed the MLL3 complex by `assembly`
  edges. This is exactly the decomposed-root structure the protocol intends.
- The old protocol pinned the **MLL3 complex itself** at 80x.
- Under entry pinning only KMT2C is pinned. `DS_ASSEMBLY_LIMITING=1` makes a
  complex its scarcest subunit: min(80, 1, 1, 1, 1) = 1. The overexpression
  cannot pass through any complex.
- **So the broad pin was compensating for a modelling choice.** It set every
  complex directly, which is how 89% of pins came to be mid-pathway, and how
  overexpression cases scored.

**Status:** the pinning is now what the protocol intends, but the model cannot
carry a single-subunit overexpression through an assembly. Adopting arm 1 alone
would make the benchmark honest and the model worse at the thing curators score
most.

## Pre-registration, arm 2 (committed before it runs)

`DS_PIN_SCOPE=entry` + `DS_ASSEMBLY_LIMITING=0`, so assembly inputs combine by
the AND mode (`hill_sat`, multiplication capped at 100) instead of min.
- A knockdown still pulls a complex down (0 × anything = 0).
- An overexpressed subunit now raises it.
- The specs/009 measurement of limiting off (−196) was taken under broad
  pinning, which masked exactly this effect, so it does not carry over.

**Reading:** compare against the production control on the same columns as
arm 1. If held-out is within the noise floor (±15) of production or better,
the honest protocol plus multiplication is adoptable as a pair, pending Adam's
call on the biology. That call is whether overexpressing one subunit should
raise a complex when its partners are at baseline. If held-out is below −15,
record it and trace.

## Result, arm 2: entry pinning + `DS_ASSEMBLY_LIMITING=0`

Same build, bench `ec81afd`. The container's own environment was verified as
`DS_ASSEMBLY_LIMITING=0` before the run.

| | production | arm 1 (entry) | **arm 2 (entry + multiply)** |
|---|---|---|---|
| curator held-out acc / mF1 | 0.8687 / 0.8299 | 0.8281 / 0.7617 | **0.8475 / 0.7952** |
| curator all acc / mF1 | 0.8477 / 0.8145 | 0.7994 / 0.7401 | 0.8281 / 0.7846 |
| held-out net vs production | — | | **−402** (88 / 490), 36 of 38 moved pathways worse |
| tuning net | — | | −69 |
| experimental net | — | −137 | **−51** (1 / 52) |

Letting an overexpressed subunit raise its complex recovers about 60% of arm
1's loss, but the honest protocol still scores well below production on both
axes. The remaining breaks are 339 OE and 151 KD on held-out. They are not yet
traced.

**What this means for every number reported since 2026-07-14.** Production
accuracy (held-out 86.87%) is measured under a protocol that sets every complex
containing the perturbed gene directly. Part of that accuracy is the pin, not
the propagation. Until the protocol question is decided, **quote the broad-pin
numbers with that caveat.** This affects docs/RESULTS.md and the manuscript
claims.

**Decision needed (Adam):**
1. **Which protocol is the benchmark of record?**
   - (a) Broad pins, as today: comparable to every result since July, but it
     sets complexes by hand.
   - (b) Entry pins: what the protocol intends, but about 2pp worse today, and
     it exposes the modelling gaps.
2. **Under (b), should an overexpressed subunit raise its complex** when its
   partners are at baseline? Limiting-reactant biology says no; the curators'
   expectations say yes.

Nothing is adopted until then. The flag stays, default `all`.
