# specs/031 — Anatomy of the severed routes left after specs/030

**Status:** analysis only. No intervention is adopted or pre-registered here.

**Question (Adam, 2026-09-26):** after the complexes are built (specs/030,
step 1), how do we fix these, in order?
- step 2: "the protein has no edge to the complex";
- step 3: "the missing edges to reaction";
- step 4: better pinning.

## Input

- Build `20260926-0917_8563de9_hier2`, arm at solver 12d6e11: the specs/030
  re-measure.
- Curator held-out failures where Reactome has a strict route and ours has
  none: **252 cases**.
- `bench/analysis/route_breaks.py --faithful` walks Reactome's route using only
  the hops the generator can represent. `case_oracle.faithful_comp_pairs` keeps
  `hasComponent` hops into root complexes and their nesting, which is what
  specs/030 decomposes.
- A route that needs a hop into a **produced** complex is
  **composition-only**.
- The oracle tools are now deterministic (commit 7238206). The route tie-break
  had followed the string hash seed, and read 186 or 166 "edge missing" on the
  same input.

## First unreached step on the faithful route

| cases | first break | route |
|---|---|---|
| 80 | entity → reaction, node exists but edge missing | faithful |
| 62 | component → complex, node exists but edge missing | faithful |
| 62 | component → complex, node exists but edge missing | composition-only |
| 36 | route starts from another form of the gene | faithful |
| 6 | whole route reached | faithful |
| 6 | other (node absent, or entity → reaction) | mixed |

## Step 2: "the protein has no edge to the complex"

Walking Reactome's own shortest route, 186 cases break at component →
complex. Classified by the complex:

- **106: the complex is never formed in Reactome.** Its only producers
  consume something that already contains it (recycling or dissociation
  steps), or it has no producer at all.
  - PDGF:Phospho-PDGF receptor dimer (R-HSA-186811), 54 cases. Its only
    producer is "Activated PLC gamma dissociates from the PDGF receptor".
  - p-2S-SMAD3:p-2S-SMAD3:SMAD4 (R-HSA-8878153), 10.
  - NODAL:p-NODAL Receptor (R-HSA-1181134), 8.
  - two SSA-repair complexes (R-HSA-9980019, R-HSA-9980022), 12.
  - HOX bivalent chromatin (HOXD1/D3/B3/C4, HOXA4), 14.
  - SLIT1:ROBO1, BOC:PTCH1, IL receptor:SHC1:SHIP1 and RIG-I/MDA5:…:Casp-8, 8.

  The curators reached the readout through `hasComponent` membership, not
  through any reaction. The generator can reach these only with composition
  edges into produced complexes, which were measured harmful twice
  (specs/016, specs/029) because they close cycles. **They are candidates for
  curation review** (a complex used by reactions and formed by none), not a
  generator lever.
- **80: the complex is formed in the pathway.** These split further:
  - 26: the copy of the component the gene reaches is **downstream** of the
    complex (HOX KMT2D/EP300, DSB RAD52). Joining it would weld a cycle, so
    the specs/030 guard is right to refuse. Traced for DSB: the RAD52
    heptamer copy the KPNA2 perturbation reaches is itself downstream of the
    complex.
  - 6: the reached copy is upstream but was not joined (DSB MUS81, Hedgehog
    PTCH). This is the case an "any eligible copy" pool would fix, together
    with the 25 STAT1 cases the root-copy preference cost (specs/030).
  - the rest reach the complex along a reaction route whose break is elsewhere.
    Walked faithfully, they appear under step 3.

## Step 3: "the missing edges to reaction"

All **80** entity → reaction breaks on the faithful route have one mechanism.
The gene reaches a copy of the reaction's curated input entity, but **not the
copy that feeds the reaction** (Hedgehog 24, DSB 22, apoptosis 12, RUNX 8,
FGFR1/3/4 12, CD28 2).

- **None is a missing curated edge.** Every reaction has its input edges;
  they come from another copy of the same entity.
- **Traced (Intrinsic Pathway for Apoptosis):** BCL2 [mitochondrial outer
  membrane] (R-HSA-50757) has 7 copies. STAT3 → "BCL2 gene expression"
  (R-HSA-6790025) → a BCL2 copy with no outgoing edge. "Sequestration of tBID
  by BCL-2" (R-HSA-114352) reads a different, root copy.
- **Traced (FGFR4):** "Activated FGFR4" is a set, represented by its members.
  "Activated FGFR4 binds FRS2" exists as 13 variant copies, each fed by one
  member form. The perturbation reaches other copies of those forms.
- The generator links an output to an input only along `precedingEvent` or
  diagram connectivity, not by shared entity.
- **This is the sink-bridge / instance-multiplicity family**
  (`project_gap_anatomy_2026_09_18`, specs/015, specs/019). Bridging
  same-entity copies was measured three times, null or negative:
  - all-pairs silo bridge −73 / −77;
  - sink bridges with a cycle guard, held-out ≈ 0, false change +171;
  - the median fan-out is 139 copies per entity.

  **Not re-tried here without a new mechanism.**

## Step 4 (pinning)

36 cases start from another form of the gene (PTK6 26). They are not analysed
here yet.

## What is left that is new

- **The any-eligible-copy pool (step 2, about 31 cases):** all eligible copies
  feed one pool node by OR (mean), and the pool feeds the complex as one AND
  input, so nothing is double-counted. Small.
- **The 106 never-formed complexes** are a list for Reactome curation, not for
  the generator.
- Everything else found in steps 2 and 3 is a mechanism already measured
  closed.
