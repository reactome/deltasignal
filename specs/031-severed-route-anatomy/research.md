# specs/031 — Anatomy of the severed routes left after specs/030

**Status:** analysis only. No intervention is adopted or pre-registered here.

**Question (Adam, 2026-09-26):** after the complexes are built (specs/030,
step 1), how do we fix these, in order?
- step 2: "the protein has no edge to the complex";
- step 3: "the missing edges to reaction";
- step 4: better pinning.

## Input and method

- Build `20260926-0917_8563de9_hier2`, arm at solver 12d6e11 (the specs/030
  re-measure). The failures are the curator held-out cases where Reactome has
  a strict route and ours has none: **252 cases**.
- `bench/analysis/route_breaks.py --faithful` walks Reactome's route using only
  the hops the generator can represent, where such a route exists.
- `case_oracle.faithful_comp_pairs` defines those hops. They go into a ROOT
  complex (used by a reaction, produced by none) and into its unproduced
  nesting; a produced nested complex is a join, not a descent.
- A route with no such path is **composition-only**: it needs a `hasComponent`
  hop into a complex that no reaction forms.
- **Determinism:** the tools are now deterministic (commit 7238206). The route
  tie-break had followed the string hash seed, and read 186 or 166 "edge
  missing" on the same input.
- **Review revision:** a first version of `faithful_comp_pairs` also descended
  into produced nested complexes. It labelled 66 cases composition-only; the
  corrected rule labels 108. The numbers below are the corrected ones.
- Every class below was assigned **by script over all its cases**. The named
  traces are illustrations.

## First unreached step

| cases | first break | route |
|---|---|---|
| 96 | entity → reaction, node exists but edge missing | faithful |
| 34 | component → complex, node exists but edge missing | faithful |
| 8 | route starts from another form of the gene | faithful |
| 6 | whole route reached (readout mapping) | faithful |
| **144** | | **faithful subtotal** |
| 70 | component → complex, node exists but edge missing | composition-only |
| 30 | route starts from another form of the gene | composition-only |
| 6 | component → complex, node absent | composition-only |
| 2 | entity → reaction | composition-only |
| **108** | | **composition-only subtotal** |

## Composition-only (108)

- Reactome connects gene and readout only through `hasComponent` membership of
  a complex that no reaction of the pathway forms.
- The largest is PDGF:Phospho-PDGF receptor dimer (R-HSA-186811), 60 cases.
  Its only producer is "Activated PLC gamma dissociates from the PDGF
  receptor", which consumes a complex that already contains it.
- PTK6 accounts for 26 of the "another form" rows here.
- The generator could reach these only with composition edges into produced
  complexes, which were measured harmful twice (specs/016, specs/029) because
  they close cycles.
- A complex used by reactions and formed by none is a **candidate for
  curation review**.

## Faithful routes (144)

**component → complex (34):**

| cases | mechanism | examples |
|---|---|---|
| 16 | the reached copy is upstream, but a different copy was joined | NODAL:p-NODAL Receptor (R-HSA-1181134) 8, IL-2 family 2, Hedgehog BOC:PTCH1 2 |
| 10 | every copy of the complex in our network is produced | p-2S-SMAD3:p-2S-SMAD3:SMAD4 (R-HSA-8878153) 8, SLIT1:ROBO1 2 |
| 8 | the reached copy is downstream of the complex | HOX bivalent chromatin (KMT2D/EP300) |

- **The 16** are what an "any eligible copy" pool would fix, together with the
  25 STAT1 cases the root-copy preference cost (specs/030): about **41 cases**.
- **The 10** are not traced further.
- **The 8:** joining would weld a cycle, so the specs/030 guard is right to
  refuse.

**entity → reaction (96):** in every case the gene reaches a copy of the
reaction's curated input entity, but **not the copy that feeds the reaction**.
By pathway: Hedgehog 24, DSB 20, apoptosis 12, RUNX 8, FGFR1/3/4 12, and
others.

- **None is a missing curated edge.** Every reaction variant has its input
  edges; they come from another copy of the same entity.
- **Traced (Intrinsic Pathway for Apoptosis):** BCL2 [mitochondrial outer
  membrane] (R-HSA-50757) has 7 copies. STAT3 → "BCL2 gene expression"
  (R-HSA-6790025) → a BCL2 copy with no outgoing edge. "Sequestration of tBID
  by BCL-2" (R-HSA-114352) reads a different, root copy.
- **Traced (FGFR4):** "Activated FGFR4" is a set, represented by its members.
  "Activated FGFR4 binds FRS2" exists as 13 variant copies, each fed by one
  member form. The perturbation reaches other copies of those forms.
- The generator links an output to an input only along `precedingEvent` or
  diagram connectivity, not by shared entity.
- **This is the instance-multiplicity class** (`project_gap_anatomy_2026_09_18`,
  specs/015, specs/019). The nearest measured interventions were null or
  negative:
  - all-pairs silo bridge −73 / −77;
  - sink-to-sink bridging with a cycle guard, held-out ≈ 0 with false change
    +171;
  - a median fan-out of 139 copies per entity.

  A targeted join (the produced copy → the specific copy feeding each
  reaction that uses the entity) **has not been measured**. It is the untried
  variant of this class.

## Step 4 (pinning)

38 cases start from another form of the gene: 30 composition-only (PTK6 26)
and 8 faithful. Not analysed here.

## What is left that is new

1. **The any-eligible-copy pool (about 41 cases).** All eligible copies feed
   one pool node by OR (mean), and the pool feeds the complex as one AND
   input, so nothing is double-counted.
2. **The targeted produced-copy → feeding-copy join (up to 96 cases).**
   Unmeasured, and near a family that has failed three times.
3. **The composition-only list** (108 cases, 60 of them one PDGF complex),
   for curation review rather than generator change.
