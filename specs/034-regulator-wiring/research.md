# specs/034 — How catalysts and regulators are wired: inventory

**Status:** analysis, a design input. Nothing is proposed for adoption here.

**Why (Adam, 2026-09-26):** "we might need to put a lot more thought to how
these regulators are connected ... there are a lot of potential
configurations." Four findings in a row traced to regulator wiring:
- set members wired as separate required inputs (specs/033);
- an arbitrary copy chosen for joins (specs/030, STAT1);
- drug complexes as inhibitors (specs/032);
- recycled catalysts closing loops (RAF).

Before designing, this inventory counts which configurations occur and which
cases route through them.

**Tool:** `bench/analysis/regulator_inventory.py --catalog BUILD --cases ARM`.
- Build `20260926-1221_590301c` (canonical); cases from arm
  `results/039ce74/poolbase`.
- The role comes from the Regulation subclass in Neo4j.

## 1. Annotations (Neo4j, 92 pathways)

| role | protein | complex | complex with set | set | set of complexes | small molecule | total |
|---|---|---|---|---|---|---|---|
| catalyst | 263 | 406 | 427 | 156 | 108 | 0 | 1,361 |
| positive gene-expression regulation | 42 | 141 | 120 | 4 | 4 | 1 | 312 |
| negative regulation | 46 | 59 | 108 | 13 | 9 | 2 | 237 |
| positive regulation | 28 | 50 | 38 | 15 | 4 | 12 | 147 |
| negative gene-expression regulation | 4 | 53 | 64 | 1 | 3 | 0 | 125 |
| requirement | 3 | 1 | 10 | 1 | 1 | 9 | 25 |

**The generator has three roles** (`_CATALYST_CYPHER`, `_POS_REG_CYPHER`,
`_NEG_REG_CYPHER`), and Neo4j's labels fold the rest into them:
`Requirement` and `PositiveGeneExpressionRegulation` both carry
`PositiveRegulation`, and `NegativeGeneExpressionRegulation` carries
`NegativeRegulation`. So a Requirement ("needed, not sufficient") is wired
exactly like a positive regulator: an AND input whose rise multiplies the
reaction.

## 2. Edges in the network, by role and wired copy

`copy` is the SOURCE node of the edge:
- **root:** no incoming edge;
- **upstream:** produced, and not reached from the reaction;
- **downstream:** produced by something the reaction itself reaches, so the edge
  closes a loop.

| role | edges | root | upstream | downstream |
|---|---|---|---|---|
| catalyst | 26,457 | 20% | 48% | **32%** (8,481) |
| negative | 1,954 | 16% | 82% | 2% |
| positive | 1,948 | 16% | 83% | 1% |
| positive gene expression | 409 | 8% | 88% | 5% |
| negative gene expression | 158 | 4% | 95% | 1% |
| requirement | 118 | 29% | 42% | 29% |

**A third of all catalyst edges are fed by a copy the reaction regenerates
downstream.** That is catalytic recycling (the loop taxonomy of 2026-06-12),
and it is rare for regulators. By shape:

| catalyst shape | edges from a downstream copy |
|---|---|
| set | 6,575 |
| complex with set | 1,435 |
| set of complexes | 369 |
| complex | 89 |
| protein | 13 |

The cells with the most edges are catalyst / complex with set / upstream
(7,210), catalyst / set / downstream (6,575) and catalyst / set / root (3,124).

## 3. Cases routed through each configuration

A case counts toward a cell if its pinned nodes reach the edge's source and the
edge's reaction reaches its readout. **A case can count in many cells, so these
accuracies are associations, not causes.** They point at where to trace.

**Low cells (n ≥ 40):**

| axis | configuration | n | accuracy |
|---|---|---|---|
| curator | positive / complex with set / downstream | 240 | **42.9%** |
| curator | catalyst / set of complexes / root | 30 | 40.0% |
| curator | negative / protein / downstream | 42 | 35.7% |
| curator | catalyst / set of complexes / downstream | 316 | 55.7% |
| curator | positive gene expression / complex with set / downstream | 280 | 58.2% |
| curator | requirement / complex with set / upstream | 46 | 63.0% |
| curator | negative / complex / upstream | 312 | 65.1% |
| curator | negative / complex with set / upstream | 948 | 67.8% |
| experimental | catalyst / set of complexes / downstream | 43 | **20.9%** |
| experimental | negative gene expression / complex with set / upstream | 56 | 37.5% |
| experimental | negative / complex / upstream | 62 | 46.8% |
| experimental | negative / complex with set / upstream | 87 | 48.3% |
| experimental | positive gene expression / complex with set / upstream | 97 | 56.7% |

**High cells**, for contrast: catalyst / set / upstream (curator 91.3%, n = 738),
catalyst / set of complexes / upstream (experimental 96.4%, n = 165),
negative / set / root (curator 100%, n = 60).

## 4. Mismatches between wiring and Reactome's meaning

| # | mismatch | size | status |
|---|---|---|---|
| M1 | Requirement wired as a multiplying positive regulator | 25 annotations, 118 edges | not tested |
| M2 | set-valued catalysts and regulators as separate AND/OR member terms | 1,881 + 56 + 146 groups | specs/033: pooling is neutral-to-positive; draft LNG #99 |
| M3 | catalyst fed by a downstream (recycled) copy | 8,481 edges, 32% of catalysts | `DS_SCC_BREAK_CATALYST` (TP53 +104, held-out −65) and specs/018 roles: not adopted |
| M4 | inputs fed by a downstream copy (RAF F-actin, CNKSR2, Ca2+) | not yet counted | specs/033 result: traced, open |
| M5 | drug-derived regulators propagate | 202 entities, 24 pathways | specs/032: inert not adopted |
| M6 | negative regulators that are complexes (with or without sets) | the low cells in §3 | not traced |
| M7 | gene-expression regulation wired like reaction regulation | 437 annotations | semantics plausible; low experimental cells, not traced |

**Next:** trace one case from each low cell in §3 (M6 and M7 first, because they
are the least explained), and count M4. Then propose a rule set for review.

## 5. Traced: M6 (negative regulators that are complexes)

`regulator_inventory.py --cell negative/complex/upstream`. The failed cases
concentrate on five regulator → reaction pairs:

| pathway, perturbation | inhibitor | inhibited reaction | mechanism |
|---|---|---|---|
| RAF (KRAS/HRAS/NRAS up, NF1 KO) | BRAP:KSR1:MARK3 | "RAS:GTP:'activator' RAF homo/heterodimerizes …" | inside the RAF loop (see §6) |
| PIP3, PDPK1 KO | **N:M:PDPK1**: PDPK1 bound by **SARS-CoV-2** N and M | "PDPK1 phosphorylates AKT at T308" | self-contained, so the KO also removes the inhibitor and cancels itself; the virus is exogenous |
| ERBB2 up | **ERBB2:trastuzumab**:ERBIN:HSP90:CDC37 | "Trans-autophosphorylation of ERBB2 heterodimers" | self-contained and a **drug** |
| Activin (INHBA, ACVR2A) | INHIBIN-A:TGFBR3:ACVR2A | "Activin … binds Activin Receptor ACVR2A,B:ACVR1B" | self-contained sequestration (specs/022) |
| WNT (WNT1, WNT5A) | pT298-NLK dimer | "TCF/LEF:CTNNB1 bind canonical WNT target promoters" | not traced |

**Exogenous entities** are not in the cell a benchmark case describes:

| kind | entities | pathways |
|---|---|---|
| drug-derived | 202 | 24 |
| contains a non-human-species protein | 66 | 9 |

The pathogens are Influenza A, RSV, Rotavirus, HCV, Measles, SARS-CoV-2,
N. meningitidis and C. trachomatis. DDX58/IFIH1 has 45 of the 66, and PIP3 has
the SARS-CoV-2 complex above.

specs/032 held drugs inert: held-out **+30 / 0 broken**. It failed its gate
only on RAF's experimental cases, where the loop in §6 dominates.

## 6. Counted: M4 (inputs supplied only by the reaction's own downstream)

| where an input edge's source comes from | edges | share |
|---|---|---|
| a root | 38,599 | 43.7% |
| upstream | 41,047 | 46.4% |
| **downstream only** (every producer is reached from the reaction) | **7,989** | **9.0%** |
| downstream, plus an outside supply | 746 | 0.8% |

- This happens in 70 pathways: Class I MHC 1,098, DNA Repair 525, DNA Damage
  Bypass 422, HOX 421, DSB 407, EGFR 407, Cell Cycle Checkpoints 399, HDR 398.
- **87% (6,948) have no other usable copy** of the entity in the network; 11%
  have a root copy and 3% an upstream one.
- The entity exists only inside its own cycle, so nothing outside the cycle
  anchors its baseline. That is RAF's F-actin, CNKSR2 and Ca2+ (supplied only by
  "Dissociation of RAS:RAF complex"), and the knife-edge loops of specs/013 and
  014.
- MP-BioPath cut such loops by hand, which made the entry entity a root.
  `DS_PIN_SCOPE=root_cycle` reproduces that only for the PERTURBED gene.

## 7. Design candidates this points to (for review, none built)

- **D1 — anchor cycle-only entities.** An input whose every producer is
  downstream of its consumer gets a baseline supply: a root copy (fold 1)
  OR-pooled with the recycled copy. This is the curated entity's initial
  supply, the pool a cell starts with.
  - Mechanism: the any-copy pool (specs/031, 033) with a synthetic root member.
  - It removes the all-zero trap without cutting a curated edge.
  - Addresses M4 and RAF, and plausibly M3, since recycled catalysts are the
    same shape.
- **D2 — exogenous entities inert.** Drugs AND pathogen-derived entities are
  held at baseline (specs/032's rule, widened to species). Addresses M5 and
  part of M6.
- **D3 — Requirement as a limiter.** A knockout blocks the reaction; a rise
  does not raise it, so min(fold, 1). Addresses M1.
- **D4 — set pools (specs/033),** already measured neutral to positive.

D1 is the largest and the most structural. D2 and D3 are small, principled and
independent. Each is to be pre-registered and measured on its own before any
combination.

## 8. RAF traced by iteration (2026-09-27)

HRAS up (80x) on the canonical build. The solve was snapshotted at iteration
budgets 1, 2, 3, 5, 10, 20, 50 and 535. Fixed-point iteration is deterministic
within a process, so each budget continues the same trajectory.

| iteration | RAS:GTP:RAF complex (input) | activated RAF:scaffold:MAP2K:MAPK | 128-node component below 0.5x |
|---|---|---|---|
| 2 | 1.04x | **0.98x** | 0 |
| 3 | 1.16x | 0.92x | 0 |
| 5 | 1.57x | 0.60x | 4 |
| 10 | 3.76x | 0.04x | 115 |
| 20 | 100x | 6e-5 | 127 |

**The up-signal turns into a down-signal at the step that forms the activated
scaffold complex.** That step is "MAP2Ks and MAPKs bind to the activated RAF
complex" (R-HSA-5672972). At iteration 3:

- **The push:** negative regulator R-HSA-5675413, p21 RAS:GTP:activated RAF1
  homo/heterodimer:PEBP1 (RKIP), rises with RAS (1.154x). It contains the
  reaction's own input, the self-contained shape of specs/022.
- **The amplifier:** in our network the reaction has **27 AND inputs**, about 20
  of them recycled from the downstream "Dissociation of RAS:RAF complex" (each
  0.994x). Their product is 0.994^20 ≈ 0.89 per pass: a loop gain of about 20
  on any step down, so it collapses instead of settling.

**Reactome curates this reaction with 4 inputs:**
- p21 RAS:GTP:activated RAF;
- MAPKs (a DefinedSet of 2);
- MAP2K homo/heterodimers (a DefinedSet of 3);
- RAF/MAPK scaffolds (a **CandidateSet of 8**).

That is at most 2 × 3 × 8 alternatives. The network has **one reaction node
that requires every member of every set at once**, including the leaves of all
8 scaffold candidates (F-actin, TLN1, VWF, fibrin, integrins, CNKSR1/2, KSR2,
IQGAP1, ARRB1/2 …).

**How common this is:**
- Counted by leaves unique to each member (members sharing a subunit are not
  counted twice), **194 reaction copies in 16 pathways** require at least 2
  alternatives of one input set: DDX58/IFIH1 96, Class I MHC 33, Neurexins 13,
  RAF 13, RET 12.
- A looser count (any shared leaf) read 3,709. That count was an artefact of
  shared subunits and is withdrawn.

**Hypothesis under test:** the variant cap (`LNG_MAX_VARIANTS`, 512) bundles an
over-cap reaction's alternatives into one node, which is an AND of all of them.
RAF is being regenerated with the cap at 512 and at 10^6 to test it.

**This also explains the three failed RAF fixes.** Drugs held inert, set
pools, and conserved participants each removed a secondary term. None of them
removed the ~20-way AND of alternatives that turns a small step down into a
collapse.

## 9. RAF: MP-BioPath's network against ours (2026-09-27)

MP-BioPath's RAF network (`~/codes_and_results/mp_biopath/V_78_direct/pathways_2/RAF_MAP_kinase_cascade.tsv`,
Reactome v78, loops cut by hand) scores **100%** experimental with DeltaSignal's
own solver (docs/RESULTS §1b). Ours scores 29%. Same solver, so the difference is
the network:

| | MP-BioPath | ours (v97) |
|---|---|---|
| size | 842 nodes, 1,077 edges | 1,571 nodes, 4,762 edges |
| largest loops | **12 and 6 nodes** | 138 and 128 |
| sets | **one node per set, members OR-ed at the curated-alternative level** ("RAF/MAPK scaffolds" ← its 8 candidates; "RAF activating kinases" ← 7 kinases) | expanded: input sets into reaction copies (and past the variant cap into all-required bundles); catalyst and regulator sets into every member on every copy |
| "MAP2Ks and MAPKs bind to the activated RAF complex" | 4 inputs (the 3 set nodes and RAS:GTP:RAF), inhibitor **free PEBP1** | 28 member-level inputs, inhibitor PEBP1 bound to RAS:GTP:RAF |
| "Phosphorylation of RAF" catalyst | **one** set node | 9 members, all required |
| "RAF phosphorylates MAP2K dimer" | the complex once, plus 2 drug inhibitors | the complex as input AND catalyst (squared) |
| "Dissociation of RAS:RAF complex" outputs | p-MAP2K dimers, MAPKs, 5672712: **no scaffold release** | releases every scaffold leaf (F-actin, CNKSR, talin, …), which recycle into the binding step |
| RAF's loop | the real MAP2K → RAF feedback, through the OR catalyst node | welded through ~20 recycled scaffold leaves |

**Implication.** Every patch tried here (specs/032, 033, 035, 036, 037) tried to
re-create inside our representation what MP-BioPath's representation has
natively: one node per curated set, OR over its alternatives, and no member-level
fan-out.

The generator's design is the opposite. Sets "should always decompose to member
species, they should NOT survive as network nodes" (set-expansion fix, July
2026). That choice is what produces copy multiplicity (specs/031), capped
bundles (specs/036) and member fan-out (specs/033).

**Whether to revisit it is a design decision for Adam, not an experiment to
run unasked.** The scaffold release on dissociation is a curation difference
(v78 against v97), to be checked in Neo4j before it is attributed to the
generator.

## 10. The MP-BioPath hand edits (2026-09-27)

**Source:**
- `~/gitroot/mp-biopath-test-neo4j-generated-pathways/PathwayAnalysis/QA/*_changes.tsv`
  (87 files; columns source, destination, interaction type, and/or, action,
  comment), applied by `bin/add_changes_to_pathway.pl`.
- Unedited: `reactome_exports/`. Edited: `pathways/`.
- Adam: "some of these loops/cycles were in reactome and were intentional, some
  others were loops caused by the fact that a single entity in a compartment has
  the same id even if it is in multiple places."

**Totals:**
- **1,800 edits in 71 pathways:** 1,171 deletions, 629 additions.
- Cyclic nodes fell from **1,007 (export) to 413 (edited)**, −59%.
- Comments are rare: "set to member link" ×4, "…with Active AKT", "receptor:ligand
  set to member links".

| action | edge (MP-BioPath node kinds) | n | inside a loop in the export |
|---|---|---|---|
| delete | entity → entity (set-to-member / same-id links) | 607 | 55 |
| delete | entity → reaction (inputs) | 163 | 12 |
| delete | **reaction → entity (outputs)** | 148 | **84** |
| delete | entity → set-decomposition node | 96 | 1 |
| delete | set-decomposition → entity | 76 | 22 |
| delete | **reaction → set-decomposition node** | 63 | **38** |
| delete | negative entity → reaction | 13 | 0 |
| add | positive entity → entity | 182 | – |
| add | **negative entity → entity** | 131 | – |
| add | **negative entity → reaction** | 87 | – |
| add | positive reaction → entity | 78 | – |
| add | positive entity → reaction | 36 | – |
| add | positive, involving an unnamed (synthetic) node | ~116 | – |

**Three separable kinds of edit:**
1. **Loop cutting.** About 120 of the 208 in-loop deletions remove a reaction's
   output into a recycled species or a set's members. RAF: "Dissociation of
   RAS:RAF complex" → the MAPK set's member-level AND/OR decomposition was
   deleted and replaced by an edge to the set node.
2. **Same-id and set-to-member artefacts** (most of the 607 entity → entity
   deletions, which are mostly outside loops). This is the class the
   generator's per-position UUIDs are meant to handle.
3. **218 negative edges added** that Reactome does not curate, 131 of them
   entity → entity. Our networks have none of them. This is a direct source of
   difference that no generator change can recover from Neo4j alone.

**Next:**
- Map the loop-cutting deletions (class 1) onto our v97 networks: do the same
  recycling edges close loops there? That would give a curator-made label set of
  "artefact loop edges" to design loop handling against, instead of guessing.
- List the 218 added inhibitions and ask whether each is biology Reactome lacks
  or a modelling patch.

## 11. Our loops are mostly curated reaction cycles, not artefacts (2026-09-27)

**MP-BioPath's loop-cutting deletions mapped onto our v97 networks** (canonical
`20260926-1221_590301c`): 122 reaction → output edges were deleted INSIDE a loop,
in 25 pathways.

| where the deleted edge lands in our network | edges |
|---|---|
| present, and **not** in a loop | **85** |
| present, and still closes a loop | 14 (DSB/HDR 5687758 → 68462; Intrinsic Apoptosis BAD / p-S99-BAD, 139904, 141643, 350870; Mitotic G1 69195 → 68377) |
| no such edge | 20 |
| pathway not in our catalog | 3 |

**So the per-position UUIDs already resolve about 86% (85 of 99) of the loop
artefacts MP-BioPath cut by hand.**

**What our 8,064 cyclic nodes are made of** (nodes that stop being cyclic when
one edge type is removed):

| edge type removed | cyclic nodes no longer cyclic |
|---|---|
| output | 8,058 (100%) |
| input | 7,878 (98%) |
| depletion | 1,426 (18%) |
| regulator | 1,033 (13%) |
| catalyst | 487 (6%) |
| dissociation | 0 |
| assembly | 0 |

**Most cyclic nodes sit in loops made of curated input/output edges alone**:
Reactome's own reaction cycles (interconversions such as RAS GDP ⇄ GTP,
phosphorylation ⇄ dephosphorylation, association ⇄ dissociation). Derived edges
contribute a minority. This supports Adam's position: keep the loops, and
handle them mathematically.

**The math they need.** A cycle A → B → A driven by opposing enzymes (GEF / GAP,
kinase / phosphatase) has a CONSERVED total whose split is set by the two
activities: B/A ∝ forward / backward. The solver instead multiplies around the
cycle: gain 1 at baseline, the knife-edge of specs/013 and 014, and the RAF
collapse. specs/017's loop pool used the product of entries, which is not
this ratio. **A ratio-based treatment of interconversion cycles is the
untested candidate.** It needs a design pass (identifying interconversion
cycles and their opposing catalysts from the curated reactions) before any
pre-registration.

## 12. The 219 negative edges MP-BioPath added (2026-09-27)

Adam: do the work first; the curators can then be asked what was done and
whether it can be replicated automatically in the new networks.

**Patterns** (by the schema classes of source and target, from Neo4j):

| pattern | examples |
|---|---|
| **reversal enzyme ⊣ the modified form, or ⊣ the forward reaction** (the largest group) | phosphatases: PPP5C, PTPN6 / PTPN11 ⊣ "Phosphorylation of INFAR1 by TYK2", myosin phosphatase ⊣ p-MRLC; demethylases: KDM1A/B, KDM2A/B, KDM4A, KDM5A-D, KDM6B, JMJD6 ⊣ histone methylation reactions; deubiquitinase OTUD5; **RAS GAPs, SPRED:NF1 ⊣ p21 RAS:GTP** |
| ubiquitin ligase ⊣ its target | ITCH, RNF125:E2, CUL1:SKP1:SKP2:CKS1B |
| sequestering trap or decoy ⊣ ligand | FST, FSTL3 ⊣ Activin; "Ligand Trap" ⊣ BMP2; IL13RA2 ⊣ "IL13 binds IL13RA:TYK2"; I-SMAD ⊣ BMP:BMPR |
| transcriptional repressor ⊣ target-gene expression | RB1 ⊣ E2F1-driven TK1, CCNA1, DHFR, TYMS expression |

**Which of them our v97 network already has** (canonical build; the source and
target entity with a negative edge between them):

| status | edges |
|---|---|
| same negative edge exists (a `depletion` edge) | 7 |
| the source is a negative source in ours, but with a different target | 9 |
| **absent** | **187** |
| pathway not in our catalog | 16 |

**Reading.** The dominant pattern is the **backward enzyme of an
interconversion cycle** (GAP, phosphatase, demethylase, deubiquitinase) used as
a brake on the forward-modified species. It is the hand-made counterpart of §11:
in a cycle A ⇄ B the split B/A follows forward / backward activity, and
MP-BioPath encoded "backward" as an inhibition of B.

- Our generator derives a small part of it: depletion edges for phosphatases
  (Pi as output) and Ub ligases (Ub as input).
- It does not derive it for demethylases, GAPs (GTP hydrolysis), deubiquitinases
  or transcriptional repressors.
- The traps and decoys are ordinary sequestration.

**Automatic equivalent to test:** for each curated reversal reaction B → A with
catalyst E_b, make E_b act against B, or better, give the A ⇄ B cycle its
ratio-based treatment (§11).
- This is derivable from Reactome: pairs of reactions whose input/output
  entities are the modified/unmodified forms of one protein.
- It is the pattern to confirm with the curators.

## 12b. Correction: how much of the added-inhibition pattern Reactome's curation supports (2026-09-27)

§12 called the added inhibitions "mostly reversal enzymes", **from their
names**. Tested formally against Neo4j (does inhibitor E catalyse a curated
reaction that CONSUMES the target X, or that consumes the target reaction's
output?):

| relationship in the curated reactions | edges |
|---|---|
| E catalyses a reaction consuming X (entity target) | 37 |
| E catalyses a reaction consuming the target reaction's output | 16 |
| **total with a curated reversal link** | **53 (24%)** |
| **no catalytic link** (entity target 89, reaction target 67) | **156 (71%)** |
| E catalyses X's production | 1 |
| not in this release | 9 |

**So only about a quarter can be derived automatically from curated
catalysis.** For the rest, the inhibitor LOOKS like a reversal enzyme
(phosphatase, demethylase, GAP), or is a trap or a repressor, but Reactome
does not curate the reverse reaction with that enzyme in the pathway. §12's
"automatic equivalent" therefore covers about 24% of what the curators did.
The remainder is exactly the question for the curators: what rule, or what
outside knowledge, did they apply?
