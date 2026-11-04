# specs/133 research

## R1: the closure (census, no solving; 2026-11-03)

`bench/analysis/closure_133.py` (Neo4j only) follows specs/132's L1 handoffs back from the 92 scored pathways until a
round adds nothing. Output (private): `~/deltasignal-catalogs/analysis/133/` (`links.tsv`, `s_q.tsv`, `closure.tsv`).

- **275 pathways** are reached: the 92, then 83, 50, 22, 15, 7, 2, 1, 1 and 2 by round (nine rounds). 216 of them
  receive at least one handoff; 2,433 handoff rows over 407 entities (input 1,952, catalyst 352, positive regulator
  115, negative regulator 14). These rounds are the superset the builder draws from; each scored pathway keeps only
  its own tree.
- **Extracellular, never bridged: 26 rows, 18 entities**, every one a signal between cells or a solute outside them:
  serotonin, noradrenaline, L-glutamate, choline, LTC4, bile salts, complement C3a and C5a, angiotensin fragments and
  alamandine, GLP-1, acyl ghrelin, IL4, IL21; and H+, Na+, H2O2 outside the cell. The rule removes the IL4 class
  without naming it.
- 33 of the 407 entities are sets (79 members; `set_members_132.py` on this table).
- 100 of the 275 pathways were not in specs/132's base build; the base is rebuilt whole with the frozen generator
  (`8b6527e`, tag `freeze-2026-10-26`) as variant `ctx133`.

## R2: the bundles (census, no solving; 2026-11-03)

Base build `20261103-2155_8b6527e_ctx133` (275 pathways, all built; not `current`). Bundles from it with
`bench/context_build.py` (commit 910c5a7), each over the 92 scored pathways:

| bundle | rule | pathways with a bridge | bridges |
|---|---|---|---|
| `.c0` | no links (each pathway alone; the control) | 0 | 0 |
| `.A` | specs/132 L1, extracellular rows removed | 25 | 72 |
| `.T` | `--transitive` (the tree) | 25 | 163 |

- **A is specs/132 L1 minus exactly the two cell-to-cell bridges** (L1: 27 pathways, 74 bridges): IL4 into *Interleukin-4
  and Interleukin-13 signaling* and IL21 into *Interleukin-2 family signaling*. Nothing else in the one-hop set is
  extracellular.
- **T joins the same 25 pathways, not more.** A pathway only gets a tree when one of its own roots receives a handoff,
  and the closure's depth adds sources to those, never new receivers. Of 2,109 handoff rows the builder examined, 270
  are bridged (by depth: 115, 76, 54, 13, 12); 1,540 land on a copy that is not a root in the receiving pathway (the
  entity is there but made inside it), 232 come from a pathway already joined at a lesser depth, 59 carry a held entity
  (cofactor, drug, pathogen), 6 are extracellular (4 at depth 5 in *Costimulation by the CD28 family*'s tree, and the
  IL4/IL21 rows).
- **Depth is shallow for most.** 16 of the 25 pathways stop at depth 1 (T equals A there). The deep trees: *Transcriptional
  Regulation by TP53* (28 pathways joined, depth 5), *PIP3 activates AKT signaling* (11, 3), *Intrinsic Pathway for
  Apoptosis* (10, 4), *Signaling by ERBB2*, *S Phase*, *RUNX2* (9 each, 4), *Mitotic G1-G1/S* (8, 4), *CD28 family
  costimulation* (7, 4).

So "all of Reactome" in practice is these nine trees; the other 83 scored pathways answer as they do alone. The ceiling
(cases with no route alone and a route through T) is read from the arm results (R3).

## R3: the first runs of A and T, and the defect they exposed (2026-11-04)

Control `s133_c0` (build `…ctx133.c0`, code 910c5a7) reproduces specs/132's control case for case on every axis
(experimental 796, curator held-out 18,481, korkut 11,062, decryptM 1,802 cases: 0 predictions and 0 routes differ), so
the rebuilt 275-pathway base changed nothing for the 92 scored pathways.

**Reporting gap.** The API passed solver diagnostics through a fixed list that lacked `context_bridge_order`; the first
runs' logs read "order unknown". The containers did run `DS_CONTEXT_BRIDGE_ORDER=last` (`server.env`), and the case
changes below confirm it. Fixed in 4800ab3.

**Q1 holds.** Against specs/132 L1 (same bridges into Cell Cycle Checkpoints, order `first`), arm A returns ATM's four
experimental cases to the correct answer (knockout DOWN, overexpression UP at both readouts). **Q2 holds**: the 91
curator and 34 korkut cases of *Interleukin-4 and Interleukin-13 signaling* that L1 moved are back.

**The defect.** A also moved five PTEN overexpression cases in *Transcriptional Regulation by TP53* away from L1, four to
wrong. Traced by solving the pathway under both orders (scratch Julia; rule edges compared): the only rule difference is
three recycled-twin edges (specs/059), from TP53's own bridged copies of p-AKT1/2/3 [cytosol] (R-HSA-198356, -202062,
-3009367) to the p-AKT copies imported from *PIP3 activates AKT signaling* that feed the bridges. With the bridges set
aside the imported copies look recycled, and the rule, which matches copies by entity alone, joined them to TP53's.
Restored, that is a loop PIP3 → TP53 → PIP3: a shared molecule joining two pathways, which specs/132's principle forbids,
and feedback into the context, which this spec's tree forbids. Registered as amendment 1 (spec.md) before the
first runs' scores were read; implemented in 4800ab3: a bundle records each node's pathway (`context_origin.csv`) and,
under `last`, the solver drops every rule edge between two pathways.

`bench/analysis/origin_check_133.jl` checks that the drop equals running the structural rules on each pathway's part
alone. It does, on all 25 bridged pathways of both bundles (0 differ). The edges the rules had drawn across pathways were
not few: Cell Cycle Checkpoints 140 (both bundles), TP53 transcription 46 under A and 116 under T, 0–9 elsewhere. Under
`first` the same identity joins can form from any unbridged starting copy, so specs/056 and specs/132 arms likely held
some too; not measured, since neither is a candidate.

First-run scores, kept as the record of the defect, not as results (rule+F against c0; gate as registered):

| arm | experimental fixed / broken | curator held-out | korkut fixed / broken | decryptM | pooled gain [95%] |
|---|---|---|---|---|---|
| A (defective) | 5 / 0 | 12 / 4 | 69 / 310 | 0 / 22 | +0.0115 [+0.0014, +0.0225] |
| T (defective) | 5 / 0 | 12 / 6 | 63 / 312 | 0 / 22 | +0.0097 [−0.0001, +0.0208] |

The corrected arms are `s133_A1` and `s133_T1` (bundles `.A1`, `.T1`, identical to `.A`, `.T` but for
`context_origin.csv`; code 3541e6a). The control is not re-run: c0 has no bridges and no origin table, so the amendment
cannot change it.

## R4: A1 and T1 (amendment 1), and the source pathways that were not themselves (2026-11-04)

Arms `s133_A1`, `s133_T1` (code 3541e6a), against c0. Logs: "order last", 3,908–7,490 cross-pathway rule edges dropped
over all solves. Rule + F (the decision):

| arm | experimental fixed / broken / wrong→wrong | curator held-out | korkut | decryptM | pooled gain [95%] |
|---|---|---|---|---|---|
| A1+F | 5 / 0 / 10 | 12 / 4 | 71 / 250 | 0 / 22 | +0.0141 [+0.0037, +0.0254] |
| T1+F | 5 / 0 / 10 | 12 / 6 | 71 / 250 | 0 / 22 | +0.0141 [+0.0037, +0.0254] |

Accuracy / macro-F1 (rule+F minus c0): experimental +0.0063 / +0.0064, curator +0.0003 / +0.0005, korkut −0.0162 /
+0.0008, decryptM −0.0122 / −0.0065. The gate passes on balanced accuracy while korkut's raw accuracy falls: by unique
(drug set, readout) korkut is 30 fixed / 115 broken. Q4 holds: 306 of the 355 korkut changes are AKT readouts, and 324
of the 337 that have a PIP3 case repeat *PIP3 activates AKT signaling*'s own answer, which is wrong there (PIP3 alone
calls 655 korkut readouts DOWN, 515 measured unchanged). T adds nothing over A on the gated sets.

**These are not results either.** Every experimental case T1+F changes is PTEN in *Transcriptional Regulation by TP53*,
and against specs/132 L1 the overexpression cases moved to wrong (BAX, NDRG1, BTG2 DOWN; BIRC5 UP; measured the
reverse) while knockout read DOWN too. Traced (scratch Julia: every rule edge compared between orders, then the route):

1. Between `first` and `last`+amendment 1 the only differences are cross-pathway identity edges present under `first`:
   RNA twins joining TP53's and PIP3's PTEN mRNA both ways, recycled twins of GDP, GTP, H+, O2, cardiolipin and PE
   between pathways, carrier and component twins from TP53 into *HDR* copies. specs/132's "PTEN 9 / 0" leaned on
   those.
2. With them gone, the bridge sources, PIP3's p-AKT1/2/3 [cytosol], sit at the 100× ceiling under PTEN knockout AND
   overexpression. PIP3 solved alone reads PTEN overexpression DOWN (0.17×) and knockout UP (83.5×), correctly. The PIP3
   in the bundle was not PIP3: the builder (specs/056's design) removes from a source pathway every reaction of the
   pathway it feeds. TP53 transcription holds PTEN's transcription and translation, so the imported PIP3 lost exactly
   the PTEN machinery specs/095 and 100 fixed. (`--no-prune` changes nothing here: `main_curated` always pruned.)

Amendment 2 (`--whole`) enters each source pathway whole. `bench/analysis/source_check_133.jl` pins three starting
copies of each source pathway in turn and compares every one of its nodes with the pathway solved alone: in `.A2` the
worst |log ratio| over 49 (bundle, source) pairs is 2.2e-4 (solver tolerance). The same check on `.A1` (reduced and
pruned): 13 of the 49 pairs differ by more than 1% on nodes both hold, the worst |log ratio| 16 (*Signaling by
Insulin receptor* into *Signaling by Type 1 IGF1R*'s bundle), 4.6 in four, 1.1 for PIP3 in TP53's. The TP53 trace under `.A2`: PTEN
knockout puts p-AKT at 50–100×, TP53 tetramers fall to 0.28–0.53×, BAX 0.44× (measured DOWN); overexpression puts p-AKT
near 0, TP53 tetramers 1.15–1.32×, BAX 1.21× (measured UP). The signal now runs the curated route, PTEN ⊣ PIP3 → AKT
(PIP3 pathway) → bridge → AKT → MDM2 ⊣ TP53 → BAX (TP53 pathway). Arms `s133_A2`, `s133_T2` (code 5869e98).
